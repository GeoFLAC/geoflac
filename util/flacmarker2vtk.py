#!/usr/bin/env python3

'''Convert the binary marker output of flac to VTK (vtp) files.

Thermochronology age conventions:
  - NaN: Unreset system (the marker was never heated above its closure temperature Tc to reset its clock).
  - -1.0: Unclosed / open system (the marker was reset, but is currently hotter than Tc).
  - >= 0.0: Closed system cooling age in Myr (time elapsed since cooling below Tc).
'''

from __future__ import print_function
import sys, os, glob
import numpy as np
import flac
from flac2vtk import vts_dataarray

max_thermochron_temp = 620.0

def parse_chron_file(filepath):
    systems = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    current_system = None
    rates = []
    temps = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith(';'):
            continue
        if line.startswith('0 (#'):
            if current_system is not None:
                systems.append((np.array(rates), np.array(temps)))
            current_system = line
            rates = []
            temps = []
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                rate = float(parts[0])
                temp = float(parts[1])
                rates.append(rate)
                temps.append(temp)
            except ValueError:
                pass
    if current_system is not None:
        systems.append((np.array(rates), np.array(temps)))
    return systems


def filter_marker(x, z, age, phase, ID, a1, a2, ntriag, xlim, zlim):
    xmin, xmax = xlim
    zmin, zmax = zlim
    ind = (xmin <= x) * (x <= xmax) * (zmin <= z) * (z <= zmax)
    x = x[ind]
    z = z[ind]
    age = age[ind]
    phase = phase[ind]
    ID = ID[ind]
    a1 = a1[ind]
    a2 = a2[ind]
    ntriag = ntriag[ind]
    return x, z, age, phase, ID, a1, a2, ntriag, ind


def main(path, start=1, end=-1, thermochron=False, filtering=False, xlim=None, zlim=None):

    # changing directory
    os.chdir(path)

    fl = flac.Flac()
    if end == -1:
        end = fl.nrec

    if start == -1:
        vtplist = sorted(glob.glob('flacmarker.*.vtp'))
        lastframe = int(vtplist[-1][11:-4]) if vtplist else 0
        start = lastframe + 1

    if thermochron:
        # Resolve path to the thermochronology database file
        script_dir = os.path.dirname(os.path.abspath(__file__))
        chron_file = os.path.join(script_dir, 'thermo_chron.dat')
        if not os.path.exists(chron_file):
            print(f"Error: thermochron reference file not found at {chron_file}")
            sys.exit(1)

        # Verify that node cooling rates are available
        if not os.path.exists('coolingrate.0'):
            print("Error: coolingrate.0 not found in the current directory.")
            sys.exit(1)

        # Parse reference cooling rates and closure temperatures for all systems
        systems = parse_chron_file(chron_file)
        nchron = len(systems)
        thermochron_system_names = ['ZFT', 'ZHe', 'AFT', 'AHE', 'orthoclase', 'biotite', 'muscovite', 'hornblende'][:nchron]

        # Initialize tracking arrays for all markers across the model's history.
        # We index these arrays using the unique marker ID.
        # Define max_markers as nx * nz * 32 (upper bound of marker ID values)
        max_markers = fl.nx * fl.nz * 32
        # chron_times stores the model time at which the marker cooled past its closure temperature.
        # Sentinel values: -1.0 = unclosed/open system, NaN = unreset (never heated enough to open).
        chron_times = [np.full(max_markers + 1, -1.0, dtype=np.float32) for _ in range(nchron)]
        # max_temp_history stores the maximum temperature (in °C) each marker has experienced in its lifetime.
        max_temp_history = np.full(max_markers + 1, -1.0, dtype=np.float32)

    # If computing thermochronology, we must start tracking from frame 1 to compile the correct history,
    # even if we only write output VTP files starting from the specified 'start' frame.
    loop_start = 1 if thermochron else start
    for i in range(loop_start, end+1):
        # Read coordinates and marker properties for the current frame
        x, z, age, phase, ID, a1, a2, ntriag = fl.read_markers(i)
        nmarkers = len(x)

        if thermochron:
            t_frame = fl.time[i - 1]

            # Read nodal temperature and cooling rate fields from binary outputs
            T_nodes = fl.read_temperature(i)
            CR_nodes = fl.read_coolingrate(i)

            if nmarkers > 0:
                # Interpolate nodal temperature and cooling rate onto the marker positions
                t_temps = flac.marker_interpolate_node(ntriag, a1, a2, fl.nz, T_nodes)
                t_rates = flac.marker_interpolate_node(ntriag, a1, a2, fl.nz, CR_nodes)

                m_ids = ID.astype(np.int32)

                # Update the historical maximum temperature for each active marker ID
                max_temp_history[m_ids] = np.maximum(max_temp_history[m_ids], t_temps)

                for c in range(nchron):
                    rates_ref, temps_ref = systems[c]

                    # Retrieve previous closure times for the active marker IDs
                    c_times = chron_times[c][m_ids]

                    # Thermal resetting: If a marker's temperature exceeds its cooling-rate-dependent
                    # closure temperature, reset its closure time to -1.0 (unclosed/open).
                    cooling_rate_val_all = np.maximum(1e-10, -t_rates)
                    t_closure_all = np.interp(cooling_rate_val_all, rates_ref, temps_ref)
                    hot_mask = (t_temps > t_closure_all)
                    c_times[hot_mask] = -1.0

                    # Evaluate closure criteria for markers that are currently unclosed (-1.0)
                    unclosed_mask = (c_times == -1.0)
                    if np.any(unclosed_mask):
                        unclosed_temps = t_temps[unclosed_mask]
                        unclosed_rates = t_rates[unclosed_mask]

                        # Ensure marker temperature is below the maximum allowed closure threshold
                        valid_temp = (unclosed_temps <= max_thermochron_temp)
                        cooling_rate_val = np.maximum(1e-10, -unclosed_rates)
                        t_closure = np.interp(cooling_rate_val, rates_ref, temps_ref)

                        # Marker closes if its temperature is <= closure temperature and is valid
                        should_close = (unclosed_temps <= t_closure) & valid_temp
                        if i == 1:
                            # Frame 1: If it starts below the closure temperature, it was never reset
                            # (initialized as NaN/unreset), otherwise it remains -1.0 (open/unclosed).
                            c_times[unclosed_mask] = np.where(should_close, np.nan, -1.0)
                        else:
                            # Subsequent frames: If it cools past the closure temp, record the closure time.
                            c_times[unclosed_mask] = np.where(should_close, t_frame, -1.0)

                    # Save updated closure times back to the global history array
                    chron_times[c][m_ids] = c_times

            # Prepare the point data arrays to be written if we have reached the output range
            if i >= start:
                if nmarkers > 0:
                    current_frame_ages = []
                    for c in range(nchron):
                        c_times = chron_times[c][m_ids]
                        ages = np.full_like(c_times, -1.0, dtype=np.float32)

                        # Reset markers: Age = current time - closure time
                        reset_mask = (c_times >= 0.0)
                        ages[reset_mask] = t_frame - c_times[reset_mask]

                        # Unreset markers: Age remains NaN
                        unreset_mask = np.isnan(c_times)
                        ages[unreset_mask] = np.nan

                        current_frame_ages.append(ages)
                    
                    # Extract the historical maximum temperatures for current markers
                    current_max_temp = max_temp_history[m_ids]
                else:
                    current_frame_ages = [np.array([], dtype=np.float32) for _ in range(nchron)]
                    current_max_temp = np.array([], dtype=np.float32)

        if i < start:
            continue

        if filtering:
            x, z, age, phase, ID, a1, a2, ntriag, ind = filter_marker(x, z, age, phase, ID, a1, a2, ntriag, xlim, zlim)
            nmarkers = len(x)
            if thermochron:
                current_frame_ages = [ages[ind] for ages in current_frame_ages]
                current_max_temp = current_max_temp[ind]

        print('Writing record #%d, model time=%.3e, %d markers' % (i, fl.time[i-1], nmarkers), end='\r')
        sys.stdout.flush()
        fvtp = open('flacmarker.%06d.vtp' % i, 'w')
        vtp_header(fvtp, nmarkers, fl.time[i-1], fl.steps[i-1])

        # point-based data
        fvtp.write('  <PointData>\n')
        vts_dataarray(fvtp, age, 'age', 1)
        vts_dataarray(fvtp, phase.astype(np.int32), 'phase', 1)
        vts_dataarray(fvtp, ID.astype(np.int32), 'ID', 1)
        vts_dataarray(fvtp, a1, 'a1', 1)
        vts_dataarray(fvtp, a2, 'a2', 1)
        vts_dataarray(fvtp, ntriag.astype(np.int32), 'ntriag', 1)

        if thermochron:
            for name, ages_array in zip(thermochron_system_names, current_frame_ages):
                vts_dataarray(fvtp, ages_array, f'age_{name}', 1)
            vts_dataarray(fvtp, current_max_temp, 'max_temp', 1)

        fvtp.write('  </PointData>\n')

        # point coordinates
        tmp = np.zeros((nmarkers, 3), dtype=x.dtype)
        tmp[:,0] = x
        tmp[:,1] = z
        fvtp.write('  <Points>\n')
        vts_dataarray(fvtp, tmp, '', 3)
        fvtp.write('  </Points>\n')

        vtp_footer(fvtp)
        fvtp.close()
    print()
    return


def vtp_header(f, npoints, time, step):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<PolyData>
<Piece NumberOfPoints="{0}">
<FieldData>
  <DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">
    {1}
  </DataArray>
  <DataArray type="Float32" Name="CYCLE" NumberOfTuples="1" format="ascii">
    {2}
  </DataArray>
</FieldData>
'''.format(npoints, time, step))
    return


def vtp_footer(f):
    f.write(
'''</Piece>
</PolyData>
</VTKFile>
''')
    return


if __name__ == '__main__':
    doc = '''usage: flacmarker2vtk.py [-t] [-f xmin,xmax,zmin,zmax] path [frame_min [frame_max]]

Processing flac marker output to VTK format.

If -t is specified, compute thermochronology closure ages of markers and save them in the vtp files.
If frame_min is -1, start from the latest vtp file.
If frame_max is not given, processing to latest frames
If both frame_min and frame_max are not given, processing all frames

-f xmin,xmax,zmin,zmax: if provided, only output markers within the domain range (in km)
'''

    if len(sys.argv) < 2:
        print(doc)
        sys.exit(1)

    args = sys.argv[1:]
    thermochron = False
    if '-t' in args:
        thermochron = True
        args.remove('-t')

    filtering = False
    xmin, xmax, zmin, zmax = 0.0, 0.0, 0.0, 0.0
    if '-f' in args:
        idx = args.index('-f')
        filtering = True
        xmin, xmax, zmin, zmax = (float(f) for f in args[idx+1].split(','))
        args.pop(idx+1)
        args.pop(idx)

    if len(args) < 1:
        print(doc)
        sys.exit(1)

    path = args[0]
    start = 1
    end = -1
    if len(args) >= 2:
        start = int(args[1])
        if len(args) >= 3:
            end = int(args[2])

    main(path, start, end, thermochron=thermochron, filtering=filtering, xlim=(xmin, xmax), zlim=(zmin, zmax))
