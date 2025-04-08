#!/usr/bin/env python3

'''Convert the binary marker output of flac to VTK (vtp) files.
'''

from __future__ import print_function
import sys, os, glob
import numpy as np
import flac
from flac2vtk import vts_dataarray


# filtering markers to only those within the domain bounds (in km)
filtering = False
xmin = 300
xmax = 700
zmin = -50
zmax = 100


def filter_marker(x, z, age, phase, ID, a1, a2, ntriag):
    # bool * bool is element-wise logical AND
    ind = (xmin <= x) * (x <= xmax) * (zmin <= z) * (z <= zmax)
    x = x[ind]
    z = z[ind]
    age = age[ind]
    phase = phase[ind]
    ID = ID[ind]
    a1 = a1[ind]
    a2 = a2[ind]
    ntriag = ntriag[ind]
    return x, z, age, phase, ID, a1, a2, ntriag


def main(path, start=1, end=-1):

    # changing directory
    os.chdir(path)

    fl = flac.Flac()
    if end == -1:
        end = fl.nrec

    if start == -1:
        vtplist = sorted(glob.glob('flacmarker.*.vtp'))
        lastframe = int(vtplist[-1][11:-4]) if vtplist else 0
        start = lastframe + 1

    for i in range(start, end+1):
        x, z, age, phase, ID, a1, a2, ntriag = fl.read_markers(i)
        if filtering:
            x, z, age, phase, ID, a1, a2, ntriag = filter_marker(x, z, age, phase, ID, a1, a2, ntriag)
        nmarkers = len(x)

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
        fvtp.write('  </PointData>\n')

        # point coordinates

        # VTK requires vector field (velocity, coordinate) has 3 components.
        # Allocating a 3-vector tmp array for VTK data output.
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
<FieldData>
  <DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">
    {1}
  </DataArray>
  <DataArray type="Float32" Name="CYCLE" NumberOfTuples="1" format="ascii">
    {2}
  </DataArray>
</FieldData>
<Piece NumberOfPoints="{0}">
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
    doc = '''usage: flacmarker2vtk.py [-f xmin,xmax,zmin,zmax] path [frame_min [frame_max]]

Processing flac marker output to VTK format.

If frame_min is -1, start from the latest vtp file.
If frame_max is not given, processing to latest frames
If both frame_min and frame_max are not given, processing all frames

-f xmin,xmax,zmin,zmax: if provided, only output markers within the domain range (in km)
'''

    if len(sys.argv) < 2:
        print(doc)
        sys.exit(1)

    try:
        n = 0
        if sys.argv[1] == '-f':
            filtering = True
            xmin, xmax, zmin, zmax = (float(f) for f in sys.argv[2].split(','))
            n = 2

        path = sys.argv[n+1]

        start = 1
        end = -1
        if len(sys.argv) >= n+3:
            start = int(sys.argv[n+2])
            if len(sys.argv) >= n+4:
                end = int(sys.argv[n+3])
    except:
        print(doc, '\n')
        raise

    main(path, start, end)
