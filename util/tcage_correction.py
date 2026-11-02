#!/usr/bin/env python3
"""
The Correction Utility: tcage_correction.py

This script takes two inputs:
1.  Sample List File: A text file containing a list of geological thermochronology samples.
2.  FLAC VTS File: A structured grid `.vts` output file representing the model at a specific time step (along with corresponding `.vtp` marker files in the same directory).

It automatically performs:
*   Valley Correction (Low Topo): Interpolates the 2D age field at the sample's actual coordinates (x_sample, z_sample) from the current time step's grid elements.
*   Peak Correction (High Topo): Identifies the marker closest to (x_sample, z_model) and tracks it back in time to find the exact frame t2 when its depth matches the peak's elevation offset d. It extracts the surface age at that frame and adds the elapsed travel time (t1 - t2).
*   NaN & Out-of-Bounds Filtering: Employs robust fallbacks to safely handle unreset ages (which are stored as NaNs or -1.0 in the VTS/VTP files) and ensure clean interpolation.

---------------------------------------------------------
Expected Format of the Sample List File:
---------------------------------------------------------
The sample list file can be comma-separated or whitespace-separated. Lines starting with '#' or ';' are treated as comments and ignored.

Columns:
1. TC_type  : The type of thermochronometer (e.g. AFT, AHe, ZFT, ZHe, Orthoclase, Biotite, Muscovite, Hornblende). Case-insensitive.
2. x_sample : Horizontal coordinate of the sample in km.
3. z_sample : Elevation of the sample in km (with positive being above sea level / datum, matching model).
4. TC_age   : Observed age of the sample in Myr` (see Note below).
5. note     : (Optional) A note or description label for the sample.

Example File Contents:
# TC_type  x_sample  z_sample  TC_age  note
AFT       60.0      7.31      6.5     peak_sample
AHe       40.0      -0.84     2.1     valley_sample
AHe       20.0      -8.14     0.5     flat_sample

---------------------------------------------------------
Expected Format of the Screen Output:
---------------------------------------------------------
The final output is printed to standard output in a clean, space-separated format with a header line beginning with '#'.

Columns:
1. TC_type             : The type of thermochronometer system (e.g. AFT, AHe).
2. x_sample            : The sample horizontal position in km.
3. z_sample            : The sample vertical elevation in km.
4. TC_age              : The observed sample age in Myr.
5. model_age_corrected : The corrected modeled thermochronology age in Myr (see Note below).
6. note                : The optional comment label copied from the sample input file.

Example Screen Output:
# TC_type  x_sample  z_sample  TC_age  model_age_corrected  note
AFT       60.00     7.31      6.50    3.45                 peak_sample
AHe       40.00     -0.84     2.10    -1.00                valley_sample

Note: NaN represents an unreset age (the marker/cell was never heated above Tc to reset its clock), while -1.0 represents an unclosed / open system (currently hotter than Tc).
"""
import sys
import os
import glob
import numpy as np
from scipy.interpolate import griddata

# Append this script's directory so we can import flac
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import flac

def get_vts_array_name(tc_type):
    t = tc_type.lower()
    if t in ['ahe', 'age_ahe', 'ahe_age']:
        return 'age_AHE'
    elif t in ['aft', 'age_aft', 'aft_age']:
        return 'age_AFT'
    elif t in ['zft', 'age_zft', 'zft_age']:
        return 'age_ZFT'
    elif t in ['zhe', 'age_zhe', 'zhe_age']:
        return 'age_ZHe'
    elif t in ['orthoclase', 'age_orthoclase', 'orthoclase_age']:
        return 'age_orthoclase'
    elif t in ['biotite', 'age_biotite', 'biotite_age']:
        return 'age_biotite'
    elif t in ['muscovite', 'age_muscovite', 'muscovite_age']:
        return 'age_muscovite'
    elif t in ['hornblende', 'age_hornblende', 'hornblende_age']:
        return 'age_hornblende'
    else:
        return f"age_{tc_type}"

def get_vtp_array_name(tc_type):
    t = tc_type.lower()
    if t in ['ahe', 'age_ahe']:
        return 'age_AHE'
    elif t in ['aft', 'age_aft']:
        return 'age_AFT'
    elif t in ['zft', 'age_zft']:
        return 'age_ZFT'
    elif t in ['zhe', 'age_zhe']:
        return 'age_ZHe'
    elif t in ['orthoclase', 'age_orthoclase']:
        return 'age_orthoclase'
    elif t in ['biotite', 'age_biotite']:
        return 'age_biotite'
    elif t in ['muscovite', 'age_muscovite']:
        return 'age_muscovite'
    elif t in ['hornblende', 'age_hornblende']:
        return 'age_hornblende'
    else:
        return f"age_{tc_type}"

# Safe sorted interpolation helper with NaN filtering
def safe_interp(x_target, x_arr, y_arr):
    valid = ~np.isnan(y_arr)
    x_val = x_arr[valid]
    y_val = y_arr[valid]
    if len(x_val) == 0:
        return np.nan
    sort_idx = np.argsort(x_val)
    return np.interp(x_target, x_val[sort_idx], y_val[sort_idx])

def valley_correction(x_s, z_s, xc, zc, age_field):
    """
    Low-topography (valley) age correction.
    
    Interpolates the 2D thermochronology age field at the sample's actual underground 
    coordinates (x_s, z_s) using a 2D unstructured mesh interpolation.
    
    Parameters:
    -----------
    x_s : float
        Horizontal coordinate of the sample (km).
    z_s : float
        Vertical elevation of the sample (km), which is below the swath topography.
    xc : numpy.ndarray of shape (nx-1, nz-1)
        2D grid of horizontal coordinates of element cell-centers (km).
    zc : numpy.ndarray of shape (nx-1, nz-1)
        2D grid of vertical coordinates of element cell-centers (km).
    age_field : numpy.ndarray of shape (nx-1, nz-1)
        2D array of interpolated thermochronology ages for the target system (Myr).
        May contain NaNs representing unreset cells.
        
    Returns:
    --------
    float
        The corrected modeled thermochronology age at (x_s, z_s) in Myr.
        Returns NaN if no valid age elements exist in the field.
    """
    points = np.column_stack((xc.flatten(), zc.flatten()))
    values = age_field.flatten()
    valid = ~np.isnan(values)
    if np.any(valid):
        age_corrected = griddata(points[valid], values[valid], (x_s, z_s), method='linear')
        if np.isnan(age_corrected):
            # Fallback to nearest neighbor
            age_corrected = griddata(points[valid], values[valid], (x_s, z_s), method='nearest')
    else:
        age_corrected = np.nan
    return age_corrected

def peak_correction(x_s, z_model, dz, target_frame_idx, fl, surf_row, vts_arr_name, age_model_uncorr):
    """
    High-topography (peak) age correction using backward marker tracking.
    
    Traces the history of the closest marker backwards in time to determine the time t2 
    when it was at a depth equal to the sample's elevation offset (dz) below the past 
    topography. The corrected age is the surface age at that past location plus the 
    elapsed time (t1 - t2).
    
    Parameters:
    -----------
    x_s : float
        Horizontal coordinate of the sample (km).
    z_model : float
        Swath topography elevation today at x_s (km).
    dz : float
        Elevation offset of the sample above today's swath topography (km).
    target_frame_idx : int
        0-based index of the current time step frame in the simulation.
    fl : flac.FlacFromVTK
        Unified parser instance for StructuredGrid and PolyData reader.
    surf_row : int
        Row index in the grid representing the top surface (usually 0 or -1).
    vts_arr_name : str
        Name of the target thermochronology age array in the VTS dataset (e.g., 'age_AFT').
    age_model_uncorr : float
        The uncorrected modeled age today at today's surface at x_s (Myr).
        Used as a fallback if marker tracking fails or early ages are unreset.
        
    Returns:
    --------
    age_corrected : float
        The corrected modeled thermochronology age for the peak in Myr.
    corr_type : str
        A descriptive label of the correction type, including the resolved t2 frame time,
        or indicating if a NaN fallback was applied.
    """
    # 1. Load marker data at t1
    x_markers, z_markers, _, _, marker_ids, _, _, _ = fl.read_markers(target_frame_idx)
    
    # 2. Find marker closest to (x_s, z_model) today
    dist = (x_markers - x_s)**2 + (z_markers - z_model)**2
    idx_closest = np.argmin(dist)
    target_id = marker_ids[idx_closest]
    
    # 3. Track backwards
    depths = []
    times = []
    x_positions = []
    valid_frame_idx = []
    
    # We loop back through VTS file frame numbers
    for f_idx in range(target_frame_idx, -1, -1):
        vtp_t2 = os.path.join(fl.wdir, 'flacmarker.%06d.vtp' % fl.frames[f_idx])
        vts_t2 = os.path.join(fl.wdir, 'flac.%06d.vts' % fl.frames[f_idx])
        if not os.path.exists(vtp_t2) or not os.path.exists(vts_t2):
            continue
            
        xm_t2, zm_t2, _, _, ids_t2, _, _, _ = fl.read_markers(f_idx)
        
        idx_t2_arr = np.where(ids_t2 == target_id)[0]
        if len(idx_t2_arr) == 0:
            continue
        idx_t2 = idx_t2_arr[0]
        xm = xm_t2[idx_t2]
        zm = zm_t2[idx_t2]
        
        # Topography at t2
        x_mesh_t2, z_mesh_t2 = fl.read_mesh(f_idx)
        z_surf_m = safe_interp(xm, x_mesh_t2[:, surf_row], z_mesh_t2[:, surf_row])
        
        depths.append(z_surf_m - zm)
        
        # Populate and read time
        fl._get_vtk_data(f_idx)
        times.append(fl.time[f_idx])
        x_positions.append(xm)
        valid_frame_idx.append(f_idx)
    
    # Interpolate frame when marker depth was d = dz
    if len(depths) >= 2:
        depths = np.array(depths)
        times = np.array(times)
        x_positions = np.array(x_positions)
        valid_frame_idx = np.array(valid_frame_idx)
        
        sort_idx = np.argsort(depths)
        depths_s = depths[sort_idx]
        times_s = times[sort_idx]
        x_pos_s = x_positions[sort_idx]
        v_frame_idx_s = valid_frame_idx[sort_idx]
        
        t1_time = fl.time[target_frame_idx]
        
        if dz <= depths_s[0]:
            age_corrected = age_model_uncorr
            t2_time = times_s[0]
        elif dz >= depths_s[-1]:
            # Max available depth limit
            t2_time = times_s[-1]
            t2_x = x_pos_s[-1]
            t2_frame_idx = v_frame_idx_s[-1]
            
            age_field_t2 = fl.read_cell_data(t2_frame_idx, vts_arr_name)
            x_mesh_t2, _ = fl.read_mesh(t2_frame_idx)
            xc_t2 = 0.25 * (x_mesh_t2[:-1, :-1] + x_mesh_t2[1:, :-1] + x_mesh_t2[:-1, 1:] + x_mesh_t2[1:, 1:])
            age_surf_t2 = age_field_t2[:, surf_row]
            age_t2 = safe_interp(t2_x, xc_t2[:, surf_row], age_surf_t2)
            
            age_corrected = age_t2 + (t1_time - t2_time)
        else:
            t2_time = np.interp(dz, depths_s, times_s)
            t2_x = np.interp(dz, depths_s, x_pos_s)
            
            idx_above = np.searchsorted(depths_s, dz)
            idx_below = idx_above - 1
            
            frame_idx_a = v_frame_idx_s[idx_below]
            frame_idx_b = v_frame_idx_s[idx_above]
            
            # Fetch age at frame_idx_a
            age_field_a = fl.read_cell_data(frame_idx_a, vts_arr_name)
            x_mesh_a, _ = fl.read_mesh(frame_idx_a)
            xc_a = 0.25 * (x_mesh_a[:-1, :-1] + x_mesh_a[1:, :-1] + x_mesh_a[:-1, 1:] + x_mesh_a[1:, 1:])
            age_a = safe_interp(t2_x, xc_a[:, surf_row], age_field_a[:, surf_row])
            
            # Fetch age at frame_idx_b
            age_field_b = fl.read_cell_data(frame_idx_b, vts_arr_name)
            x_mesh_b, _ = fl.read_mesh(frame_idx_b)
            xc_b = 0.25 * (x_mesh_b[:-1, :-1] + x_mesh_b[1:, :-1] + x_mesh_b[:-1, 1:] + x_mesh_b[1:, 1:])
            age_b = safe_interp(t2_x, xc_b[:, surf_row], age_field_b[:, surf_row])
            
            time_a = times_s[idx_below]
            time_b = times_s[idx_above]
            
            age_t2 = age_a + (age_b - age_a) * (t2_time - time_a) / (time_b - time_a)
            age_corrected = age_t2 + (t1_time - t2_time)
            
        if np.isnan(age_corrected):
            age_corrected = age_model_uncorr
            corr_type = f"Peak (t2={t2_time:.2f}Ma) - NaN Fallback"
        else:
            corr_type = f"Peak (t2={t2_time:.2f}Ma)"
    else:
        age_corrected = age_model_uncorr
        corr_type = "Failed (No history)"
        
    return age_corrected, corr_type

def process_sample(s, x_surf, z_surf, fl, target_frame_idx, surf_row, xc_surf, xc, zc, target_frame_num, t1_time, z_mesh=None):
    """
    Processes a single geological sample by calculating topography offset and applying age corrections.
    
    Computes today's model topography at the sample's horizontal position, determines if the 
    sample is situated in a valley (below swath topo) or on a peak (above swath topo), applies 
    the appropriate correction algorithm (2D linear interpolation for valleys; marker tracking 
    for peaks), prints the formatted correction summary to standard output, and returns the result dictionary.
    
    Parameters:
    -----------
    s : dict
        A dictionary containing sample parameters:
        - 'type': str, thermochronometer system (e.g. 'AFT', 'AHe').
        - 'x': float, horizontal position of the sample (km).
        - 'z': float, vertical elevation of the sample (km).
        - 'age': float, observed sample age (Myr).
        - 'note': str, comment label.
    x_surf : numpy.ndarray of shape (nx,)
        Horizontal coordinates of top surface nodes (km).
    z_surf : numpy.ndarray of shape (nx,)
        Vertical elevations of top surface nodes (km).
    fl : flac.FlacFromVTK
        Initialized reader instance.
    target_frame_idx : int
        Index of current target VTS frame.
    surf_row : int
        Row index corresponding to the surface (0 or -1).
    xc_surf : numpy.ndarray of shape (nx-1,)
        Horizontal cell center coordinates of the top surface layer.
    xc : numpy.ndarray of shape (nx-1, nz-1)
        Cell center horizontal coordinates.
    zc : numpy.ndarray of shape (nx-1, nz-1)
        Cell center vertical coordinates.
    target_frame_num : int
        Frame number of the target frame.
    t1_time : float
        Model time of the target frame (Myr).
    z_mesh : numpy.ndarray of shape (nx, nz), optional
        Grid node vertical coordinates (km) used to calculate local element thickness.

    Returns:
    --------
    dict or None
        Dictionary containing corrected sample data, or None if reading failed.
    """
    tc_type = s['type']
    x_s = s['x']
    z_s = s['z']
    age_obs = s['age']
    
    # Calculate model topography at sample's x coordinate
    z_model = safe_interp(x_s, x_surf, z_surf)
    dz = z_s - z_model # positive for peak, negative for valley
    
    # Calculate local surface element thickness threshold (10% of local cell height, min 10m)
    if z_mesh is not None:
        sub_row = 1 if surf_row == 0 else -2
        cell_heights = np.abs(z_mesh[:, surf_row] - z_mesh[:, sub_row])
        dz_cell = safe_interp(x_s, x_surf, cell_heights)
        dz_thresh = max(0.01, 0.1 * dz_cell)
    else:
        dz_thresh = 0.05

    vts_arr_name = get_vts_array_name(tc_type)
    try:
        age_field = fl.read_cell_data(target_frame_idx, vts_arr_name)
    except Exception as err:
        print(f"Warning: Could not read array '{vts_arr_name}' for type '{tc_type}': {err}")
        return None

    # Get surface age profile for uncorrected modeled age
    age_surf = age_field[:, surf_row]
    age_model_uncorr = safe_interp(x_s, xc_surf, age_surf)
    
    age_corrected = age_model_uncorr
    corr_type = "None"

    # Apply correction logic using resolution-adaptive threshold
    if dz < -dz_thresh:
        age_corrected = valley_correction(x_s, z_s, xc, zc, age_field)
        corr_type = "Valley (2D)"
    elif dz > dz_thresh:
        try:
            age_corrected, corr_type = peak_correction(
                x_s, z_model, dz, target_frame_idx, fl,
                surf_row, vts_arr_name, age_model_uncorr
            )
        except Exception as e:
            print(f"Warning: Marker tracking failed for sample at x={x_s}: {e}")
            age_corrected = age_model_uncorr
            corr_type = "Failed (No marker)"
            
    print(f"{tc_type:<6} | {x_s:<6.1f} | {z_s:<6.2f} | {z_model:<8.2f} | {dz:<+6.2f} | {age_obs:<8.1f} | {age_model_uncorr:<9.1f} | {age_corrected:<8.1f} | {corr_type:<10}")
    
    return {
        'type': tc_type,
        'x': x_s,
        'z': z_s,
        'age_obs': age_obs,
        'age_uncorr': age_model_uncorr,
        'age_corr': age_corrected,
        'corr_type': corr_type
    }

def load_samples(filename):
    samples = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith(';'):
                continue
            if ',' in line:
                parts = [p.strip() for p in line.split(',')]
            else:
                parts = line.split()
            if len(parts) >= 4:
                tc_type = parts[0]
                x_sample = float(parts[1])
                z_sample = float(parts[2])
                age_sample = float(parts[3])
                note = parts[4] if len(parts) > 4 else ""
                samples.append({
                    'type': tc_type,
                    'x': x_sample,
                    'z': z_sample,
                    'age': age_sample,
                    'note': note
                })
    return samples

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 plot_corrected_tc.py <sample_file> <vts_file>")
        sys.exit(1)

    sample_file = sys.argv[1]
    vts_file = sys.argv[2]

    if not os.path.exists(sample_file):
        print(f"Error: Sample file '{sample_file}' not found.")
        sys.exit(1)
    if not os.path.exists(vts_file):
        print(f"Error: VTS file '{vts_file}' not found.")
        sys.exit(1)

    # 1. Load samples
    samples = load_samples(sample_file)
    print(f"Loaded {len(samples)} samples from '{sample_file}'.")

    vts_dir = os.path.dirname(os.path.abspath(vts_file))
    if not vts_dir:
        vts_dir = '.'
    vts_filename = os.path.basename(vts_file)
    
    # Initialize FlacFromVTK
    fl = flac.FlacFromVTK(wdir=vts_dir)
    
    # Identify target frame index
    target_frame_idx = -1
    for idx, fn in enumerate(fl.frames):
        fn_str = 'flac.%06d.vts' % fn
        if fn_str == vts_filename:
            target_frame_idx = idx
            break
            
    if target_frame_idx == -1:
        print(f"Error: Could not resolve frame index for VTS file '{vts_filename}'")
        sys.exit(1)

    target_frame_num = fl.frames[target_frame_idx]
    fl._get_vtk_data(target_frame_idx)
    t1_time = fl.time[target_frame_idx]
    print(f"Target frame number: {target_frame_num} (index: {target_frame_idx}) representing time: {t1_time:.2f} Myr")

    # Read current topography and grid
    x_mesh, z_mesh = fl.read_mesh(target_frame_idx)
    nx, nz = fl.nx, fl.nz
    if z_mesh[0, 0] > z_mesh[0, -1]:
        surf_row = 0
    else:
        surf_row = -1
    x_surf = x_mesh[:, surf_row]
    z_surf = z_mesh[:, surf_row]

    # Element-centered coordinates for cell data
    xc = 0.25 * (x_mesh[:-1, :-1] + x_mesh[1:, :-1] + x_mesh[:-1, 1:] + x_mesh[1:, 1:])
    zc = 0.25 * (z_mesh[:-1, :-1] + z_mesh[1:, :-1] + z_mesh[:-1, 1:] + z_mesh[1:, 1:])
    xc_surf = xc[:, surf_row]

    corrected_samples = []

    print("\nProcessing and correcting samples:")
    print(f"{'Type':<6} | {'x_s':<6} | {'z_s':<6} | {'z_model':<8} | {'dz':<6} | {'Age_obs':<8} | {'Age_model':<9} | {'Age_corr':<8} | {'Correction':<10}")
    print("-" * 90)

    for s in samples:
        res = process_sample(
            s, x_surf, z_surf, fl, target_frame_idx,
            surf_row, xc_surf, xc, zc, target_frame_num, t1_time, z_mesh=z_mesh
        )
        if res is not None:
            corrected_samples.append(res)

    print("\n# TC_type  x_sample  z_sample  TC_age  model_age_corrected  note")
    for cs, s in zip(corrected_samples, samples):
        note = s.get('note', '')
        print(f"{cs['type']:<8}  {cs['x']:<8.2f}  {cs['z']:<8.2f}  {cs['age_obs']:<6.2f}  {cs['age_corr']:<19.2f}  {note}")

if __name__ == '__main__':
    main()
