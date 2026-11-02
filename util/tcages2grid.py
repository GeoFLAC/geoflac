#!/usr/bin/env python3

'''Interpolate marker thermochron ages onto grid elements and save them in grid VTS files.

Thermochronology age conventions:
  - NaN: Unreset system (the marker/cell was never heated above its closure temperature Tc to reset its clock).
  - -1.0: Unclosed / open system (the marker/cell was reset, but is currently hotter than Tc).
  - >= 0.0: Closed system thermochronology age in Myr (time elapsed since cooling below Tc).
'''

from __future__ import print_function
import sys
import os
import glob
import xml.etree.ElementTree as ET
import numpy as np
import zlib
import base64
import io
from scipy.interpolate import griddata
import pandas as pd
import flac
import flac2vtk

def read_data_array(el):
    dtype_str = el.get("type", "Float32")
    if dtype_str == "Float32":
        dtype = np.float32
    elif dtype_str == "Int32":
        dtype = np.int32
    else:
        dtype = np.float32 # fallback
        
    fmt = el.get("format", "ascii")
    
    text = el.text
    if text is None:
        text = ""
        
    if fmt == "ascii":
        return np.fromstring(text, sep=' ', dtype=dtype)
    elif fmt == "binary":
        text = text.strip()
        if not text:
            return np.array([], dtype=dtype)
        uncompressed_data = flac.FlacFromVTK._unpack_vtk(None, text)
        return np.frombuffer(uncompressed_data, dtype=dtype).copy()
    else:
        raise ValueError(f"Unknown format: {fmt}")


def make_vts_dataarray_element(data, name, format_type):
    orig_output_in_binary = flac2vtk.output_in_binary
    flac2vtk.output_in_binary = (format_type == "binary")
    
    f_str = io.StringIO()
    flac2vtk.vts_dataarray(f_str, data, name)
    
    flac2vtk.output_in_binary = orig_output_in_binary
    
    xml_fragment = f_str.getvalue()
    return ET.fromstring(xml_fragment)


def interpolate_vtp_to_vts(vtp_path, vts_path):
    print(f"Interpolating {os.path.basename(vtp_path)} -> {os.path.basename(vts_path)}   ", end='\r')
    sys.stdout.flush()
    
    # 1. Parse VTS to get grid coordinates and dimensions
    vts_tree = ET.parse(vts_path)
    vts_root = vts_tree.getroot()
    
    grid_el = vts_root.find(".//StructuredGrid")
    extent_str = grid_el.get("WholeExtent")
    extent_parts = [int(x) for x in extent_str.split()]
    nex = extent_parts[1]
    nez = extent_parts[3]
    nx = nex + 1
    nz = nez + 1
    
    # Get grid coordinates from VTS Points
    points_el = vts_root.find(".//Points/DataArray")
    points_data = read_data_array(points_el)
    points_data = points_data.reshape((nz, nx, 3))
    grid_x = points_data[:, :, 0]
    grid_z = points_data[:, :, 1]
    
    # Calculate cell centers
    xc = 0.25 * (grid_x[:-1, :-1] + grid_x[1:, :-1] + grid_x[:-1, 1:] + grid_x[1:, 1:])
    zc = 0.25 * (grid_z[:-1, :-1] + grid_z[1:, :-1] + grid_z[:-1, 1:] + grid_z[1:, 1:])
    
    # Check format of VTS data arrays
    first_da = vts_root.find(".//PointData/DataArray")
    if first_da is not None:
        format_type = first_da.get("format", "ascii")
    else:
        format_type = "binary" # fallback
        
    # 2. Parse VTP to get marker coordinates and thermochron ages
    vtp_tree = ET.parse(vtp_path)
    vtp_root = vtp_tree.getroot()
    
    piece_el = vtp_root.find(".//Piece")
    nmarkers = int(piece_el.get("NumberOfPoints", 0))
    
    if nmarkers == 0:
        # If there are no markers, we write all NaNs for any system found in PointData
        point_data_el = vtp_root.find(".//PointData")
        for da in point_data_el.findall("DataArray"):
            name = da.get("Name")
            if name and (name.startswith("age_") or name == "max_temp"):
                elem_age = np.full((nez, nex), np.nan, dtype=np.float32)
                new_da = make_vts_dataarray_element(elem_age, name, format_type)
                
                cell_data_el = vts_root.find(".//CellData")
                found = False
                for cda in cell_data_el.findall("DataArray"):
                    if cda.get("Name") == name:
                        idx = list(cell_data_el).index(cda)
                        cell_data_el[idx] = new_da
                        found = True
                        break
                if not found:
                    cell_data_el.append(new_da)
                    
        # Write VTS back
        xml_str = ET.tostring(vts_root, encoding='utf-8').decode('utf-8')
        with open(vts_path, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write(xml_str)
        return

    # Extract marker coordinates
    vtp_points_el = vtp_root.find(".//Points/DataArray")
    vtp_points_data = read_data_array(vtp_points_el).reshape((nmarkers, 3))
    x_markers = vtp_points_data[:, 0]
    z_markers = vtp_points_data[:, 1]
    
    # Extract ntriag
    point_data_el = vtp_root.find(".//PointData")
    ntriag_el = None
    for da in point_data_el.findall("DataArray"):
        if da.get("Name") == "ntriag":
            ntriag_el = da
            break
            
    if ntriag_el is None:
        raise ValueError("ntriag array not found in VTP PointData!")
        
    ntriag_markers = read_data_array(ntriag_el)
    
    # Compute element index for each marker
    ntriag_clean = np.maximum(1, ntriag_markers)
    k_m = (ntriag_clean - 1) % 2 + 1
    j_el = np.clip(((ntriag_clean - k_m) // 2) % (nz - 1), 0, nz - 2)
    i_el = np.clip(((ntriag_clean - k_m) // 2) // (nz - 1), 0, nx - 2)
    
    points = np.column_stack((x_markers, z_markers))
    
    # Find all thermochronology systems (starting with age_ or max_temp)
    cell_data_el = vts_root.find(".//CellData")
    for da in point_data_el.findall("DataArray"):
        name = da.get("Name")
        if name and (name.startswith("age_") or name == "max_temp"):
            ages = read_data_array(da)
            
            # Filter out young outliers from newly inserted/replenished markers.
            # Newly inserted markers are assigned large IDs (> max ID of initial markers).
            # If they are created in a cold element, they instantly record a cooling age of 0.0 Myr.
            if name.startswith("age_"):
                # Try to find the maximum initial marker ID from the first frame
                max_initial_id = None
                dirname = os.path.dirname(vtp_path)
                initial_vtp_path = os.path.join(dirname, "flacmarker.000001.vtp")
                if os.path.exists(initial_vtp_path):
                    try:
                        initial_tree = ET.parse(initial_vtp_path)
                        initial_root = initial_tree.getroot()
                        initial_point_data_el = initial_root.find(".//PointData")
                        if initial_point_data_el is not None:
                            for tmp_da in initial_point_data_el.findall("DataArray"):
                                if tmp_da.get("Name") == "ID":
                                    initial_ids = read_data_array(tmp_da)
                                    if len(initial_ids) > 0:
                                        max_initial_id = np.max(initial_ids)
                                    break
                    except Exception:
                        pass
                
                # Try to read the ID array of the current frame
                current_ids = None
                for tmp_da in point_data_el.findall("DataArray"):
                    if tmp_da.get("Name") == "ID":
                        current_ids = read_data_array(tmp_da)
                        break
                
                if max_initial_id is not None and current_ids is not None and len(current_ids) == len(ages):
                    # Precise filter: exclude any markers with IDs larger than the initial marker count
                    outlier_mask = (current_ids > max_initial_id)
                    if np.any(outlier_mask):
                        ages = ages.copy()
                        ages[outlier_mask] = np.nan
                else:
                    # Robust fallback: filter outliers using a local median absolute/relative threshold
                    elem_markers = {}
                    valid_indices = np.where(ages > -0.5)[0]
                    for idx in valid_indices:
                        key = (j_el[idx], i_el[idx])
                        if key not in elem_markers:
                            elem_markers[key] = []
                        elem_markers[key].append(idx)
                    
                    outlier_indices = []
                    for key, indices in elem_markers.items():
                        if len(indices) >= 3:
                            element_ages = ages[indices]
                            median_age = np.median(element_ages)
                            if median_age > 0.1:
                                for idx in indices:
                                    if ages[idx] < 0.15 * median_age:
                                        outlier_indices.append(idx)
                                        
                    if len(outlier_indices) > 0:
                        ages = ages.copy()
                        ages[outlier_indices] = np.nan
            
            # Determine state of each cell that contains markers
            # Priority: Reset (0) > Unclosed (1) > Unreset (2)
            has_reset = np.zeros((nez, nex), dtype=bool)
            has_unclosed = np.zeros((nez, nex), dtype=bool)
            has_unreset = np.zeros((nez, nex), dtype=bool)
            
            has_reset[j_el[ages > -0.5], i_el[ages > -0.5]] = True
            has_unclosed[j_el[ages == -1.0], i_el[ages == -1.0]] = True
            has_unreset[j_el[np.isnan(ages)], i_el[np.isnan(ages)]] = True
            
            elem_has_marker = has_reset | has_unclosed | has_unreset
            
            cell_state = np.zeros((nez, nex), dtype=np.int32)
            cell_state[has_unreset] = 2
            cell_state[has_unclosed] = 1
            cell_state[has_reset] = 0
            
            # Interpolate state for cells that have no markers using nearest neighbor
            points_with_marker = np.column_stack((xc[elem_has_marker], zc[elem_has_marker]))
            states_with_marker = cell_state[elem_has_marker]
            elem_state = griddata(points_with_marker, states_with_marker, (xc, zc), method='nearest')
            
            # Initialize elem_age to NaN
            elem_age = np.full((nez, nex), np.nan, dtype=np.float32)
            
            # For elements that are unclosed (State 1), set to -1.0
            elem_age[elem_state == 1] = -1.0
            
            # For elements that are reset (State 0), interpolate the actual ages
            reset_mask = (elem_state == 0)
            if np.any(reset_mask):
                ages_reset = ages[ages > -0.5]
                points_reset = points[ages > -0.5]
                if len(ages_reset) > 0:
                    ages_interp = griddata(points_reset, ages_reset, (xc, zc), method='linear')
                    
                    # Fill any NaNs in the linear interpolation (e.g. outside convex hull)
                    # using nearest neighbor interpolation of the ages
                    nan_in_interp = np.isnan(ages_interp)
                    if np.any(nan_in_interp):
                        ages_nearest = griddata(points_reset, ages_reset, (xc, zc), method='nearest')
                        ages_interp[nan_in_interp] = ages_nearest[nan_in_interp]
                    
                    elem_age[reset_mask] = ages_interp[reset_mask]
                else:
                    elem_age[reset_mask] = -1.0
                    
            # Mask out elements if they do not contain any markers
            elem_age[~elem_has_marker] = np.nan
                
            new_da = make_vts_dataarray_element(elem_age.astype(np.float32), name, format_type)
            
            # Update/Insert into CellData of VTS
            found = False
            for cda in cell_data_el.findall("DataArray"):
                if cda.get("Name") == name:
                    idx = list(cell_data_el).index(cda)
                    cell_data_el[idx] = new_da
                    found = True
                    break
            if not found:
                cell_data_el.append(new_da)
                
    # Write VTS back
    xml_str = ET.tostring(vts_root, encoding='utf-8').decode('utf-8')
    with open(vts_path, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write(xml_str)


def main(args):
    doc = '''usage: tcages2grid.py path [frame_min [frame_max]]
       tcages2grid.py vtp_file vts_file

Interpolate marker thermochron ages onto grid elements and save them in grid VTS files.
'''
    if len(args) < 1:
        print(doc)
        sys.exit(1)
        
    path_or_file = args[0]
    
    if path_or_file.endswith('.vtp'):
        # Mode: Single VTP / VTS pair
        vtp_file = path_or_file
        if len(args) >= 2:
            vts_file = args[1]
        else:
            # Infer VTS file name
            # Replace flacmarker with flac and .vtp with .vts
            base = os.path.basename(vtp_file)
            dirname = os.path.dirname(vtp_file)
            if base.startswith("flacmarker."):
                vts_base = base.replace("flacmarker.", "flac.", 1).replace(".vtp", ".vts", 1)
            else:
                vts_base = base.replace(".vtp", ".vts", 1)
            vts_file = os.path.join(dirname, vts_base)
            
        if not os.path.exists(vtp_file):
            print(f"Error: {vtp_file} not found!")
            sys.exit(1)
        if not os.path.exists(vts_file):
            print(f"Error: {vts_file} not found!")
            sys.exit(1)
            
        interpolate_vtp_to_vts(vtp_file, vts_file)
        print()
        
    else:
        # Mode: Directory with start/end frames
        path = path_or_file
        start = 1
        end = -1
        if len(args) >= 2:
            start = int(args[1])
            if len(args) >= 3:
                end = int(args[2])
                
        # Find matching VTP files
        vtplist = sorted(glob.glob(os.path.join(path, 'flacmarker.*.vtp')))
        if not vtplist:
            print(f"Error: No flacmarker.*.vtp files found in {path}")
            sys.exit(1)
            
        if end == -1:
            end = int(os.path.basename(vtplist[-1])[11:-4])
            
        for i in range(start, end + 1):
            vtp_file = os.path.join(path, f'flacmarker.{i:06d}.vtp')
            vts_file = os.path.join(path, f'flac.{i:06d}.vts')
            if os.path.exists(vtp_file) and os.path.exists(vts_file):
                interpolate_vtp_to_vts(vtp_file, vts_file)
            else:
                print(f"Skipping frame {i}: VTP or VTS file not found.")
        print()

if __name__ == '__main__':
    main(sys.argv[1:])
