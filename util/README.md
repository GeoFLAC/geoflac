# Geoflac Utility Scripts (`util/`)

This directory contains Python libraries, post-processing tools, visualization scripts, and helper utilities for the **geoflac** model.

## Core Libraries

* **[flac.py](flac.py)**: The central library providing Python interface classes for geoflac outputs:
  * `Flac`: Reads binary model output files (e.g., meshes, temperatures, densities, velocities, phases).
  * `FlacFromVTK`: Parses StructuredGrid (`.vts`) output files.
  * Helper functions such as `make_uniform_grid()`, `nearest_neighbor_interpolation2d()`, and `printing()`.
* **[flac_interpolate.py](flac_interpolate.py)**: Interpolates irregularly spaced node/element data from the model mesh onto a uniform grid. Mainly used for plotting.

## Plotting & Visualization

* **[plot_flac.py](plot_flac.py)**: A unified command-line plotting utility that consolidates older standalone plotting scripts.
  * **Usage**:
    * Matplotlib Topography Plot: `python3 plot_flac.py -m topo <frame>`
    * GMT Phase Plot: `python3 plot_flac.py -m phase <frame>`
    * GMT Strain Rate Plot: `python3 plot_flac.py -m strain_rate <frame>`
    * GMT Stress Plot: `python3 plot_flac.py -m stress <frame>`
    * GMT Stacked Multi-panel Plot: `python3 plot_flac.py -m all <frame>`
  * Relies on local Color Palette Table (`.cpt`) files: [phase15.cpt](phase15.cpt), [strain_rate.cpt](strain_rate.cpt), and [stress.cpt](stress.cpt).
* **[visc_profile.py](visc_profile.py)**: Extracts and plots vertical profiles of viscosity at specific x-coordinates.

## VTK Exporters & Post-processing

* **[flac2vtk.py](flac2vtk.py)**: Converts raw binary grid outputs (mesh, velocity, temperature, phase, stress, etc.) into VTK StructuredGrid (`.vts`) files.
* **[flacmarker2vtk.py](flacmarker2vtk.py)**: Converts raw binary Lagrangian marker files to VTK PolyData (`.vtp`) files.
  * Includes the `-t` option to compute thermochronological cooling/closure ages and track historical maximum temperature (`max_temp`) on markers offline.
* **[tcages2grid.py](tcages2grid.py)**: Reads marker thermochronological ages and historical maximum temperature (`max_temp`) from a `.vtp` file and maps/interpolates them onto the cell elements of the corresponding `.vts` grid file, resolving overlapping marker states (e.g. reset vs. unreset).

## Geophysical Calculations

* **[flac_gravity.py](flac_gravity.py)**: Computes surface topography and gravity anomalies from model density distributions.
  * Supports options to fill sedimentary basins to a custom depth (`-s`), erode forearc topography above sea level (`-e`), and choose the reference gravity frame (`-r`).
* **[seismic_velocities.py](seismic_velocities.py)**: Computes seismic velocities ($V_p$, $V_s$) based on density and stress fields.

## Model Configuration & Control

* **[update_inp.py](update_inp.py)**: A comprehensive script to inspect, format, modify, and migrate geoflac input parameter files (`.inp`).
* **[save2rs.py](save2rs.py)**: Converts regular model output frames to restart (`.rs`) files.
* **[rs-extract.py](rs-extract.py)**: Extracts data arrays from restart files.
* **[restart2init.py](restart2init.py)**: Modifies model state/restart metadata to reset the iteration/frame count to 0, useful for using a matured frame as a starting point.

## Helper & Developer Utilities

* **[check-acc-decl.py](check-acc-decl.py)**: Developer diagnostic utility that parses Fortran files to validate that OpenACC GPU declarations match the declared local variables.
* **[flacbackproj.py](flacbackproj.py)**: Back-projects coordinates of Lagrangian markers to map tracking points through time.
* **[readallmarkers.py](readallmarkers.py)**: Diagnostic reader to parse binary marker output files.

## Support Files

* **[phase15.cpt](phase15.cpt)**: GMT color palette table for lithological phases.
* **[strain_rate.cpt](strain_rate.cpt)**: GMT color palette table for strain rate.
* **[stress.cpt](stress.cpt)**: GMT color palette table for stress fields.
* **[thermo_chron.dat](thermo_chron.dat)**: Constants and parameter database for thermochronometer systems.

## Deprecated / Scratch Scripts

* **[marker2vtk.py](marker2vtk.py)**: Outdated marker-to-VTK converter, replaced by [flacmarker2vtk.py](flacmarker2vtk.py).
* **[intp_mark.py](intp_mark.py)**: Obsolete test script for marker interpolation.
* **[irtemp.py](irtemp.py)**: Legacy scratch file for testing temperature reads.
