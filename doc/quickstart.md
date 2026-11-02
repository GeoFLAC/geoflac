# GeoFLAC Quick-Start Guide

This guide will walk you through compiling the solver, running the standard subduction case, and visualizing the numerical results.

---

## 1. Directory Layout

The GeoFLAC repository is structured as follows:

*   **`src/`**: The core Fortran 90 solver source code, including headers, math libraries, and the compilation Makefile.
*   **`util/`**: Post-processing and analysis utility scripts (e.g., VTK converters, plotting scripts).
*   **`examples/`**: Step-by-step benchmark tutorials (directories starting with `tutorial*`) and input templates.
*   **`doc/`**: Documentation split into user references (root level) and internal code guides ([`doc/developer/`](./developer/)).

---

## 2. Compilation

1. Navigate to the source directory:
   ```bash
   cd src
   ```
2. Build the optimized parallel binary with your Fortran compiler (e.g., `gfortran`):
   ```bash
   make clean
   make FC=gfortran omp=1
   ```
   (Note: Both `FC=gfortran` and `F90=gfortran` are supported.)
   This will generate a compiled binary named `flac` in the `src/` directory.

---

## 3. Running Your First Simulation

We provide a standard **Ocean-Ocean Subduction** model under the `examples/` directory.

1. Navigate to the examples folder:
   ```bash
   cd ../examples
   ```
2. Clear out any old output files:
   ```bash
   rm -f *.0 sys.msg output.asc _contents.rs
   ```
3. Set the number of threads for parallel execution (OpenMP):
   ```bash
   export OMP_NUM_THREADS=15
   ```
4. Run the subduction simulation using the input file:
   ```bash
   ../src/flac subduction.inp
   ```

You will see progress messages indicating step count and elapsed execution time. Once the simulation completes, it will display `Congratulations !`.

---

## 4. Visualizing Results with VisIt or ParaView

GeoFLAC writes its outputs in custom format files. We provide a utility script to convert these outputs into standard VTK files.

1. Convert the generated binary outputs to VTK files:
   ```bash
   ../util/flac2vtk.py .
   ```
   This will generate a series of `.vts` structured grid files (e.g., `flac.000001.vts`, `flac.000002.vts`, etc.).
2. Launch **VisIt** or **ParaView**.
3. Open the `flac*.vts` file series.
4. Select the variables you want to display:
    * **`Phase`**: Visualizes different material phases.
    * **`Strain rate`**: Visualizes the deformation. 
    * **`Temperature`**: Visualizes the temperature field. It is usually plotted as contours (isotherms).
    * **`Velocity`**: Visualizes velocity vectors.
5. Remember to fix the min/max range of the colorbar and the magnitude scale of vectors.

