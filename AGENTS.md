# Guidelines for AI Agents (GeoFLAC Workspace)

## 1. Codebase Overview
* **Language**: Free-form Fortran 90/95/2008 (`.f90`).
* **Parallelization**: OpenMP (CPU multi-threading) and OpenACC (GPU acceleration).
* **Model**: 2D finite-difference geomechanical solver (`flac` binary).

---

## 2. Directory & Documentation Overview
* **`src/`**: Fortran solver source files (`.f90`) and standard `Makefile`.
* **`doc/`**: Documentation ([`quickstart.md`](doc/quickstart.md), [`codebase_structure.md`](doc/developer/codebase_structure.md), [`parallel_architecture.md`](doc/developer/parallel_architecture.md), [`troubleshooting.md`](doc/troubleshooting.md)).
* **`examples/`**: Tutorial setups (e.g. [`examples/tutorial13-subduction/tutorial.md`](examples/tutorial13-subduction/tutorial.md)).
* **`util/`**: Post-processing tools (`plot_flac.py`, `flac2vtk.py`, `flacmarker2vtk.py`, `flac.py`, see [`util/README.md`](util/README.md)).
* **`benchmark/`**: Single-threaded regression testing suite (`make cmp`).

---

## 3. Fortran Coding Standards
* **`implicit none`**: Mandatory across all routines. Explicitly type all variables, loop counters, and external function return values.
* **Explicit Intents**: Mandatory `intent(in|out|inout)` on all dummy arguments.
* **No `goto`**: Prohibited. Use `if-else`, named `do-end do` loops with `exit`/`cycle`, or `select case`.
* **Dynamic I/O**: Use `open(newunit=u, ...)` and `inquire(file=..., exist=...)`.
* **Preserve Comments**: Do not modify or remove existing code comments.

---

## 4. Build & Run Guidelines
* **Build**: Run `cd src/ && make clean && make omp=1` (`FC=gfortran`, `-O2` default; parallel builds `make -j4` supported).
* **Sandboxed Executions**: Simulation runs require elevated/unsandboxed execution (`BypassSandbox: true`) to access files and hardware threads.
  * Example: `src/flac collision.inp`

---

## 5. Model Setup & Verification
* **Mesh & Zones**: Element counts (`nx-1`, `nz-1`) must match `nelem_per_zone` sum to avoid `STOP 17/19` crashes.
* **Benchmarking**: Run `cd benchmark/ && make cmp` with `BypassSandbox: true`. Correct output requires exact `0.0 0.0` relative difference.
* **Git**: Do not commit automatically unless requested. Use descriptive topic branches for new models.
