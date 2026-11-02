# DATA ARRAYS IN GEOFLAC

This document describes the primary data arrays used in the `geoflac` Fortran engine, including their shape, dimension meanings, and where they are updated.

---

## `cord`

*   **Meaning:** Eulerian coordinates of the grid nodes.
*   **Shape:** `(nz, nx, 2)`
*   **Dimension 1 (`nz`):** Vertical index (Z direction), from 1 to `nz`. Usually 1 is the surface and `nz` is the bottom.
*   **Dimension 2 (`nx`):** Horizontal index (X direction), from 1 to `nx`. From left to right.
*   **Dimension 3 (2):** Coordinate component:
    *   `1`: X-coordinate [m]
    *   `2`: Z-coordinate [m] (typically negative for depth)
*   **Update Locations:**
    *   **Initialization:** `src/init_cord.f90` (generates the initial grid or reads from `coordfile`).
    *   **Main Update:** `src/fl_move.f90` (updates coordinates using nodal velocities: `cord = cord + vel * dt`).
    *   **Topography:** `src/fl_move.f90` also adjusts the surface (`cord(1,:,2)`) for erosion, sedimentation, and extrusion.
    *   **Remeshing:** `src/rem_cord.f90` and `src/remesh.f90` recalculate coordinates during grid reorganization.

---

## `xoriginal`

*   **Meaning:** Nodal X-coordinates at the beginning of the experiment (the initial reference state). Used to calculate the cumulative horizontal displacement of each grid node since the start of the simulation (`cord(:,:,1) - xoriginal`) [m].
*   **Shape:** `(nz, nx)`
*   **Dimension 1 (`nz`):** Vertical index (node), from 1 (surface) to `nz` (bottom).
*   **Dimension 2 (`nx`):** Horizontal index (node), from 1 (left) to `nx` (right).
*   **Update Locations:**
    *   **Initialization:** `src/setflac.f90` (stores the initial grid coordinates: `xoriginal = cord(:,:,1)`).
    *   **Remeshing:** `src/remesh.f90` (interpolates these reference coordinates onto the new grid so that displacement tracking is preserved across remeshing events).
    *   **Restart/State:** Written in `src/saveflac.f90` and loaded in `src/rsflac.f90`.

---

## `zoriginal`

*   **Meaning:** Nodal Z-coordinates at the beginning of the experiment (the initial reference state). Used to calculate the cumulative vertical displacement of each grid node since the start of the simulation (`cord(:,:,2) - zoriginal`) [m].
*   **Shape:** `(nz, nx)`
*   **Dimension 1 (`nz`):** Vertical index (node), from 1 (surface) to `nz` (bottom).
*   **Dimension 2 (`nx`):** Horizontal index (node), from 1 (left) to `nx` (right).
*   **Update Locations:**
    *   **Initialization:** `src/setflac.f90` (stores the initial grid coordinates: `zoriginal = cord(:,:,2)`).
    *   **Remeshing:** `src/remesh.f90` (interpolates these reference coordinates onto the new grid so that displacement tracking is preserved across remeshing events).
    *   **Restart/State:** Written in `src/saveflac.f90` and loaded in `src/rsflac.f90`.

---

## `temp`

*   **Meaning:** Temperature at grid nodes [Celsius].
*   **Shape:** `(nz, nx)`
*   **Dimension 1 (`nz`):** Vertical index (Z direction), from 1 (surface) to `nz` (bottom).
*   **Dimension 2 (`nx`):** Horizontal index (X direction), from 1 (left) to `nx` (right).
*   **Update Locations:**
    *   **Initialization:** `src/init_temp.f90` (calculates initial geotherms or reads from `tempfile`).
    *   **Initial Perturbation:** `src/init_temp.f90` (adds Gaussian temperature anomalies defined by `inhom` parameters).
    *   **Main Update:** `src/fl_therm.f90` (solves the heat conduction equation).
    *   **Phase Change / Shear Heating:** Indirectly affects temperature via source terms in `fl_therm.f90`.
    *   **Remeshing:** Interpolated in `src/remesh.f90`.

---

## `temp0`

*   **Meaning:** Nodal temperature values at the start of the current time step [Celsius]. It serves as a reference temperature for thermal stress / expansion calculations.
*   **Shape:** `(nz, nx)`
*   **Dimension 1 (`nz`):** Vertical index (node), from 1 (surface) to `nz` (bottom).
*   **Dimension 2 (`nx`):** Horizontal index (node), from 1 (left) to `nx` (right).
*   **Update Locations:**
    *   **Initialization:** `src/setflac.f90` (set equal to `temp`).
    *   **Main Update:** `src/fl_therm.f90` (sets `temp0 = temp` at the very beginning of the thermal solver step, before calculating new temperatures).
    *   **Remeshing:** `src/remesh.f90` (updates `temp0` to match the new interpolated temperatures).
    *   **Restart:** `src/rsflac.f90` (restored from checkpoint files).
    *   **Usage:** Used in `src/fl_rheol.f90` to compute thermal stress increments: $\sigma_{\text{therm}} = -\alpha \cdot K_b \cdot (T - T_0)$ where $T_0$ is the element-averaged temperature from `temp0`.

---

## `flux`

*   **Meaning:** Thermal heat flux components computed inside the grid elements. Unlike the mechanical solver which divides quadrilateral elements into 4 sub-triangles, the heat conduction solver uses a simpler subdivision of 2 triangles per element.
*   **Shape:** `(2, 2, nz-1, nx-1)`
*   **Dimension 1 (2):** Triangle index within the quadrilateral element (1: Triangle A, 2: Triangle B).
*   **Dimension 2 (2):** Heat flux component direction:
    *   `1`: Horizontal heat flux $q_x$ [W/m²]
    *   `2`: Vertical heat flux $q_z$ [W/m²]
*   **Dimension 3 (`nz-1`):** Vertical index (element).
*   **Dimension 4 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Main Update:** `src/fl_therm.f90` (calculated using temperature gradients and isotropic thermal conductivity: $\mathbf{q} = -k \nabla T$).
    *   **Usage:** `src/fl_therm.f90` integrates these fluxes along element boundaries to update nodal temperatures.

---

## `vel`

*   **Meaning:** Nodal velocities.
*   **Shape:** `(nz, nx, 2)`
*   **Dimension 1 (`nz`):** Vertical index (Z direction), from 1 to `nz`.
*   **Dimension 2 (`nx`):** Horizontal index (X direction), from 1 to `nx`.
*   **Dimension 3 (2):** Velocity component:
    *   `1`: Vx (horizontal velocity) [m/s]
    *   `2`: Vz (vertical velocity) [m/s]
*   **Update Locations:**
    *   **Main Update:** `src/fl_node.f90` (calculated from forces using momentum balance and mass scaling).
    *   **Boundary Conditions:** Overridden in `src/fl_node.f90` using the `bc` array.

---

## `stress0`

*   **Meaning:** Stress components in the elements.
*   **Shape:** `(nz-1, nx-1, 4, 4)`
*   **Dimension 1 (`nz-1`):** Vertical index (element), from 1 to `nz-1` (`nez`).
*   **Dimension 2 (`nx-1`):** Horizontal index (element), from 1 to `nx-1` (`nex`).
*   **Dimension 3 (4):** Stress components:
    *   `1`: Sxx (horizontal normal stress) [Pa]
    *   `2`: Szz (vertical normal stress) [Pa]
    *   `3`: Sxz (shear stress) [Pa]
    *   `4`: Syy (out-of-plane stress) [Pa] (if `ndim`=3 is used)
*   **Dimension 4 (4):** Sub-triangle index (each quadrilateral element is divided into 4 triangles in the FLAC method).
*   **Update Locations:**
    *   **Initialization:** `src/init_stress.f90` (typically lithostatic stress).
    *   **Constitutive Law:** `src/fl_rheol.f90` (updated based on strain rate and material rheology).
    *   **Rotation:** `src/fl_move.f90` (stresses are rotated to account for grid motion and rotation).

---

## `iphase`

*   **Meaning:** Dominant phase (material type) ID of each element.
*   **Shape:** `(nz-1, nx-1)`
*   **Dimension 1 (`nz-1`):** Vertical index (element).
*   **Dimension 2 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Initialization:** `src/init_phase.f90` (initial distribution or reading from `phasefile`).
    *   **Marker Update:** `src/marker_data.f90` (via `count_phase_ratio`, determined by the majority phase of markers within the element).
    *   **Phase Changes:** `src/change_phase.f90` (calls `count_phase_ratio` to re-evaluate elemental phase after markers change properties).
    *   **Remeshing:** Re-evaluated in `src/marker2elem.f90` after markers are redistributed.

---

## `phase_ratio`

*   **Meaning:** Volume fraction of each material phase within an element.
*   **Shape:** `(nphase, nz-1, nx-1)`
*   **Dimension 1 (`nphase`):** Phase ID index, from 1 to `nphase`.
*   **Dimension 2 (`nz-1`):** Vertical index (element).
*   **Dimension 3 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Main Update:** `src/marker_data.f90` (calculated in `count_phase_ratio` based on markers in the element).
    *   **Phase Changes:** `src/change_phase.f90` (directly modifies ratios for processes like serpentinization or eclogitization).
    *   **Remeshing:** Re-evaluated in `src/marker2elem.f90` after markers are redistributed.

---

## `aps`

*   **Meaning:** Accumulated Plastic Strain.
*   **Shape:** `(nz-1, nx-1)`
*   **Dimension 1 (`nz-1`):** Vertical index (element).
*   **Dimension 2 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Initialization:** `src/init_phase.f90` (can set initial weak zones).
    *   **Main Update:** `src/fl_rheol.f90` (accumulates when the stress state exceeds the Mohr-Coulomb yield criterion).
    *   **Healing:** `src/fl_rheol.f90` also implements an exponential decay (healing) if `tau_heal > 0`.
    *   **Remeshing:** the plastic strain of incoming material is reset to 0.

---

## `visn`

*   **Meaning:** Effective viscosity of the element [Pa-s].
*   **Shape:** `(nz-1, nx-1)`
*   **Dimension 1 (`nz-1`):** Vertical index (element).
*   **Dimension 2 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Initialization:** `src/init_visc.f90`.
    *   **Main Update:** `src/fl_rheol.f90` (calculated based on temperature, strain rate, and phase-specific rheological parameters).

---

## `force`

*   **Meaning:** Net nodal force [Nt].
*   **Shape:** `(nz, nx, 2)`
*   **Dimension 1 (`nz`):** Vertical index (node).
*   **Dimension 2 (`nx`):** Horizontal index (node).
*   **Dimension 3 (2):** Force component:
    *   `1`: Fx (horizontal force)
    *   `2`: Fz (vertical force, includes gravity)
*   **Update Locations:**
    *   **Main Update:** `src/fl_node.f90` (sum of internal stress divergence, gravity, and external boundary conditions).
    *   **Boundary Conditions:** `src/bc_update.f90` (applies forces from stress-based BCs).
    *   **Damping:** `src/fl_node.f90` applies numerical damping (`demf`) to the net force.

---

## `amass`

*   **Meaning:** Nodal inertial mass used for scaling the explicit time-stepping scheme (mass scaling) to ensure numerical stability for a given time step.
*   **Shape:** `(nz, nx)`
*   **Dimension 1 (`nz`):** Vertical index (node).
*   **Dimension 2 (`nx`):** Horizontal index (node).
*   **Update Locations:**
    *   **Inertial Mass (`amass`):** `src/dt_mass.f90` (calculated to ensure stability for a given time step, often using "mass scaling" to speed up static simulations).

---

## `rmass`

*   **Meaning:** Physical (real) nodal mass integrated from the densities of surrounding elements [kg].
*   **Shape:** `(nz, nx)`
*   **Dimension 1 (`nz`):** Vertical index (node).
*   **Dimension 2 (`nx`):** Horizontal index (node).
*   **Update Locations:**
    *   **Real Mass (`rmass`):** `src/rmasses.f90` (integrated from the densities of surrounding elements).

---

## `strain`

*   **Meaning:** Accumulated strain components in the element.
*   **Shape:** `(nz-1, nx-1, 3)`
*   **Dimension 1 (`nz-1`):** Vertical index (element).
*   **Dimension 2 (`nx-1`):** Horizontal index (element).
*   **Dimension 3 (3):** Strain components:
    *   `1`: Exx (horizontal normal strain)
    *   `2`: Ezz (vertical normal strain)
    *   `3`: Exz (shear strain)
*   **Update Locations:**
    *   **Main Update:** `src/fl_rheol.f90` (integrated from the element-averaged strain rate: `strain = strain + dt * strainr_avg`).
    *   **Rotation:** `src/fl_move.f90` (strains are rotated to account for grid rotation).

---

## `dvol`

*   **Meaning:** Relative volume change (volumetric strain increment) of each sub-triangle during a single time step ($\Delta V / V_0$).
*   **Shape:** `(nz-1, nx-1, 4)`
*   **Dimension 1 (`nz-1`):** Vertical index (element), from 1 to `nz-1` (`nez`).
*   **Dimension 2 (`nx-1`):** Horizontal index (element), from 1 to `nx-1` (`nex`).
*   **Dimension 3 (4):** Sub-triangle index within the element.
*   **Update Locations:**
    *   **Main Update:** `src/fl_move.f90` (calculated as `det * area_old - 1`, where `det` is twice the new area of the sub-triangle, and `area_old` is the inverse of twice the area at the start of the step).
    *   **Usage:** Read in the constitutive equations in `src/fl_rheol.f90` (variable `dv`) to calculate volumetric deformation and isotropic pressure changes.

---

## `e2sr`

*   **Meaning:** Second invariant of the deviatoric strain rate, averaged over time [1/s].
*   **Shape:** `(nz-1, nx-1)`
*   **Update Locations:**
    *   **Main Update:** `src/fl_rheol.f90` (averaged over `ifreq_avgsr` time steps).

---

## `se2sr`

*   **Meaning:** Accumulator array that integrates strain rate components over the averaging time window (`ifreq_avgsr` time steps).
*   **Shape:** `(nz-1, nx-1, 3)`
*   **Dimension 1 (`nz-1`):** Vertical index (element).
*   **Dimension 2 (`nx-1`):** Horizontal index (element).
*   **Dimension 3 (3):** Strain rate component accumulator:
    *   `1`: Integrated horizontal normal strain rate $\int \dot{\epsilon}_{xx} \, dt$
    *   `2`: Integrated vertical normal strain rate $\int \dot{\epsilon}_{zz} \, dt$
    *   `3`: Integrated shear strain rate $\int \dot{\epsilon}_{xz} \, dt$
*   **Update Locations:**
    *   **Main Accumulation:** `src/fl_rheol.f90` (adds the current step's strain increment: `se2sr = se2sr + sr_loc * dt`).
    *   **Averaging & Reset:** `src/fl_rheol.f90` computes the second invariant of strain rate `e2sr` at the end of the averaging interval (`dtavg`) using the accumulated components in `se2sr`, then resets `se2sr` to `0.d0` (or initialized to `1d-16` in `setflac.f90`).

---

## `area`

*   **Meaning:** Geometric property of sub-triangles. Specifically, it stores the **inverse of twice the area** ($1 / (2 \times \text{Area})$) of each sub-triangle [1/m^2].
*   **Shape:** `(nz-1, nx-1, 4)`
*   **Dimension 1 (`nz-1`):** Vertical index (element).
*   **Dimension 2 (`nx-1`):** Horizontal index (element).
*   **Dimension 3 (4):** Sub-triangle index.
*   **Update Locations:**
    *   **Initialization:** `src/init_areas.f90`.
    *   **Main Update:** `src/fl_move.f90` (re-calculated as the grid nodes move).

---

## `source`

*   **Meaning:** Internal heat source term (e.g., radiogenic heating) [W/kg].
*   **Shape:** `(nz-1, nx-1)`
*   **Update Locations:**
    *   **Initialization:** `src/init_temp.f90` (calculates depth-dependent radiogenic heating using `hs` and `hr`).

---

## `shrheat`

*   **Meaning:** Shear heating (viscous dissipation) term, averaged over time.
*   **Shape:** `(nz-1, nx-1)`
*   **Update Locations:**
    *   **Main Update:** `src/fl_rheol.f90` (calculated as $\sigma_{ij} \dot{\epsilon}_{ij}$ and averaged over `ifreq_avgsr` time steps).

---

## `sshrheat`

*   **Meaning:** Accumulator array that integrates viscous shear heating (dissipation) [W/m³] over the averaging time window (`ifreq_avgsr` time steps).
*   **Shape:** `(nz-1, nx-1)`
*   **Dimension 1 (`nz-1`):** Vertical index (element).
*   **Dimension 2 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Main Accumulation:** `src/fl_rheol.f90` (adds viscous dissipation power density: `sshrheat = sshrheat + shear_heating * dt`).
    *   **Averaging & Reset:** `src/fl_rheol.f90` computes the time-averaged shear heating `shrheat` at the end of the averaging interval as `sshrheat / dtavg`, then resets `sshrheat` to `0.d0`.

---

## `dtopo`

*   **Meaning:** Change in topography at the surface nodes due to surface processes during a time step [m].
*   **Shape:** `(nx)`
*   **Dimension 1 (`nx`):** Horizontal index (node).
*   **Update Locations:**
    *   **Main Update:** `src/fl_move.f90` (calculated from the diffusion-based erosion/sedimentation model: `dtopo = dt * diffusion_term`).

---

## `dhacc`

*   **Meaning:** Accumulated non-tectonic surface elevation change ($\Delta h$) for each surface element since the last marker resurfacing (addition/deletion) [m]. Saved in `dhacc.rs`.
*   **Shape:** `(nx-1)`
*   **Dimension 1 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Main Accumulation:** `src/fl_move.f90` (accumulates element-averaged `dtopo` and `extrusion`).
    *   **Marker Addition/Deletion & Reset:** `src/fl_move.f90` (`resurface` adds/deletes markers when `dhacc(i)` exceeds threshold and resets `dhacc(i) = 0.d0`).

---

## `dhacc_correction`

*   **Meaning:** Accumulated non-tectonic surface elevation displacement ($\Delta h$) at each surface node since the last remeshing (from erosion, sedimentation, and volcanic extrusion) [m]. Saved in `dhacc_corr.rs`.
*   **Shape:** `(nx)`
*   **Dimension 1 (`nx`):** Horizontal index (node).
*   **Update Locations:**
    *   **Main Accumulation:** `src/fl_move.f90` (accumulates nodal `dtopo` and `extrusion`).
    *   **Reset After Remeshing:** `src/remesh.f90` (resets `dhacc_correction(:) = 0.d0` after remeshing and marker re-mapping complete).

---

## `extrusion`

*   **Meaning:** Thickness of volcanic material added to the surface elements (e.g., arc volcanism) [m].
*   **Shape:** `(nx-1)`
*   **Dimension 1 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Main Update:** `src/fl_move.f90` (calculated from the total melt produced in the mantle wedge during a time step).

---

## `extr_acc`

*   **Meaning:** Accumulated thickness of extruded (volcanic) material added to the surface elements [m]. Tracks total magmatic extrusion volume per column.
*   **Shape:** `(nx-1)`
*   **Dimension 1 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Initialization:** `src/init_cord.f90` (initialized to `0.d0`).
    *   **Main Accumulation:** `src/fl_move.f90` (adds current step's extrusion: `extr_acc(i) = extr_acc(i) + extrusion(i)`).
    *   **Remeshing:** Reset to `0.d0` in `src/fl_move.f90` when topography is updated, or interpolated onto the new grid in `src/remesh.f90`.
    *   **Restart:** Saved in `src/saveflac.f90` and loaded in `src/rsflac.f90`.

---

## `fmelt`

*   **Meaning:** Volume fraction of melt produced in an element during a time step.
*   **Shape:** `(nz-1, nx-1)`
*   **Update Locations:**
    *   **Main Update:** `src/change_phase.f90` (calculated in the mantle wedge based on P-T conditions and the presence of hydrated phases).

---

## `fmagma`

*   **Meaning:** Local volume fraction of magma residing within the mantle element.
*   **Shape:** `(nz-1, nx-1)`
*   **Update Locations:**
    *   **Main Update:** `src/fl_therm.f90` (increases due to melt migration from `fmelt` and decreases due to freezing or extraction to the surface).

---

## `bc`

*   **Meaning:** Prescribed boundary conditions for nodal velocities [m/s].
*   **Shape:** `(nz, nx, 2)`
*   **Dimension 1 (`nz`):** Vertical index (node).
*   **Dimension 2 (`nx`):** Horizontal index (node).
*   **Dimension 3 (2):** Boundary condition component:
    *   `1`: Prescribed Vx (horizontal velocity)
    *   `2`: Prescribed Vz (vertical velocity)
*   **Update Locations:**
    *   **Initialization:** `src/init_bc.f90` (calculates and assigns velocity values based on the `.inp` parameter file's boundary conditions).
    *   **Application:** `src/fl_node.f90` (overrides calculated nodal velocities with these prescribed values).

---

## `ncod`

*   **Meaning:** Nodal boundary condition indicators for velocities/degrees of freedom. It specifies if a node has a prescribed velocity boundary condition applied.
*   **Shape:** `(nz, nx, 2)`
*   **Dimension 1 (`nz`):** Vertical index (node).
*   **Dimension 2 (`nx`):** Horizontal index (node).
*   **Dimension 3 (2):** Nodal velocity component direction:
    *   `1`: Horizontal velocity indicator ($V_x$)
    *   `2`: Vertical velocity indicator ($V_z$)
*   **Values:**
    *   `0`: Free node (velocities are calculated using the momentum balance equations).
    *   `1` (or other non-zero flags like `10`, `30`): Prescribed boundary velocity (momentum solver is bypassed, and velocity is forced to the value in `bc`).
*   **Update Locations:**
    *   **Initialization:** `src/init_bc.f90` (reads `.inp` parameters and sets indicator to `1` for nodes within boundary conditions).
    *   **Usage:** Checked in `src/fl_node.f90` when solving the equations of motion to decide whether to apply the boundary velocity `bc`.

---

## `bcstress`

*   **Meaning:** Prescribed stress (normal or shear) boundary conditions applied on domain boundary segments [Pa].
*   **Shape:** `(nbou, 3)` where `nbou = ((nz-1)+(nx-1))*2`.
*   **Dimension 1 (`nbou`):** Index of the boundary segment where stress is applied (up to `nopbmax` active segments).
*   **Dimension 2 (3):** Stress component:
    *   `1`: Prescribed normal stress component [Pa]
    *   `2`: Prescribed shear stress component [Pa]
    *   `3`: Prescribed out-of-plane shear component (unused/commented)
*   **Update Locations:**
    *   **Initialization:** `src/init_bc.f90` (assigns stress values based on prescribed boundary functions in the `.inp` file).
    *   **Application:** `src/bc_update.f90` (converts these prescribed stresses into equivalent nodal forces in the `force` array).

---

## `nopbou`

*   **Meaning:** Segment endpoint mapping array for boundaries with applied stress (force) boundary conditions.
*   **Shape:** `(nbou, 4)` where `nbou = ((nz-1)+(nx-1))*2`.
*   **Dimension 1 (`nbou`):** Boundary segment index (up to `nopbmax` active segments).
*   **Dimension 2 (4):** Indices of the endpoint nodes of the boundary segment:
    *   `1`: Horizontal index ($i_1$) of node 1
    *   `2`: Horizontal index ($i_2$) of node 2
    *   `3`: Vertical index ($j_1$) of node 1
    *   `4`: Vertical index ($j_2$) of node 2
*   **Update Locations:**
    *   **Initialization:** `src/init_bc.f90` (maps boundary segments and endpoints based on boundary regions defined in the `.inp` parameter file).
    *   **Usage:** Used in `src/bc_update.f90` to compute segment length, orientation, and distribute stress-induced forces to the endpoints.

---

## `ncodbou`

*   **Meaning:** Applied stress boundary condition component indicators.
*   **Shape:** `(nbou, 3)` where `nbou = ((nz-1)+(nx-1))*2`.
*   **Dimension 1 (`nbou`):** Boundary segment index.
*   **Dimension 2 (3):** Indicator flag for active stress components on the segment:
    *   `1`: Prescribed normal stress component (flag = `1`).
    *   `2`: Prescribed shear stress component (flag = `1`).
    *   `3`: Prescribed out-of-plane shear component (unused/commented).
*   **Update Locations:**
    *   **Initialization:** `src/init_bc.f90` (sets the flags based on stress boundary types in `.inp` file).
    *   **Usage:** Reference indicator in `src/bc_update.f90`.

---

## `jmoho`

*   **Meaning:** Vertical index (`j`) of the element representing the Moho (crust-mantle boundary).
*   **Shape:** `(nx-1)`
*   **Dimension 1 (`nx-1`):** Horizontal index (element).
*   **Update Locations:**
    *   **Main Update:** `src/marker2elem.f90` (determined by finding the shallowest element in each column where the `phase_ratio` of mantle-type phases exceeds 0.5).

---

## `nmark_elem`

*   **Meaning:** The count/number of active Lagrangian markers located within each grid element.
*   **Shape:** `(nz-1, nx-1)`
*   **Dimension 1 (`nz-1`):** Vertical element index.
*   **Dimension 2 (`nx-1`):** Horizontal element index.
*   **Update Locations:**
    *   **Re-indexing:** `src/marker2elem.f90` (markers are re-binned into their containing elements after advection).
    *   **Creation/Deletion:** `src/marker_data.f90` (updated via `add_marker` or when markers are removed).

---

## `mark_id_elem`

*   **Meaning:** Global indices (IDs) of the Lagrangian markers located within each grid element.
*   **Shape:** `(max_markers_per_elem, nz-1, nx-1)`
*   **Dimension 1 (`max_markers_per_elem`):** Local list of markers inside the element (up to 32).
*   **Dimension 2 (`nz-1`):** Vertical element index.
*   **Dimension 3 (`nx-1`):** Horizontal element index.
*   **Update Locations:**
    *   **Re-indexing:** `src/marker2elem.f90` (markers are re-binned into their containing elements after advection).
    *   **Creation/Deletion:** `src/marker_data.f90` (updated via `add_marker` or when markers are removed).

---

## `mark_x`

*   **Meaning:** Global Eulerian X-coordinates of Lagrangian markers [m]. Accurate only during remeshing.
*   **Shape:** `(max_markers)`
*   **Dimension 1:** Total marker index.
*   **Update Locations:**
    *   **Initialization:** `src/init_marker.f90`.
    *   **Main Update:** `src/fl_move.f90` (advection using interpolated nodal velocities).
    *   **Creation:** `src/marker2elem.f90` (fills element gaps).

---

## `mark_y`

*   **Meaning:** Global Eulerian Z-coordinates (depth/vertical) of Lagrangian markers [m]. Accurate only during remeshing. Note that in the codebase `mark_y` corresponds to the vertical coordinate.
*   **Shape:** `(max_markers)`
*   **Dimension 1:** Total marker index.
*   **Update Locations:**
    *   **Initialization:** `src/init_marker.f90`.
    *   **Main Update:** `src/fl_move.f90` (advection using interpolated nodal velocities).
    *   **Creation:** `src/marker2elem.f90` (fills element gaps).

---

## `mark_a1`

*   **Meaning:** First local barycentric coordinate of markers within a specific element (relative distance to horizontal element boundaries).
*   **Shape:** `(max_markers)`
*   **Dimension 1:** Total marker index.
*   **Update Locations:**
    *   **Main Update:** `src/marker2elem.f90` (calculated whenever markers move or grid is reorganized).

---

## `mark_a2`

*   **Meaning:** Second local barycentric coordinate of markers within a specific element (relative distance to vertical element boundaries).
*   **Shape:** `(max_markers)`
*   **Dimension 1:** Total marker index.
*   **Update Locations:**
    *   **Main Update:** `src/marker2elem.f90` (calculated whenever markers move or grid is reorganized).

---

## `mark_phase`

*   **Meaning:** Material phase ID assigned to each marker.
*   **Type:** `integer(1)` (signed 1 byte). Phase IDs must stay within 1..127, so `maxph` (`src/params.f90`) cannot exceed 127.
*   **Shape:** `(max_markers)`
*   **Dimension 1:** Total marker index.
*   **Update Locations:**
    *   **Initialization:** `src/init_marker.f90`.
    *   **Phase Changes:** `src/change_phase.f90` (based on P-T history and composition).

---

## `mark_age`

*   **Meaning:** Creation time or geological age of each marker [second].
*   **Shape:** `(max_markers)`
*   **Update Locations:**
    *   **Creation:** `src/marker_data.f90` (assigned in `add_marker`).
    *   **Initialization:** `src/init_marker.f90`.

---

## `mark_dead`

*   **Meaning:** Status flag for marker activity.
*   **Type:** `integer(1)` (signed 1 byte).
*   **Value:** `1` for active, `0` for dead/inactive.
*   **Shape:** `(max_markers)`
*   **Update Locations:**
    *   **Activation:** `src/marker_data.f90` (set to `1` in `add_marker`).
    *   **Deactivation:** `src/fl_move.f90` (set to `0` if a marker moves outside the domain boundaries).

---

## `mark_ntriag`

*   **Meaning:** Global flattened triangle index identifying which sub-triangle of which element currently contains the marker: `2*(nx-1)*(nz-1)` possible values, decoded as `k = mod(n-1,2)+1` (sub-triangle 1 or 2), `j = mod((n-k)/2, nz-1)+1` (element row), `i = (n-k)/2/(nz-1)+1` (element column). See e.g. `src/change_phase.f90` or `src/lpeuler2bar.f90` for the decode.
*   **Type:** Default `integer` (4 bytes). Kept full width deliberately -- its range scales with mesh size and can exceed 127 (or even 32767) on realistic meshes, so it cannot be narrowed the way `mark_dead`/`mark_phase` were.
*   **Shape:** `(max_markers)`
*   **Update Locations:**
    *   **Main Update:** `src/marker2elem.f90` (re-calculated during the search process to find which element/triangle contains the marker).

---

## `vel_flex_old`

*   **Meaning:** Flexure-induced vertical velocity from the previous timestep [m/s]. Used to prevent cumulative velocity inflation when applying flexural deflection shifts to the entire grid column.
*   **Shape:** `(nx)`
*   **Dimension 1 (`nx`):** Horizontal index (node), from 1 (left) to `nx` (right).
*   **Update Locations:**
    *   **Initialization:** `src/arrays.f90` (allocated and initialized to `0.d0` in `allocate_arrays`).
    *   **Main Update:** `src/fl_flexure.f90` (subtracted from node velocities at the start of the step, recalculated as `vel_flex = dflex / dt` and applied to the column, and then stored back into `vel_flex_old` at the end of the step).
    *   **Restart/State:** Written to `vel_flex_old.rs` in `src/saveflac.f90` and restored from it in `src/rsflac.f90` (with fallback to `0.d0` for backwards compatibility).

---

## Remeshing and Temporary Workspace Arrays

The following arrays are allocated to manage dynamic grid restructuring (remeshing) or as short-term workspace inside subroutines. They are generally not part of the persistent physical state.

### Remeshing Grid Arrays

*   **`cordo(nz, nx, 2)`:** Stores the coordinates of the old grid before restructuring.
*   **`cold(nz, nx, 2)` / `cnew(nz, nx, 2)`:** Temporary arrays tracking coordinate mapping during interpolation.
*   **`pt((nz-1)*(nx-1)*2, 2, 3)`:** Sub-triangle coordinates used for Delaunay triangulation or element interpolation.
*   **`barcord(nz, nx, 3)`:** Barycentric interpolation weights on the restructured grid.
*   **`numtr(nz, nx)`:** Auxiliary indexing array for grid triangles.
*   **`dhnew(nx-1)` / `extnew(nx-1)`:** Workspace for interpolating topography and extrusion history onto the new elements.

### Temporary Subroutine Workspace

*   **`dummyn(nz, nx)`:** Nodal real-valued workspace buffer.
*   **`dummye(nz-1, nx-1)`:** Element real-valued workspace buffer.
*   **`stmpn(max(nx,nz))`:** Single-dimension real-valued temp array.
*   **`itmp(nz, nx)`:** Nodal integer-valued workspace buffer.
