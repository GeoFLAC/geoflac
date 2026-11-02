# GeoFLAC Developer Guide: Mixed Discretization & Element Triangulation

This guide explains the mathematical and numerical details of why and how each quadrilateral element in GeoFLAC is split into 4 sub-triangles, and how the shape function derivative arrays (`shpdx_loc` and `shpdz_loc`) are formulated and used.

---

## 1. Why Quadrilaterals are Split into 4 Triangles

In standard 2D finite-difference and finite-element modeling, simple 4-node quadrilateral elements with one integration point suffer from **hourglassing** (also known as grid-stepping or zero-energy modes). This is a numerical phenomenon where the grid can deform in an alternating "hourglass" pattern without changing volume or shear angles, producing zero strain energy and causing grid instability.

GeoFLAC solves this using the **Mixed Discretization** technique (Cundall, 1989). Instead of treating the quadrilateral as a single element, the code splits each quadrilateral into **two overlapping pairs of constant-strain triangles**:

```
Node numbering of quadrilateral:
 1 (j, i)   -------  3 (j, i+1)
            |     |
            |     |
 2 (j+1, i) -------  4 (j+1, i+1)

Diagonal / splitting (Pair 1):       Diagonal \ splitting (Pair 2):
       1-------3                            1
       | /                                  | \
       2                                    2---3
    Triangle A (1)                       Triangle C (3)

               3                            1---3
             / |                             \  |
           2---3                                2
    Triangle B (2)                       Triangle D (4)
```

### Key Numerical Benefits:
1.  **Overcoming Volumetric Locking**: Triangles have constant strain. In incompressible or plastic flows (like geology), quadrilaterals easily lock up. The mixed formulation allows independent dilatation and shear tracking.
2.  **Symmetry**: By overlaying both diagonals (Pair 1 split by `/` and Pair 2 split by `\`), the model is perfectly symmetric and does not favor shear deformation along one specific diagonal.
3.  **Hourglass Suppression**: Since any hourglass deformation mode in one pair of triangles is resisted by the other pair, zero-energy modes are naturally suppressed without needing artificial viscous damping forces.

---

## 2. Triangle Node Mapping (`kk` and `mm` parameters)

To map the vertices of the 4 sub-triangles (A, B, C, D) back to the 4 corners of the parent quadrilateral, the solver uses two static lookup arrays defined in `src/fl_node.f90` and `src/fl_rheol.f90`:

```fortran
integer, parameter :: kk(3, 4) = reshape((/ 2,3,4, 1,2,3, 1,2,4, 1,3,4 /), (/3, 4/))
integer, parameter :: mm(3, 4) = reshape((/ 3,3,2, 2,2,2, 3,1,3, 1,1,1 /), (/3, 4/))
```

*   **`kk(m, k)`**: Maps the local node index $m$ (1, 2, or 3) of sub-triangle $k$ (1, 2, 3, or 4) to the parent quadrilateral corner node (1, 2, 3, or 4).
*   **`mm(m, k)`**: Used in force integration loops to associate local nodal updates with corresponding corners.

---

## 3. Shape Function Derivative Arrays (`shpdx_loc` & `shpdz_loc`)

In constant-strain triangles, the displacement and velocity fields are linear functions of coordinates $x$ and $z$. The shape functions $N_m$ (for local node $m \in \{1, 2, 3\}$) are linear:

$$N_m(x, z) = a_m x + b_m z + c_m$$

The spatial derivatives of these shape functions are constant within the triangle and are defined as:

$$\frac{\partial N_m}{\partial x} = \frac{y_b - y_c}{2 \cdot \text{Area}} = (y_b - y_c) \cdot \text{area\_inv}$$

$$\frac{\partial N_m}{\partial z} = \frac{x_c - x_b}{2 \cdot \text{Area}} = (x_c - x_b) \cdot \text{area\_inv}$$

Where $b$ and $c$ are the other two nodes of the triangle in counter-clockwise order, and `area_inv` is the inverse double-area array `area(j, i, k)` (calculated in `src/init_areas.f90` as $1 / (2 \cdot \text{Area})$).

### Code Discretization (e.g., Triangle A / index 1):
For Triangle A (nodes 1, 2, 3):
```fortran
shpdx_loc(1, 1) = (y2 - y3) * area(je, ie, 1)  ! dN1/dx
shpdx_loc(2, 1) = (y3 - y1) * area(je, ie, 1)  ! dN2/dx
shpdx_loc(3, 1) = (y1 - y2) * area(je, ie, 1)  ! dN3/dx

shpdz_loc(1, 1) = (x3 - x2) * area(je, ie, 1)  ! dN1/dz
shpdz_loc(2, 1) = (x1 - x3) * area(je, ie, 1)  ! dN2/dz
shpdz_loc(3, 1) = (x2 - x1) * area(je, ie, 1)  ! dN3/dz
```

---

## 4. Usage in the Solver

These derivative components are used in two critical physical calculation steps:

### A. Strain Rate Calculation (`src/fl_rheol.f90`)
The components of the strain-rate tensor for each sub-triangle $k$ are calculated by summing node velocities multiplied by shape function derivatives:

$$\dot{\epsilon}_{xx}^k = \sum_{m=1}^3 V_x^m \cdot \text{shpdx\_loc}(m, k)$$

$$\dot{\epsilon}_{zz}^k = \sum_{m=1}^3 V_z^m \cdot \text{shpdz\_loc}(m, k)$$

$$\dot{\epsilon}_{xz}^k = \frac{1}{2} \sum_{m=1}^3 \left[ V_x^m \cdot \text{shpdz\_loc}(m, k) + V_z^m \cdot \text{shpdx\_loc}(m, k) \right]$$

### B. Nodal Force Integration (`src/fl_node.f90`)
After stresses $\boldsymbol{\sigma}^k$ are resolved on the triangles, the equivalent nodal forces $f_x$ and $f_z$ are integrated and assembled by multiplying stresses by shape function derivatives:

$$f_x^m = \sum_{k=1}^4 \left( \sigma_{xx}^k \cdot \text{shpdx\_loc}(m, k) + \sigma_{xz}^k \cdot \text{shpdz\_loc}(m, k) \right) \cdot \text{factor}$$

$$f_z^m = \sum_{k=1}^4 \left( \sigma_{xz}^k \cdot \text{shpdx\_loc}(m, k) + \sigma_{zz}^k \cdot \text{shpdz\_loc}(m, k) \right) \cdot \text{factor}$$

Where `factor` scales the force contributions to nodes according to sub-element partitions.
