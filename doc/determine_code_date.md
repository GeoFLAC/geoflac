# Guide: Determining the Target Code Date from `snapshot.diff`

When using the parameter upgrade utility `util/update_inp.py`, you must provide the approximate target date of the GeoFLAC engine version that the old parameter file was written for. If you only have a `snapshot.diff` (or patch file), this guide explains the two systematic methods to extract the correct code date.

---

## Method 1: The Direct Commit Hash Lookup (Most Precise)

If your `snapshot.diff` was generated using `git diff` or `git format-patch`, the diff headers contain baseline commit hashes.

### Step 1: Locate baseline commit hashes in the diff
Open `snapshot.diff` and look for the file diff index headers:
```diff
diff --git a/src/read_params.f90 b/src/read_params.f90
index 5b4b263..c2af4d2 100644
--- a/src/read_params.f90
+++ b/src/read_params.f90
```
In the `index 5b4b263..c2af4d2` line:
* `5b4b263` is the baseline (old) git commit SHA-1 prefix before changes.
* `c2af4d2` is the new commit SHA-1 prefix.

### Step 2: Query git for the commit date
Run the following git command in the repository to extract the commit's calendar date:
```bash
git show -s --format=%ad --date=short 5b4b263
```
* **Expected Output:** A date string in `YYYY-MM-DD` format (e.g., `2026-05-27`).
* Use this exact date string as the second parameter in `util/update_inp.py`!

---

## Method 2: Syntax-Based Feature Matching (Structural Heuristic)

If the commit hashes are not in your repository (or if `snapshot.diff` is an ordinary unified diff file without git metadata), you can identify the epoch by looking at which variables are read or written in `snapshot.diff` and comparing them against our reference mapping:

### 1. Check the Rheology parameters block
Examine the `src/read_params.f90` section of the diff:
* **Has `vactiv`:** If the rheology line reads `vactiv` (25 parameters), the codebase date is **on or after 2025-10-17**.
* **Does not have `vactiv`:** If the rheology line does not read `vactiv` (24 parameters), the codebase date is **before 2025-10-17**.

### 2. Check the Thermal parameters block
* **Has `extra_pres`:** If `extra_pres` is read right after `i_prestress`, the date is **on or after 2025-10-17**.
* **Does not have `extra_pres`:** The date is **before 2025-10-17**.

### 3. Check the Magma / Melting parameters block
Find the melting lines inside `src/read_params.f90` in the diff:
* **Has `prod_magma2`:** If the line has `itype_melting, nelem_serp, nelem_dike, prod_magma, prod_magma2, rho_magma` (6 values), the date is **on or after 2026-07-30**.
* **Has `nelem_dike` but NOT `prod_magma2`:** If the line has `itype_melting, nelem_serp, nelem_dike, prod_magma, rho_magma` (5 values), the date is **between 2026-05-27 and 2026-07-30**.
* **Has `rho_magma` but NOT `nelem_dike`:** If it reads `itype_melting, nelem_serp, prod_magma, rho_magma` (4 values), the date is **between 2023-02-01 and 2026-05-27**.
* **Has `prod_magma` but NOT `rho_magma`:** If it reads `itype_melting, nelem_serp, prod_magma` (3 values), the date is **between 2021-08-10 and 2023-02-01**.
* **Has `arc_extrusion_rate`:** The date is **before 2021-08-10**.

### 4. Check the `nzone_age` column format
* **Split 3-Line layout with `nph_layer` on Line 2:** The date is **on or after 2023-02-01**.
* **Split 3-Line layout with `nph_layer` on Line 1:** The date is **between 2023-01-07 and 2023-02-01**.
* **Flat 1-Line layout with 15 values (`ictherm` present):** The date is **between 2023-01-05 and 2023-01-07**.
* **Flat 1-Line layout with 12 values (`ictherm` absent):** The date is **before 2023-01-05**.

---

## Epoch Quick-Reference Date Summary
Once you match the structural changes, select a date within the identified epoch range:

* **Before 2021-07-15**: Pre-tracer removal build.
* **2021-07-15**: Tracer block replaced by dummy lines.
* **2021-08-10**: Chamber max and extrusion rate renamed.
* **2023-01-05**: `ictherm` added to `nzone_age` columns.
* **2023-01-07**: `nzone_age` split into 3-line format.
* **2023-01-17**: Obsolete phase layers and tracer variables removed.
* **2023-02-01**: `nph_layer` moved to line 2; added `rho_magma` and `latent_heat_magma`.
* **2023-02-04**: `angle_mzone` renamed and simplified.
* **2025-10-17**: Added `extra_pres` and rheology `vactiv` parameter.
* **2026-05-27**: Added crustal dike width `nelem_dike`.
* **2026-07-30**: Added hydrous-mantle magma migration rate `prod_magma2`.
