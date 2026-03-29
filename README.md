# README Scientific Computing — Assignment Set 3

## Challenge A — Kármán Vortex Street

---

### Overview

This notebook implements the **Kármán vortex street** benchmark for 2D
incompressible flow past a circular cylinder using three independent numerical
methods:

1. **Finite Difference (FD)** — explicit Euler, first-order upwind advection, Jacobi pressure solver
2. **Lattice Boltzmann (LBM)** — D2Q9 BGK, Zou-He inlet, half-way bounce-back
3. **Finite Element (FEM)** — NGSolve, Taylor-Hood P2/P1, semi-implicit time stepping

Each method is validated at Re=100 and pushed to its maximum stable Reynolds
number via a coarse scan followed by binary search.

---

### Requirements

Install all dependencies with:

```bash
pip install numpy matplotlib ngsolve
```

| Package       | Purpose                       |
|---------------|-------------------------------|
| `numpy`       | Array operations (FD, LBM)    |
| `matplotlib`  | All plots                     |
| `ngsolve`     | FEM solver and mesh generation|

> **Note:** NGSolve includes `netgen` (mesh generator) as a dependency.
---

### Problem Setup

| Parameter         | Value |
| Channel dimensions| 2.2 m × 0.41 m |
| Cylinder centre   | (0.20, 0.20) m |
| Cylinder diameter | D = 0.10 m |
| Inlet velocity    | U_in = 1.0 m/s |
| Fluid density     | $\rho$ = 1.0 kg/$m^3$ |
| Reynolds number   | Re = U_in · D / ν |

Boundary conditions:
- **Inlet** (x=0): uniform flow u=U_in, v=0
- **Outlet** (x=L): zero-gradient
- **Top/bottom walls**: no-slip (u=v=0)
- **Cylinder surface**: no-slip (u=v=0)

---

### Structure

The notebook is organized in three self-contained sections, one per method.
Each section follows the same structure:

```
1. Parameters and grid/mesh setup
2. Core solver functions
3. Validation run at Re=100  ->  flow field plot + wake probe plot
4. Re_max sweep              ->  coarse scan + binary search
5. Final run at Re_max       ->  flow field plot + wake probe plot
```

---

### How to Run

Run the notebook cells **top to bottom**. Each method section is independent.


---

### Key Results

| Method    | Re_max    | Stability limit                   | Notes                     |
|-----------|-----------|-----------------------------------|---------------------------|
| FD        | N/A       | Artificial (numerical diffusion)  | Excluded from competition |
| LBM       | **323**   | BGK hard limit $\tau$ > 0.5       | Genuine physical limit    |
| FEM       | 203       | Explicit advection CFL            | Soft limit; improvable    |

---

## Output Figures
 
For each method, the code generates the following plots:
 
- **u-velocity field** (contour plot over the full domain): shows the
  recirculation zone behind the cylinder and the alternating wake structure.
  Produced at Re=100 and at Re_max.
- **Velocity magnitude field** (contour plot over the full domain): shows
  the speed distribution across the channel, highlighting the wake and
  the acceleration around the cylinder.
- **Wake probe signal** (time series of the transverse velocity component
  at a point two diameters downstream of the cylinder centre): periodic
  oscillations confirm vortex shedding onset.
 
All three plots are produced twice per method — once at Re=100 (validation)
and once at Re_max (stability limit). The FEM section additionally reports
the $L^2$ norm of the velocity field at regular intervals during the time loop.

---

### Notes

- The FD and LBM solvers use a **staircase approximation** for the cylinder
  boundary. The FEM solver uses **exact curved P3 elements** via `mesh.Curve(3)`.
- The FD solver is intentionally excluded from the Re_max competition because
  its first-order upwind scheme introduces numerical diffusion
  $ν_num = U_in·\Delta x/2$, which artificially prevents divergence at all tested Re.
- The LBM Re_max of 323 is a physical stability limit of the BGK model
  $(\tau → 0.5^+)$. Pushing beyond this would require an MRT or entropic LBM scheme.
- The FEM Re_max of 203 is a soft limit set by the explicit advection CFL
  (dt < h_min/U_in). With a more relaxed scheme e.g. Crank-Nicolson scheme we should get
  an even higher Re_max



## Challenge B — Optimal WiFi Router Placement

### Overview

This module finds the best position for a WiFi router in a 10×8 m floor plan
by solving the 2D Helmholtz equation and maximising the total signal strength
across four measurement points. The governing equation is:

```
$\Delta u + k^2 u = f$
```

where `u` is the complex wave field, `k` is the local wavenumber (material-dependent),
and `f` is a Gaussian source centred at the router position.

---


### Requirements

Install dependencies with:

```bash
pip install numpy scipy matplotlib
```

| Package    | Purpose                        |
|------------|--------------------------------|
| numpy      | Array operations, grid setup   |
| scipy      | Sparse matrix, LU factorisation|
| matplotlib | All plots                      |

---

### Problem Setup

#### Floor Plan
- Domain: 10 x 8 m
- Wall thickness: 0.15 m uniformly
- Internal walls: as defined in Figure 5 of the assignment

#### Material Properties

| Region | Refractive Index  | Effect                        |
|--------|-------------------|-------------------------------|
| Air    | 1.0               | Standard wave propagation     |
| Walls  | 2.5 + 0.5j        | Slower propagation + absorption|

#### Measurement Points

| Room        | X [m] | Y [m] |
|-------------|-------|-------|
| Living Room | 1     | 5     |
| Kitchen     | 2     | 1     |
| Bathroom    | 9     | 1     |
| Bedroom 1   | 9     | 7     |

Signal strength at each point is computed as the mean of $|u|^2$ over a
circle of radius r = 5 cm centred on that point.

#### Hard Constraint
The router must be placed at least **0.5 m** from every measurement point
and must not coincide with a wall cell.

---

### Method

#### Wavenumber Scaling
The physical WiFi frequency (2.4 GHz, $\lambda$ ≈ 12.5 cm) would require a mesh
of ~500,000 nodes to resolve accurately. Instead, the frequency is scaled
by 1/3 to **f_eff = 0.8 GHz** ($\lambda$ = 37.5 cm), reducing the problem to
**32,361 unknowns** (201 × 161 grid, h = 5 cm = $\lambda$/7.5). All wavenumbers
(k0, $k^2$_air, $k^2$_wall) are scaled consistently to preserve material contrast
and wall attenuation.

#### Solver
The 5-point finite difference stencil is assembled into a sparse complex
matrix **A** $\in$ $C^{(N×N)}$ with 159,641 non-zero entries. An impedance
(absorbing) boundary condition $\frac{\partial u}{\partial n} - i k_0 u = 0$ is applied on all outer
walls. The matrix is factorised **once** with `scipy.sparse.linalg.splu`
(LU decomposition); every subsequent solve for a new router position is a
cheap forward/back substitution — no re-factorisation needed.

#### Optimisation: Three-Level Hierarchical Grid Search

| Level      | Step   | Search region                               | Candidates |
|------------|--------|---------------------------------------------|------------|
| Coarse     | 0.5 m  | Full domain                                 | 248        |
| Fine       | 0.1 m  | $\pm$ 1.5 m window around coarse optimum    | 338        |
| Ultra-fine | 0.01 m | $\pm$ 0.3 m window around fine optimum      | 1759       |

Total solves: **2345** (vs ~100,000 for brute force at 1 cm resolution).

---

### Results

| Parameter         | Value                  |
|-------------------|------------------------|
| Optimal position  | (9.05, 7.75) m         |
| Total signal S*   | 1721.97                |
| vs. house centre  | +44.5 dB (x28,000)     |

#### Per-room signal at optimum

| Room        | ⟨$\|u\|^2$⟩         | Relative [dB] |
|-------------|--------------------|---------------|
| Living Room | 7.37 x $10^{-1}$   | −33.7         |
| Kitchen     | 1.28 x $10^{-1}$   | −41.3         |
| Bathroom    | 1.40 x $10^{0}$    | −30.9         |
| Bedroom 1   | 1.72 x $10^{3}$    |   0.0         |

---

### How to Run

```bash
python challenge_b.py
```

The script runs end-to-end and produces the following plots in sequence:

1. **Floor plan** — wall geometry and measurement point exclusion zones
2. **Coarse search map** — signal landscape at 0.5 m resolution
3. **Fine search map** — zoomed in at 0.1 m resolution
4. **Ultra-fine search map** — final zoom at 0.01 m resolution
5. **Optimal signal field** — dB heatmap at the optimal router position
6. **Baseline signal field** — dB heatmap at house centre (5.0, 4.0) m
7. **Bar chart** — per-room signal comparison: baseline vs optimal

---

### Key Parameters (top of script)

```python
freq_eff   = 2.4e9 / 3.0    # Scaled frequency [Hz] 
h          = 0.05           # Grid spacing [m]
A_source   = 1.0e4          # Gaussian source amplitude
sigma      = 0.2            # Gaussian source width [m]
meas_r     = 0.05           # Averaging radius at measurement points [m]
router_exclusion = 0.5      # Minimum distance from measurement points [m]
```

---

### Notes and Limitations

- **Wavenumber scaling** is the dominant approximation. The optimal position
  at full 2.4 GHz frequency may differ due to denser interference fringes.
- **Averaging circle resolution**: at h = 5 cm, the r = 5 cm measurement
  circle contains only ~3–5 grid nodes, making the signal estimate sensitive
  to exact router placement on the grid.
- **Objective imbalance**: the unweighted sum S is dominated by Bedroom 1
  (>99.9% of S*). A weighted formulation would give more balanced
  coverage across all rooms.