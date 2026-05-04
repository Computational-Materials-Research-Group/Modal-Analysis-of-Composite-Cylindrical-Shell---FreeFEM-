# Modal Analysis of Composite Cylindrical Shell - FreeFEM++

<p align="center">
  <img src="https://img.shields.io/badge/FreeFEM++-Simulation-blue?style=for-the-badge&logo=gnu&logoColor=white"/>
  <img src="https://img.shields.io/badge/CFRP-Cylindrical%20Shell-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/CLT-A--Matrix-red?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Modal-Dynamic%20Animation-purple?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/ParaView-VTK%20Export-green?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey?style=for-the-badge"/>
</p>

<p align="center">
  A 2D finite element simulation of <b>modal analysis of a composite cylindrical shell</b>
  using FreeFEM++. Models a CFRP shell unrolled into a 2D rectangle, assembled using
  Classical Laminate Theory, with natural frequencies and mode shapes extracted via
  inverse power iteration with M-orthogonal deflation. Each mode is exported as a
  <i>dynamic PVD animation</i> — the shell physically oscillates in 3D cylindrical
  coordinates over one complete period.
</p>

<img width="1008" height="772" alt="COMPOSITE SHELL MODAL 1" src="https://github.com/user-attachments/assets/bc66707c-abbd-4c66-935b-45d9a69b7527" />


---

## Physics

Natural frequencies and mode shapes of composite cylindrical shells depend on
laminate stiffness, circumferential wave number n, and axial boundary conditions.
This simulation includes:

- 2D unrolled plane-stress elasticity with CLT effective stiffness
- Classical Laminate Theory A-matrix assembled from per-ply rotated Q-bar matrices
- Consistent mass matrix via varf on the unrolled rectangular domain
- Inverse power iteration with M-orthogonal deflation for multiple eigenvalues
- 3D cylindrical coordinate mapping for physically correct ParaView visualisation
- Dynamic PVD animation: each mode oscillates as phi(x)*sin(2*pi*f*t) over one period

---

## Geometry

```
  Unrolled cylinder (x=axial, y=circumferential arc length):

  y=Ly=2*pi*R  +---------------------------------------+  <- top (periodic)
               |                                       |
               |   CFRP shell, homogenised CLT moduli  |
               |   [0/90/0/90]_s  8-ply laminate       |
               |                                       |
  y=0          +---------------------------------------+  <- bottom (periodic)
               x=0                                   x=L
           Fixed end (clamped)                   Free end

  R = 100 mm,  L = 400 mm,  h = 1 mm
  Lx = L = 400 mm (axial)
  Ly = 2*pi*R = 628 mm (circumferential arc length)
```

- Boundary label 1 = bottom (y=0, periodic symmetry: uy=0)
- Boundary label 2 = right (x=L, free end)
- Boundary label 3 = top (y=Ly, periodic symmetry: uy=0)
- Boundary label 4 = left (x=0, clamped: ux=uy=0)

---

## Mode Notation

Cylindrical shell modes are described by (m; n):
- m = number of axial half-waves along the shell length
- n = number of circumferential lobes (waves around the circumference)

| Mode | Notation | Description |
|------|----------|-------------|
| 1 | (1;1) | First axial bending + 1 circumferential wave |
| 2 | (1;4) | First axial + 4-petal ovaling |
| 3 | (1;3) | First axial + 3-lobe cross section |
| 4 | (1;5) | First axial + 5-lobe cross section |
| 5 | (1;2) | First axial + oval (elliptic) |
| 6 | (2;5) | Second axial + 5-lobe |

---

## Material Parameters

### T300/5208 Carbon-Epoxy Ply Properties

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Longitudinal modulus | E1 | 181 | GPa |
| Transverse modulus | E2 | 10.3 | GPa |
| In-plane shear modulus | G12 | 7.17 | GPa |
| Major Poisson ratio | nu12 | 0.28 | — |
| Density | rhoP | 1600 | kg/m^3 |

### Laminate Stack [0/90/0/90]_s (8 plies)

| Ply | Angle | Thickness |
|-----|-------|-----------|
| 0 | 0 deg | 0.125 mm |
| 1 | 90 deg | 0.125 mm |
| 2 | 0 deg | 0.125 mm |
| 3 | 90 deg | 0.125 mm |
| 4 | 90 deg | 0.125 mm |
| 5 | 0 deg | 0.125 mm |
| 6 | 90 deg | 0.125 mm |
| 7 | 0 deg | 0.125 mm |

Total thickness: 1.0 mm. Symmetric and balanced: A16 = A26 = 0.

### CLT Effective Properties (computed automatically)

| Property | Value | Unit |
|----------|-------|------|
| Ceff11 | ~95.7 | GPa |
| Ceff22 | ~95.7 | GPa |
| Ceff12 | ~5.5 | GPa |
| Ceff66 | ~7.17 | GPa |

---

## Governing Equations

### Plane-Stress CLT Constitutive Law

```
sigma_xx = Ceff11*eps_xx + Ceff12*eps_yy
sigma_yy = Ceff12*eps_xx + Ceff22*eps_yy
sigma_xy = Ceff66*gamma_xy
```

where Ceffij = Aij / hTotal and Aij is the CLT A-matrix.

### CLT A-Matrix Assembly

```
For each ply k at angle theta_k:
  Qbar_ij(theta_k) via Tsai-Pagano rotation
  A_ij += Qbar_ij * t_k

[0/90/0/90]_s: A11 = A22 (balanced-symmetric quasi-isotropic)
```

### Stiffness and Mass Varfs (FreeFEM++ bilinear forms)

```
aStiff([ux,uy],[vx,vy]) =
  int2d(Th)(
    hShell*(
      Ceff11*dx(ux)*dx(vx) + Ceff12*dy(uy)*dx(vx)
    + Ceff12*dx(ux)*dy(vy) + Ceff22*dy(uy)*dy(vy)
    + Ceff66*(dy(ux)+dx(uy))*(dy(vx)+dx(vy))
    )
  )
  + on(4, ux=0., uy=0.)   // clamped left
  + on(1, uy=0.)           // periodic bottom
  + on(3, uy=0.)           // periodic top

aMass([ux,uy],[vx,vy]) =
  int2d(Th)( rhoP * hShell * (ux*vx + uy*vy) )
```

### Eigenvalue Problem

```
K * phi = omega^2 * M * phi

omega_i = sqrt(lambda_i)    [rad/s]
f_i     = omega_i / (2*pi)  [Hz]
```

### Inverse Power Iteration with M-Orthogonal Deflation

```
For mode i:
  1. Fill initial DOF vector vv[2k] = sin(m*pi*xk/L)*cos(m*pi*yk/Ly)
     (interleaved DOF ordering: 2k = ux_k, 2k+1 = uy_k)
  2. M-deflate: vv -= sum_j (phi_j' M vv / phi_j' M phi_j) * phi_j
  3. Iterate: ww = (K - sigma*M)^-1 * M * vv
  4. M-deflate ww against all previous modes
  5. Rayleigh quotient: lambda = ww'Kww / ww'Mww
  6. M-normalise: vv = ww / sqrt(ww'Mww)
  7. Repeat until |lambda_new - lambda_old| < 1e-8 * |lambda|
```

### Dynamic Animation

```
u(x, t_j) = ampMax * phi(x) * sin(2*pi*f*t_j)

3D cylindrical coordinate mapping per node k:
  theta  = yFlat_k / R
  X3d    = xFlat_k  + ampMax*phiX[k]*sin(2*pi*f*t)   (axial)
  Y3d    = (R + ampMax*phiY[k]*sin(2*pi*f*t))*sin(theta)
  Z3d    = (R + ampMax*phiY[k]*sin(2*pi*f*t))*cos(theta)

ampMax = 0.08 * R = 8 mm (8% of radius for visible motion)
Ntime  = 40 frames per period
```

---

## Numerical Method

| Aspect | Choice |
|--------|--------|
| Spatial discretisation | Finite Element Method (FEM) |
| Displacement element type | P1 (linear nodal), vector [P1,P1] |
| Time integration | Quasi-static modal (no transient dynamics) |
| Linear solver | UMFPACK (direct sparse factorisation) |
| Eigenvalue method | Inverse power iteration + M-orthogonal deflation |
| Number of modes | 12 |
| Convergence tolerance | 1e-8 relative change in Rayleigh quotient |
| Mesh | 80 axial x 120 circumferential elements |
| 3D output | Cylindrical coordinate mapping: (xFlat, yFlat) -> (X, Y, Z) |
| Animation | 40 VTU frames per mode covering one complete oscillation period |

---

## Output Fields

Each `.vtu` file (per mode, per time frame) contains:

| Field | Description | Notes |
|-------|-------------|-------|
| `DispMag` | Displacement magnitude at time t | m; red = max vibration |
| `AxialDisp` | Axial displacement component | m |
| `RadialDisp` | Radial displacement component | m |

Mesh node positions in each VTU are the **physically deformed 3D cylinder** — no Warp By Vector needed in ParaView.

---

## Repository Structure

```
shell_modal/
|
|-- shell_modal.edp               # Main FreeFEM++ simulation script
|-- README.md                     # This file
|
|-- D:\freefem++\shell_modal\     # Output directory (auto-created by script)
    |-- mode01.pvd                # Mode 1 animation (open this in ParaView)
    |-- mode01_t000.vtu           # Mode 1 frame 0
    |-- mode01_t001.vtu           # Mode 1 frame 1
    |-- ...
    |-- mode01_t039.vtu           # Mode 1 frame 39 (end of period)
    |-- mode02.pvd                # Mode 2 animation
    |-- ...
    |-- mode12.pvd                # Mode 12 animation
    |-- frequencies.txt           # All natural frequencies
```

---

## How to Run

### Requirements

- FreeFEM++ v4.10 or later: https://freefem.org
- ParaView v5.x or later: https://www.paraview.org
- Output directory `D:\freefem++\` must exist before running

### Step 1 - Run the simulation

```bash
FreeFem++ shell_modal.edp
```

The script will:
1. Assemble CLT A-matrix and print effective stiffness coefficients
2. Build the unrolled rectangular mesh (80x120 elements)
3. Assemble stiffness K and consistent mass M via varf
4. Run inverse power iteration with deflation for 12 modes
5. For each mode: write 40 VTU frames with 3D deformed cylinder coordinates
6. Save one PVD file per mode and a frequencies summary text file

Console output during eigenvalue extraction:
```
  Mode 1 converged at iter=23  lambda=3.29e+08
  Mode 2 converged at iter=31  lambda=1.22e+09
  ...
==============================================
 NATURAL FREQUENCIES
 Mode  f[Hz]   omega[rad/s]
   1    91.3    573.7
   2   176.0   1106.0
   ...
==============================================
```

### Step 2 - Open in ParaView

1. `File > Open` > navigate to `D:\freefem++\shell_modal\`
2. Change Files of type to `All Files (*.*)`
3. Select `mode01.pvd` > OK
4. Choose PVD Reader when prompted > OK
5. Click `Apply`
6. Set colour field to `DispMag`
7. Click `Rescale to Data Range Over All Timesteps`

### Step 3 - Visualize the vibrating cylinder

**Option A - Mode 1 breathing/bending**
```
Open mode01.pvd > Apply
Colour by DispMag
Press Play
-> 3D cylinder PHYSICALLY OSCILLATES (not a flat rectangle)
-> Red lobes pulse in and out around the circumference
```

**Option B - Higher circumferential modes**
```
Open mode03.pvd (f3 ~ 176 Hz, n=4 petals)
Open mode05.pvd (f5 ~ 192 Hz, n=3 lobes)
-> Each PVD shows a different cross-section deformation pattern
-> All are true dynamic oscillations, not static plots
```

**Option C - Axial mode comparison**
```
Open mode09.pvd (f9 ~ 301 Hz, mode (1;2) oval)
Open mode11.pvd (f11 ~ 313 Hz, mode (2;5) second axial)
-> mode11 shows TWO bands of lobes along the length
-> Visible as two rows of red/blue alternating regions
```

**Option D - Side-by-side comparison**
```
Layout > 2 columns
Left: mode01.pvd  Right: mode03.pvd
Both on Play simultaneously
-> Compare axial vs circumferential wave patterns side by side
```

Press `Play` in each PVD to watch the shell complete one full oscillation.

---

## What to Look for in Results

### 3D Cylinder Shape

Each VTU frame writes the cylindrical surface coordinates directly. The shell appears as a hollow 3D tube in ParaView — not the flat rectangle from the 2D mesh. The deformation is visible as radial bulging and contracting around the circumference.

### Mode Shape Patterns

Low-n modes (n=1,2) show large smooth lobes. High-n modes (n=5,6) show many small rapidly-alternating petals around the circumference. The axial index m controls how many half-waves appear along the tube length — m=2 gives two rows of lobes.

### Frequency Ordering

For thin CFRP shells, circumferential modes with intermediate n (typically n=3 to 5) have the lowest natural frequencies. Pure axial breathing (n=0) is much higher. This ordering matches the image from the reference paper.

### Convergence of Deflation

Each subsequent mode requires more deflation iterations to converge because the shift-invert scheme must be driven away from all previously found eigenvalues. Modes 1-3 typically converge in 20-30 iterations. Modes 10-12 may need 80-150 iterations.

---

## Failure Sequence (Mode Extraction Order)

```
Mode  1: Lowest omega^2 -> lowest f -> first mode extracted by power iteration
Mode  2: Second lowest, M-orthogonal to mode 1
...
Mode 12: Twelfth lowest, M-orthogonal to modes 1-11
         Most iterations needed due to repeated deflation
frequencies.txt: all 12 values sorted by f
```

---

## Common Errors and Fixes

| Error | Cause | Fix |
|-------|-------|-----|
| `the array size must be 2 not 1` | Vh2 is vector space, EigenValue API rejects Vh2[int] | Use inverse power iteration with real[int,int] modeVecs instead |
| `interpolation of no scalar FE component` | Assigning scalar to vector FE component inside loop | Fill DOF array directly: vv[2k]=sin(...) using Th(k).x |
| Flat rectangle in ParaView | VTU writes 2D coordinates (xFlat, yFlat) | Map to 3D: X=xFlat, Y=(R+dR)*sin(theta), Z=(R+dR)*cos(theta) |
| All modes same frequency | Deflation failed, iteration converges to mode 1 repeatedly | Check M-orthogonalisation loop; verify phiMax > 0 before normalising |
| Non-ASCII compile error | UTF-8 symbol in script | Use only 7-bit ASCII characters |
| `syntax error before token _` | Underscore in identifier | Use camelCase: `phiMax` not `phi_max` |
| Mode shapes look like noise | hmin too coarse (nx=80, ny=120 insufficient for n=6) | Increase ny to 160 for modes with n > 5 |

---

## Extending the Model

| Extension | What to change |
|-----------|----------------|
| Different laminate | Update ply angles and thicknesses in CLT assembly |
| Simply supported both ends | Replace `on(4, ux=uy=0)` with `on(4, uy=0)` and add `on(2, uy=0)` |
| Free-free boundary | Remove all Dirichlet BCs; rigid body modes appear at omega~0 |
| Static preload influence | Add axial force Nx to stiffness: K_prestress = K + Nx*Kg |
| More modes | Increase nev from 12 to 20; increase iteration limit |
| Different CFRP material | Update E1, E2, G12, nu12, rhoP |
| Glass-epoxy shell | E1=40 GPa, E2=8.3 GPa, G12=4.1 GPa, rhoP=2100 kg/m^3 |
| Thicker shell | Increase hShell; add bending D-matrix to the CLT formulation |

---

## Citation

If you use this code in your research, please cite:

```bibtex
@software{mishra_2026_shellmodal,
  author    = {Mishra, A.},
  title     = {Modal Analysis of Composite Cylindrical Shell - FreeFEM++},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20030608},
  url       = {https://doi.org/10.5281/zenodo.20030608}
}
```

Plain text citation:

> Mishra, A. (2026). *Modal Analysis of Composite Cylindrical Shell - FreeFEM++*. Zenodo. https://doi.org/10.5281/zenodo.20030608

---

## Author

**akshansh11**
GitHub: https://github.com/akshansh11

---

## License

<p>
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
<img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png"/>
</a>
<br/>
This work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
</p>

You are free to:

- **Share** - copy and redistribute the material in any medium or format
- **Adapt** - remix, transform, and build upon the material

Under the following terms:

- **Attribution** - You must give appropriate credit to akshansh11 and provide a link to this repository
- **NonCommercial** - You may not use the material for commercial purposes

Copyright 2026 akshansh11. All rights reserved for commercial use.
