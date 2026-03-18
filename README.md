# Spatial PDE Modelling of α-Synuclein Aggregation in Parkinson's Disease

**BENG 227 — Computational Neuroscience | UC San Diego**  
**Group 8: Nitin Shreyes Venkatesan, Bivas Talukdar**

---

## Overview

This project extends the Yang et al. (2023) well-mixed ODE model of
α-synuclein (aSyn) aggregation and clearance in Parkinson's disease (PD)
into a **1D spatial reaction-advection-diffusion PDE system** along a
dopaminergic axon. The central hypothesis is that the selective
vulnerability of dopaminergic neurons — which have the longest unmyelinated
axons in the brain — emerges from the **spatial geometry of the axon**
combined with **length-dependent oxidative burden**, rather than from
kinetic parameters alone.

The model tracks six species along the axon axis x ∈ [0, L], with x = 0
at the soma and x = L at the synaptic terminal:

| Species | Role |
|---------|------|
| ROS | Reactive oxygen species — oxidative stress driver |
| aSyn\* | Aggregated α-synuclein — immobile, synapse-local |
| ERS | ER stress signal — diffuses from synapse |
| mTOR | Autophagy suppressor — fast retrograde transport |
| Beclin1 | Autophagy initiator — slow anterograde transport |
| Caspases | Apoptosis effector — anterograde transport |

---

## Key Scientific Results

### 1. Tristability Validation
The spatial PDE under homogeneous initial conditions exactly recovers
Yang's three stable attractors (Healthy, Intermediate, Disease) to
machine precision (error < 10⁻¹⁵ µM). Bifurcation diagrams across S₁,
S₂, and S₃ reproduce Yang et al. Figure 2 quantitatively.

```
State          aSyn* (µM)   Beclin1 (µM)   Caspases (µM)
─────────────────────────────────────────────────────────
Healthy          1.815         0.762           0.118
Intermediate     3.708         0.356           0.924
Disease         16.483         0.434           0.939
```

### 2. Numerical Stability Analysis
Von Neumann analysis on the linearized PDE system reveals a **stiffness
ratio of 3.7 × 10⁵** controlled by ROS diffusion (D = 6.5 × 10⁶ µm²/h).
Explicit solvers require ~13 million steps per 100h simulation — proven
impractical. The coupled spectral radius analysis confirms coupling adds
negligible instability at stable attractors (Δρ < 10⁻⁵), while the
saddle-point coupling tightens the explicit Δt bound by 741×.

### 3. Soma-Only Beclin1 Failure Mode Analysis
The originally proposed soma-only Beclin1 production was shown to be
**structurally degenerate** through four compounding failure modes:

- **FM1 — Analytical:** Effective decay length L_eff = 251 µm << L = 1000 µm → B1(synapse) ≈ 0 by construction
- **FM2 — Numerical:** All ICs converge to identical B1 = 6.8 × 10⁻⁵ µM at synapse regardless of starting point (even with full pool at synapse)
- **FM3 — Simulation:** Aagg grows linearly (R² = 0.9998 vs zero-clearance prediction) — no bistability possible
- **FM4 — Artifact:** δ(x−L) = 1/Δx creates 7.5× aggregation rate artifact across axon lengths

**Fix:** Distributed pool-bounded activation (Yang's ODE formulation applied pointwise).

### 4. Spatial Disease Initiation
Under asymmetric ICs, the spatial model produces a synaptic aSyn*
concentration ~1.9× higher than the soma at steady state. The disease
risk index [aSyn*]/[Beclin1] is highest at the synaptic terminal —
a result impossible to obtain from the well-mixed ODE.

### 5. Fixed-S₁ Axon Length Experiment
Under identical kinetic parameters, longer axons show **paradoxical
synaptic protection** through Caspase pool depletion:

```
L (µm)   aSyn*(syn)   Beclin1(syn)   Risk index
────────────────────────────────────────────────
200        2.984         0.775          3.85
1000       2.866         0.820          3.50
1500       2.807         0.845          3.32
```

The Caspase pool decays as C_T(x) = e^(−x/500), depleting apoptotic
signalling at the synapse of longer axons and releasing Beclin1
inactivation — a spatially emergent mechanism.

### 6. S₁(L) Model — Length-Dependent Vulnerability
The length-dependent oxidative stress model S₁(L) = S₁_base + α·L
(α = 0.001 µm⁻¹) reveals discrete regime transitions with axon length:

| Critical Length | Transition |
|----------------|------------|
| L < 900 µm | Bistable only — no intermediate state accessible |
| 900–1660 µm | Tristable — fate IC-dependent (therapeutic window) |
| 1660–2040 µm | Bistable I+D — healthy attractor weakening |
| L > 2040 µm | Monostable disease — healthy attractor gone |

The tristable window (900–1660 µm) coincides with the typical range of
human nigrostriatal dopaminergic axons.



## Installation

```bash
git clone https://github.com/yourusername/asyn-spatial-pde.git
cd asyn-spatial-pde

# Create environment
conda create -n beng227 python=3.10
conda activate beng227

# Install dependencies
pip install numpy scipy matplotlib pandas jupyter
```

### Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | ≥ 1.24 | Array operations |
| scipy | ≥ 1.10 | ODE/PDE solver (Radau), fsolve |
| matplotlib | ≥ 3.7 | All visualisations |
| pandas | ≥ 2.0 | Summary tables |
| jupyter | ≥ 1.0 | Notebook execution |


## Model Equations

The spatial PDE system (Yang kinetics + transport):

$$\frac{\partial R}{\partial t} = D_R \nabla^2 R + k_1\!\left(1 + S_1 + d_\alpha \frac{(A/k_\alpha)^n}{1+(A/k_\alpha)^n}\right) - k_2 R S_2$$

$$\frac{\partial A}{\partial t} = k_3 R S_3 - k_4 k_5 A B M$$

$$\frac{\partial E}{\partial t} = D_E \nabla^2 E + k_6 k_7 A(E_T - E) - k_8 E$$

$$\frac{\partial M}{\partial t} = D_M \nabla^2 M + v_\ell \frac{\partial M}{\partial x} + (k_9 + k_{10}E)(M_T - M) - (k_{11} + k_{12}B)M$$

$$\frac{\partial B}{\partial t} = D_B \nabla^2 B - v_B \frac{\partial B}{\partial x} + \frac{(k_{13}+k_{14}E)(B_T-B)}{J_{be}+B_T-B} - \frac{(k_{15}+k_{16}C+k_{17}M)B}{J_{be}+B}$$

$$\frac{\partial C}{\partial t} = D_C \nabla^2 C - v_C \frac{\partial C}{\partial x} + \frac{(k_{18}+k_{19}E+k_{20}M)(C_T-C)}{J_{ca}+C_T-C} - \frac{(k_{21}+k_{22}B)C}{J_{ca}+C}$$

### Transport Parameters

| Species | D (µm²/h) | v (µm/h) | Direction |
|---------|-----------|----------|-----------|
| ROS | 6.5 × 10⁶ | 0 | Diffusion only |
| aSyn\* | 0 | 0 | Immobile |
| ERS | 3.6 × 10⁴ | 0 | Diffusion only |
| mTOR | 36 | 7200 | Retrograde (−x) |
| Beclin1 | 5.4 × 10⁴ | 36 | Anterograde (+x) |
| Caspases | 7.2 × 10⁴ | 36 | Anterograde (+x) |

### Numerical Scheme

| Component | Method |
|-----------|--------|
| Diffusion | Second-order central differences |
| Advection | First-order upwind |
| Time integration | Method of Lines → scipy Radau (L-stable implicit) |
| Boundary conditions | Zero-flux Neumann at x=0 and x=L |
| Grid | Nx = 100, Δx = L/100 µm |
| Tolerances | rtol = 1e-5, atol = 1e-7 |

---

## Numerical Stability Summary

```
Species       D (µm²/h)    v (µm/h)    Δt_max (Nx=100)    Controls?
──────────────────────────────────────────────────────────────────────
ROS           6.50e+06        0.0        7.69e-06 h          ← YES
Asol          1.44e+05       36.0        3.47e-04 h
ERS           3.60e+04        0.0        1.39e-03 h
mTOR          3.60e+01     7200.0        1.39e-03 h  (CFL)
Beclin1       5.40e+04       36.0        9.26e-04 h
Caspases      7.20e+04       36.0        6.94e-04 h

Stiffness ratio : 3.7 × 10⁵
Steps for 100h  : ~13,000,000  (explicit — IMPRACTICAL)
Solver verdict  : Radau mandatory
```

---

## State Classification

States are classified at each spatial grid point by normalized
Euclidean distance to the Yang ODE reference attractors in
[aSyn\*, Beclin1, Caspases] space:

```python
# Classification uses S1=1.462 reference FPs regardless of simulation S1
stable_fps = find_stable_fps({**CFG, 'S1': 1.462})
# Returns: Healthy (1.815 µM), Intermediate (3.708 µM), Disease (16.483 µM)
```

---

## Key Design Decisions

**Why distributed Beclin1 activation?**  
Soma-only production predicts B1(synapse) ≈ 0 analytically (L_eff = 251 µm << 1000 µm), confirmed numerically to be IC-independent. This structurally disables autophagy at the synapse by construction, making the model degenerate. Distributed pool-bounded activation (Yang's ODE formulation applied pointwise) restores biologically meaningful clearance gradients.

**Why Radau?**  
Stiffness ratio 3.7 × 10⁵ at Nx=100 (scales as Δx⁻²). Explicit solvers require ~13M steps per 100h. Radau is L-stable and handles all stiffness ratios unconditionally.

**Why S₁(L) = S₁_base + α·L?**  
Longer axons contain proportionally more mitochondria → more total ROS production. Fixed-S₁ experiments showed geometry alone protects longer axons (Caspase depletion); vulnerability requires the additional oxidative burden from mitochondrial length scaling.

---

## Biological Findings

1. **The synapse is structurally the most vulnerable location** — highest ROS (mitochondrial source at terminal), lowest Beclin1 clearance capacity (distributed activation cannot fully compensate at the distal end under high stress), and the site of aSyn* boundary accumulation.

2. **The intermediate state is clinically distinct from disease** — Caspases are near-maximal (0.924 µM) at the intermediate attractor, meaning apoptotic commitment precedes overwhelming aggregate accumulation. The intermediate state is where therapeutic intervention is most relevant.

3. **Axon length creates paradoxical protection at fixed stress** — Caspase pool depletion along long axons reduces Beclin1 inactivation, lowering the disease risk index by 14% from L=200 to L=1500 µm. This result is invisible to the ODE model.

4. **Vulnerability arises from length-stress coupling** — The S₁(L) model predicts a critical axon length of ~900 µm below which the intermediate state is inaccessible. Above ~2040 µm the healthy attractor disappears entirely. The tristable window (900–1660 µm) coincides with human nigrostriatal axon lengths.

5. **S₂ (antioxidant defense) has the sharpest protective threshold** — Tristable window spans only ΔS₂ = 0.17 units, predicting abrupt disease onset from gradual age-related antioxidant decline.



---

## License

This project is for academic purposes (BENG 227, UC San Diego).
All kinetic parameters are reproduced from Yang et al. (2023) under
academic fair use.

---

*For questions, contact: Nitin Shreyes Venkatesan (UC San Diego), Bivas Talukdar (UC San Diego)*
