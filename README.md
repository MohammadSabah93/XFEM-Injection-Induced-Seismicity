# XFEM Simulator for Injection-Induced Seismicity

[![MATLAB](https://img.shields.io/badge/MATLAB-research%20code-e16737)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Article DOI](https://img.shields.io/badge/Article-10.1016%2Fj.compgeo.2025.107803-blue)](https://doi.org/10.1016/j.compgeo.2025.107803)
[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--7260--6138-A6CE39)](https://orcid.org/0000-0002-7260-6138)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Mohammad%20Sabah-0A66C2)](https://www.linkedin.com/in/mohammad-sabah)

A MATLAB research framework developed by **Mohammad Sabah** for computational geomechanics and numerical simulation of **injection-induced seismicity**, combining coupled poromechanics, the **extended finite element method (XFEM)**, nonlinear fault contact, **rate-and-state friction**, inertia, and dynamic rupture.

This repository accompanies the published *Computers and Geotechnics* formulation by Sabah et al. (2026) and provides two alternative contact implementations for reproducible numerical experimentation. The work is relevant to induced seismicity, fault reactivation, reservoir geomechanics, geothermal stimulation, coupled hydromechanical modeling, and earthquake-rupture simulation.

## Repository structure

The source code is maintained in two implementation branches:

| Branch | Contact formulation | Main driver |
|---|---|---|
| [`Contact_Lagrange-Multiplier`](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/tree/Contact_Lagrange-Multiplier) | Stabilized Lagrange multiplier | `X_FEM_PoroElastic_Lagrange.m` |
| [`Contact_Penalty`](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/tree/Contact_Penalty) | Penalty method | `X_FEM_PoroElastic_V5_penalty.m` |

The `main` branch is intentionally used as the documentation, citation, and project landing page.

## Numerical capabilities

- two-dimensional coupled poroelastic deformation and fluid flow;
- XFEM representation of an embedded fracture/fault;
- matrix–fracture hydraulic exchange;
- rate-and-state friction;
- nonlinear fault contact using penalty or stabilized Lagrange-multiplier formulations;
- inertia and dynamic boundary damping;
- Newton–Raphson nonlinear solution;
- adaptive time stepping across aseismic and seismic slip; and
- seismicity measures including slip, stress drop, seismic moment, and moment magnitude.

## Quick start

### Requirements

- MATLAB;
- all `.m` files from the selected implementation branch available on the MATLAB path.

A minimum MATLAB release and toolbox compatibility matrix have not yet been formally established.

### Clone the repository

```bash
git clone https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity.git
cd XFEM-Injection-Induced-Seismicity
```

### Stabilized Lagrange-multiplier formulation

```bash
git checkout Contact_Lagrange-Multiplier
```

Run in MATLAB:

```matlab
X_FEM_PoroElastic_Lagrange
```

### Penalty formulation

```bash
git checkout Contact_Penalty
```

Run in MATLAB:

```matlab
X_FEM_PoroElastic_V5_penalty
```

Before running a case, review geometry and time-stepping controls in the selected driver, material and frictional properties in `defineModelParameters.m`, and boundary conditions in `defineBoundaryConditions.m`.

## Recommended reproducibility checks

For quantitative interpretation, document and test at least:

- mesh resolution;
- time-step limits and adaptive stepping parameters;
- nonlinear convergence tolerance;
- rate-and-state friction parameters;
- fault-contact parameters;
- damping parameters; and
- mechanical and hydraulic boundary conditions.

Mesh- and time-step-convergence studies should be performed before comparing physical outcomes across cases.

## Research-software status

This code is intended for scientific development, verification, and numerical experimentation. It is **not** an operational seismic-hazard forecasting tool. Users should independently verify the implementation and calibrate model parameters before site-specific interpretation.

## Publication

If this software supports your research, please cite:

> Sabah, M., Hofmann, H., Cacace, M., Jalali, M. R., & Kivi, I. R. (2026). Modeling injection-induced seismicity using a fully coupled poroviscoelasto-dynamic extended finite element approach with stabilized contact and rate-and-state friction. *Computers and Geotechnics, 191*, 107803. https://doi.org/10.1016/j.compgeo.2025.107803

Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## Related project

For the hybrid implicit–explicit time-integration implementation, see:

[`Hybrid-implicit-explicit-XFEM-simulation-of-injection-induced-seismicity`](https://github.com/MohammadSabah93/Hybrid-implicit-explicit-XFEM-simulation-of-injection-induced-seismicity)

## License

Released under the [MIT License](LICENSE).

## Author

**Mohammad Sabah, PhD**  
Computational geomechanics · induced seismicity · coupled multiphysics · XFEM · rate-and-state friction  
Technische Universität Berlin  
[ORCID](https://orcid.org/0000-0002-7260-6138) · [LinkedIn](https://www.linkedin.com/in/mohammad-sabah) · [GitHub](https://github.com/MohammadSabah93)

Questions, reproducibility requests, and bug reports are welcome through [GitHub Issues](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/issues).
