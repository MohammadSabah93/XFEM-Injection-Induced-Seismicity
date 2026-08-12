# XFEM Simulator for Injection-Induced Seismicity

[![MATLAB](https://img.shields.io/badge/MATLAB-research%20code-e16737)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.compgeo.2025.107803-blue)](https://doi.org/10.1016/j.compgeo.2025.107803)

A MATLAB research framework for simulating injection-induced fault reactivation using a coupled poromechanical extended finite element method (XFEM), rate-and-state friction, dynamic rupture, and nonlinear fault contact.

## Choose a contact formulation

The source code is maintained in two implementation branches:

| Branch | Contact formulation | Main driver |
|---|---|---|
| [`Contact_Lagrange-Multiplier`](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/tree/Contact_Lagrange-Multiplier) | Stabilized Lagrange multiplier | `X_FEM_PoroElastic_Lagrange.m` |
| [`Contact_Penalty`](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/tree/Contact_Penalty) | Penalty method | `X_FEM_PoroElastic_V5_penalty.m` |

The `main` branch serves as the documentation and citation landing page.

## Capabilities

- two-dimensional coupled poroelastic deformation and fluid flow;
- XFEM representation of an embedded fracture;
- matrix–fracture hydraulic exchange;
- rate-and-state friction;
- alternative penalty and Lagrange-multiplier contact formulations;
- inertia and dynamic boundary damping;
- Newton–Raphson nonlinear solution;
- adaptive time stepping across aseismic and seismic slip; and
- seismicity measures including slip, stress drop, seismic moment, and moment magnitude.

## Quick start

### Requirements

- MATLAB;
- all `.m` files from the selected branch available on the MATLAB path.

A minimum MATLAB release and toolbox compatibility matrix have not yet been established.

### Run a formulation

```bash
git clone https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity.git
cd XFEM-Injection-Induced-Seismicity
```

For the Lagrange-multiplier version:

```bash
git checkout Contact_Lagrange-Multiplier
```

Then run `X_FEM_PoroElastic_Lagrange.m` in MATLAB.

For the penalty version:

```bash
git checkout Contact_Penalty
```

Then run `X_FEM_PoroElastic_V5_penalty.m` in MATLAB.

Before running, review the geometry and time-stepping settings in the selected driver, material and frictional properties in `defineModelParameters.m`, and boundary conditions in `defineBoundaryConditions.m`.

## Research-software status

This code is intended for scientific development and numerical experimentation, not operational seismic-hazard forecasting. Users should independently verify the implementation and conduct mesh, time-step, and parameter-sensitivity studies before physical interpretation.

## Reference

If this software supports your research, please cite:

> Sabah, M., Hofmann, H., Cacace, M., Jalali, M. R., & Kivi, I. R. (2026). Modeling injection-induced seismicity using a fully coupled poroviscoelasto-dynamic extended finite element approach with stabilized contact and rate-and-state friction. *Computers and Geotechnics, 191*, 107803. https://doi.org/10.1016/j.compgeo.2025.107803

Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## License

Released under the [MIT License](LICENSE).

## Contact

Questions, reproducibility requests, and bug reports are welcome through [GitHub Issues](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/issues).
