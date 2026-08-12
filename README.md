# XFEM Injection-Induced Seismicity — Lagrange-Multiplier Formulation

[![MATLAB](https://img.shields.io/badge/MATLAB-source%20code-e16737)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This branch contains the stabilized Lagrange-multiplier implementation of a coupled poromechanical XFEM framework for injection-induced fault reactivation.

## Features

- two-dimensional poroelastic deformation and fluid flow;
- Heaviside and crack-tip XFEM enrichment;
- fracture–matrix hydraulic exchange;
- rate-and-state friction;
- exact normal-contact enforcement through Lagrange multipliers;
- Newton–Raphson nonlinear iterations;
- inertia, dynamic damping, and adaptive time stepping; and
- seismicity measures and visualization utilities.

## Run the example

1. Open this branch in MATLAB with all `.m` files on the path.
2. Review `defineModelParameters.m` and `defineBoundaryConditions.m`.
3. Run:

```matlab
X_FEM_PoroElastic_Lagrange
```

The principal contact routines are `Lagrange.m`, `StabLagrange.m`, and `StiffnessInterface_Lagrange.m`.

## Alternative formulation

The penalty-contact implementation is available in the [`Contact_Penalty`](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/tree/Contact_Penalty) branch.

## Research-software notice

This is research code for numerical experimentation. Verify mesh and time-step convergence and calibrate model parameters before interpreting results physically.

## Citation

See the project’s [citation metadata](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/blob/main/CITATION.cff) and [main documentation](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity).

## License

Released under the [MIT License](LICENSE).
