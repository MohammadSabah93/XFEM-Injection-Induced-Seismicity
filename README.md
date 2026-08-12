# XFEM Injection-Induced Seismicity — Penalty Formulation

[![MATLAB](https://img.shields.io/badge/MATLAB-source%20code-e16737)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This branch contains the penalty-contact implementation of a coupled poromechanical XFEM framework for injection-induced fault reactivation.

## Features

- two-dimensional poroelastic deformation and fluid flow;
- Heaviside and crack-tip XFEM enrichment;
- fracture–matrix hydraulic exchange;
- rate-and-state friction;
- computationally efficient normal-contact enforcement through a penalty method;
- Newton–Raphson nonlinear iterations;
- inertia, dynamic damping, and adaptive time stepping; and
- seismicity measures and visualization utilities.

## Run the example

1. Open this branch in MATLAB with all `.m` files on the path.
2. Review `defineModelParameters.m` and `defineBoundaryConditions.m`.
3. Run:

```matlab
X_FEM_PoroElastic_V5_penalty
```

The principal contact routine is `StiffnessInterface_Penalty.m`.

## Alternative formulation

The stabilized Lagrange-multiplier implementation is available in the [`Contact_Lagrange-Multiplier`](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/tree/Contact_Lagrange-Multiplier) branch.

## Research-software notice

This is research code for numerical experimentation. Verify mesh and time-step convergence and calibrate model parameters before interpreting results physically.

## Citation

See the project’s [citation metadata](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/blob/main/CITATION.cff) and [main documentation](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity).

## License

Released under the [MIT License](LICENSE).
