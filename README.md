# XFEM Simulator for Injection-Induced Seismicity

[![MATLAB](https://img.shields.io/badge/MATLAB-research%20code-e16737)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.compgeo.2025.107803-blue)](https://doi.org/10.1016/j.compgeo.2025.107803)

A research framework for simulating injection-induced fault reactivation using a fully coupled poroviscoelastodynamic extended finite element method (XFEM).

## Scientific scope

The associated formulation resolves the coupled processes that connect reservoir pressurization to fault slip and dynamic rupture:

- two-dimensional poromechanics with inertial effects;
- XFEM representation of embedded fractures;
- fracture flow and matrix–fracture hydraulic exchange;
- rate-and-state friction;
- nonlinear normal contact;
- adaptive time stepping across reservoir and rupture timescales; and
- seismicity metrics including slip, stress drop, seismic moment, and moment magnitude.

Two contact formulations are considered: a computationally efficient penalty method and an exact Lagrange-multiplier method.

## Repository status

> **Documentation release.** The current `main` branch contains the project documentation and MIT license, but it does not yet contain the MATLAB source files or a reproducible example. Until those files are added, this repository should not be treated as a runnable software release.

For the available hybrid implicit–explicit MATLAB implementation, see [Hybrid IMEX XFEM for Injection-Induced Seismicity](https://github.com/MohammadSabah93/Hybrid-implicit-explicit-XFEM-simulation-of-injection-induced-seismicity).

## Reference

If this formulation supports your research, please cite:

> Sabah, M., Hofmann, H., Cacace, M., Jalali, M. R., & Kivi, I. R. (2026). Modeling injection-induced seismicity using a fully coupled poroviscoelasto-dynamic extended finite element approach with stabilized contact and rate-and-state friction. *Computers and Geotechnics, 191*, 107803. https://doi.org/10.1016/j.compgeo.2025.107803

Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## License

Released under the [MIT License](LICENSE).

## Contact

Questions, reproducibility requests, and bug reports are welcome through [GitHub Issues](https://github.com/MohammadSabah93/XFEM-Injection-Induced-Seismicity/issues).
