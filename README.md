## XFEM Injection-Induced Seismicity Simulator

This repository provides a MATLAB-based, fully coupled poromechanical XFEM framework developed to model injection-induced seismicity in fractured reservoirs. The numerical formulation accounts for fracture slip governed by rate-and-state friction, fluid pressure evolution within the fracture, and dynamic contact constraints via Lagrange multipliers. Adaptive time stepping is used to accurately track the transition between aseismic and seismic slip.

## Features
- 2D poroelastic domain using standard finite elements
- Heaviside and crack tip enrichment for intersected elements
- Rate-and-state friction behavior for seismic slip
- Dynamic fully coupled hydro-mechanical solution using Newton–Raphson iterations
- Adaptive time stepping based on slip-velocity criteria
- Contact mechanics enforced through Lagrange multipliers
- Fluid exchange between fracture and surrounding matrix
- Computation of seismicity metrics (Mw, stress drop, slip distribution, etc.)
- Stress recovery and visualization tools included
