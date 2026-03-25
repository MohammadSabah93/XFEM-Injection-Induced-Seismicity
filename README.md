## XFEM Injection-Induced Seismicity Simulator

XFEM Injection-Induced Seismicity Simulator
Overview

This repository provides a MATLAB-based, fully coupled poromechanical extended finite element method (XFEM) framework for modeling injection-induced seismicity in fractured reservoirs. The numerical formulation captures fracture slip governed by rate-and-state friction, fluid pressure evolution within fractures, and hydro-mechanical coupling with the surrounding porous medium.

The simulator is designed to study the transition between aseismic and seismic slip under fluid injection, with particular emphasis on fault reactivation and induced seismic events.

Key Features:

2D poroelastic domain using standard finite elements, 
XFEM enrichment (Heaviside and crack-tip functions) for fracture representation, 
Fully coupled hydro-mechanical formulation
Rate-and-state friction law for fault slip
Dynamic contact constraints along fracture interfaces
Newton–Raphson iterative solution scheme
Adaptive time stepping based on slip-velocity criteria
Fluid exchange between fracture and surrounding matrix
Computation of seismicity-related metrics (moment magnitude, stress drop, slip distribution, etc.)
Contact Algorithms

The codebase includes two alternative formulations for enforcing normal contact constraints along fracture interfaces:

🔹 Penalty Method

The penalty-based formulation enforces contact constraints approximately by introducing a penalty stiffness that penalizes interpenetration between fracture surfaces. This approach is computationally efficient and avoids introducing additional degrees of freedom.

🔹 Lagrange Multiplier Method

The Lagrange multiplier formulation enforces contact constraints exactly by introducing additional unknowns representing contact tractions. This method ensures strict satisfaction of the impenetrability condition, at the cost of increased system size.


