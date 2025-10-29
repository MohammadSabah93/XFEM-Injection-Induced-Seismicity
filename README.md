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

## Code Overview
This MATLAB script performs the full hydro-mechanical simulation of injection-induced seismicity using an XFEM formulation. Below is a summary of the major stages in the code:
  
| Section                   | Description                                                                                                                                 |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| **Mesh generation**       | Creates a structured 2D finite element mesh and defines the fracture using level-set enrichment.                                            |
| **Fault geometry**        | Computes unit normal and tangential vectors along the fracture for interface calculations.                                                  |
| **Model parameters**      | Loads material, hydraulic, and frictional properties from `defineModelParameters()`.                                                        |
| **Boundary conditions**   | Applies mechanical and fluid constraints using `defineBoundaryConditions()`.                                                                |
| **Initialization**        | Allocates solution vectors and sets initial values for displacement, pressure, aperture, slip state, etc.                                   |
| **Matrix assembly**       | Assembles global stiffness, mass, damping, storage, conductance, and coupling matrices, including interface terms and Lagrange multipliers. |
| **Nonlinear solver**      | Uses Newton-Raphson iterations and adaptive time stepping to solve the coupled system at each time step.                                    |
| **Seismicity parameters** | Computes moment magnitude, stress drop, and slip metrics using `SeismicityParameters()`.                                                    |
| **Stress recovery**       | Calculates element stress components and principal values via `elemStress2()`.                                                              |
| **Post-processing**       | Visualizes displacement, pressure, and stress fields using `postProcess()`.                                                                 |

| Function                                  | Purpose                                                             |
| ----------------------------------------- | ------------------------------------------------------------------- |
| `MeshGeneration_2D()`                     | Generates nodes, connectivity, and domain discretization.           |
| `levelSet()`                              | Detects enriched nodes and interface elements for XFEM formulation. |
| `calcDOF()`                               | Computes total DOFs for displacement and pore pressure.             |
| `defineModelParameters()`                 | Defines poromechanical and friction law parameters.                 |
| `defineBoundaryConditions()`              | Specifies fixed DOFs, injection rates, and surface tractions.       |
| `MassMatrix()` / `StiffnessMatrix()`      | Builds mechanical mass and stiffness matrices.                      |
| `ConductanceMatrix()` / `StorageMatrix()` | Assembles hydraulic conductivity and storage matrices.              |
| `CouplingMatrix()`                        | Assembles Biot coupling terms between pressure and deformation.     |
| `CouplingInterface()`                     | Fluid-solid coupling along the fracture.                            |
| `Lagrange()` / `StabLagrange()`           | Constraint enforcement and stabilization for normal contact.        |
| `StiffnessInterface_Lagrange()`           | Interface constitutive tangent (rate-and-state friction & contact). |
| `FaultForce()`                            | Assembles fault traction contributions.                             |
| `interface_flow()`                        | Governs fluid flow and aperture evolution along the interface.      |
| `updateInterface()`                       | Updates slip, slip-velocity, fluid variables, and friction state.   |
| `SeismicityParameters()`                  | Computes seismic moment, Mw, static stress drop, etc.               |
| `elemStress2()`                           | Recovers stress tensors and principal directions.                   |
| `postProcess()`                           | Creates figures for pressure, stresses, and displacements.          |

## Notes
* The code uses adaptive time stepping based on slip velocity evolution for stability and efficiency.
* Normal contact is enforced using a Lagrange multiplier method with stabilization.
* Fracture behavior follows a rate-and-state friction model transitioning between stick and dynamic slip.
