# Performance Optimization of Numerical Solution Algorithm for Microchannel Flow Boiling

## Overview


## Objective

The objective of this study [thesis] is to optimize the performance of a numerical solution algorithm for micro-channel flow boiling problems (code used in previous work). The performance indices for this study are (1) solution accuracy, (2) computational duration [runtime], and (3) stability of solution algorithm [does the algorithm find convergence]. The study involves the following steps:
- Identify the bottlenecks in the existing microchannel flow boiling solution algorithm that significantly impact computational runtime.
  - Identify bottlenecks in literature.
  - Identify bottlenecks in actual code.
- Identify potential improvements to significant bottlenecks from literature review.
- Perform a comparative / parametric study of potential improvements (linear solvers) on a simplified problem (3D diffusion).
  - Perform a comparative analysis of the different improvement methods.
  - Test difference linear system solution methods on a 3D diffusion problem (parametric study).
- Determine the cost of implementation of each potential improvement method.
  - How big of a change is it to implement new methods in CFD?
- Implement (? #) of potential improvements in the microchannel flow boiling solution algorithm.
- Perform an analysis of the impact made by improving bottleneck areas.
- Identify additional areas related to this study that can be investigated in future research.

## Checklist (TODO)
- Linear Solvers
  - TDMA
  - BiCG
  - GMRES
  - BiCGStab
  - Pardiso
- Preconditioners
  - Jacobi
  - SSOR
  - ILU
  - Multigrid
- Tools (Vector / Restructuring)
  - Index Calculation
  - Variable Number
  - Sparse Matrix Creation
- Convection Integration / Testing
  - Integration with Existing Code
  - Testing with Existing Code
- Batch Scripts
  - Diffusion Tests
  - Convection Tests
  
## Directory Map
