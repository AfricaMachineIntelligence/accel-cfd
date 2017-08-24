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

### Framework

The computational framework for improving legacy CFD code is intended to be a general guideline for academic researchers to implement performance improvements into older CFD code. This is basically the first chapters / sections of the thesis where everything is explained. My idea is to provide a thorough description of how the linear system solution aspect of legacy CFD codes (similar to yours) can be improved by swapping the solution algorithms.

### Optimization Study

The optimization study consists of the comparison between different linear system solution algorithms, difference parallel protocols, and different hardware configurations. For this study, I plan to (1) minimize the computational duration (total time taken to solve a particular linear system) and (2) minimize the total iterations required to reach a particular error threshold. I plan on considering the following: 

- TDMA
- Biconjugate Gradients
- GMRES
- Biconjugate Gradients Stabilized

This objective is where the primary work comes through. Each of the solution algorithms will be compared with one another on different hardware configurations with different preconditioners. The parallel schemes will apply to particular solution methods based on the limitations of the algorithm. 

### Cost Determination

The cost determination of the solution algorithms considers the scalability (how well the code performs on more powerful machines), the complexity (how complicated the algorithm is to implement), and the stability of the solution method when different problem types are considered.

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
