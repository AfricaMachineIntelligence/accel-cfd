# Performance Optimization in Computational Fluid Dynamics with Parallelism for Linear System Solvers

## Overview


## Objective

The scope of this study is to provide a framework for improving the computational performance of legacy CFD code, perform an optimization study of the current solutions to CFD problems, and determine the costs associated with implementing new technologies in legacy CFD code. The framework for improving the computational performance of legacy CFD code considers two classes of examples. The  first class of examples, considers simple diffusion problems solved using the finite volume method and the centralized differencing scheme. This class of example is primarily used due to its simplicity as implementing parallelism on a simple problem provide a higher adherence of understanding when compared to a problem with increased complexity. The second class of examples focused on a two-phase flow boiling problem that incorporates a solution to the full, incompressible Navier-Stokes equations. This example is used to demonstrate the effectiveness of implementing parallelism in a complex, CFD problem.

The optimization study considers the rate of convergence for the solutions to sparse linear systems, the computational resource requirement for different solutions to sparse linear systems, and the duration of computation for given solutions on a variety of hardware platforms. The study also considers different means of implementation for each solution algorithm and how that implementation affects the performance metrics of the optimization study. The optimization study considers the same two classes of CFD problems as mentioned previously and considers cases for a 1-dimensional, 2-dimensional, and 3-dimensional environment.

Finally, the determination of cost associated with implementing each particular CFD solution is examined for both the simple diffusion problem and the two-phase boiling problem. The performance indications of cost are considered for the 1-dimensional, 2-dimensional, and 3-dimensional environments and assess the required revision and restructuring of legacy code necessary to implement new solution types.

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
