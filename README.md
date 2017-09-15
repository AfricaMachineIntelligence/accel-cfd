# Performance Optimization of Numerical Solution Algorithm for Microchannel Flow Boiling

## Abstract

In recent years, there have been significant advancements to computational hardware that allow for massive parallelism in computational fluid dynamics algorithms. However, as academia and industry increase the complexity of these algorithms computational runtimes can extend to days and weeks with state of the art approaches. In this study, convolutional neural networks are used as a novel approach to preconditioning computational fluid dynamics simulations in order to significantly reduce the number of iterations required for convergence. A parallel enabled, computational fluid dynamics simulation is created to match state of the art performance for a microchannel, two-phase, convection problem. A convolutional neural network is then trained on simulation data from the simulation using a regression approximation to estimate the effects of
individual forces and constraints from the Navier-Stokes equations. The trained network is then implemented as a feed-forward network to create an initial guess for the pressure, velocity, and temperature fields used in the SIMPLER algorithm. This hybrid algorithm is then compared with state of the art serial and parallel implantations of the microchannel, two-phase, convection simulation.

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

## Directory Map

TBC.
