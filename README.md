# Performance Optimization of Numerical Solution Algorithm for Microchannel Flow Boiling

## Abstract

In recent years, there have been significant advancements to computational hardware that allow for massive parallelism in computational fluid dynamics algorithms. However, as academia and industry increase the complexity of these algorithms computational runtimes can extend to days and weeks with state of the art approaches. In this study, convolutional neural networks are used as a novel approach to preconditioning computational fluid dynamics simulations in order to significantly reduce the number of iterations required for convergence. A parallel enabled, computational fluid dynamics simulation is created to match state of the art performance for a microchannel, two-phase, convection problem. A convolutional neural network is then trained on simulation data from the simulation using a regression approximation to estimate the effects of
individual forces and constraints from the Navier-Stokes equations. The trained network is then implemented as a feed-forward network to create an initial guess for the pressure, velocity, and temperature fields used in the SIMPLER algorithm. This hybrid algorithm is then compared with state of the art serial and parallel implantations of the microchannel, two-phase, convection simulation.

## Objective

Typical computational fluid dynamics simulations require a significant amount of computational resources to run in a reasonable amount of time. As the complexity of these simulations increase, however, it is common for the computational duration of a solution to take days or weeks. While computational hardware increases in ability at an exponential rate, researchers in academia and industry often find the bottleneck in the design process to be running multiple CFD simulations. 

One potential improvement to reduce the computational duration of a CFD simulation is to utilize a highly parallel solution algorithm on parallel computation hardware, like a graphics processing unit (GPU). A number of recent studies have shown this to be a highly-successful method for reducing the computational duration of serial CFD codes by factors of up to 10 times. For a serial simulation that takes approximately 140 hours, a parallel simulation could achieve results of 14 hours. While this is a significant improvement, 14 hours is still a long time to wait.

The biggest contributor to the long duration for CFD simulations is the time spent in the iterative solution process. A number of hybrid solution algorithms have been developed in recent years to improve the performance of convergence of these solvers, but a close approximate precondition method (the initial guess) is by far the most significant aspect of time spent in iterations. 

A sub-field of Artificial Intelligence called Machine Learning may offer a potential for significantly improving the computational duration of CFD simulations. A technique in Machine Learning, called Deep Learning, utilizes multi-layer convolutional neural networks to classify objects based on image inputs. This technique of using CNNs can be applied to CFD simulations as well. Instead of classifying the CFD simulation, a regression model can be used to estimate the parameters of pressure, velocity, and temperature given boundary conditions, initial conditions, and geometry. By developing and training a CNN model to generate the initial guess for the SIMPLER algorithm, a parallel enabled CFD simulation can outperform a serial CFD simulation by over 30 times. Therefore, the 140 hour simulation can be reduced to less than 5 hours.

I propose to develop a CFD simulation algorithm that solves the microchannel flow boiling problem presented in my previous work. This algorithm will be written to utilize highly-parallel, computational architectures (GPUs) and will incorporate a trained feed forward neural network to precondition the SIMPLER solution algorithm. The performance of this algorithm will them be directly compared with the performance of the previous, serial implementation used during previous work. 

## Directory Map

TBC.
