# Skyrme-model
In this repository you will find codes which are related to the 3-dimensional Skyrme model. There are simple and introductory codes to the Skyrme model as well as more advanced calculations which are focused in the study of infinite nuclear matter. Concretely:

B1.py solves the B = 1 skyrmion in spherical coordinates using a shooting method. This is the simplest code and the first which one should be able to write when working with the Skyrme model. In this code the generalized Skyrme lagrangian may be considered.

Spectral_Skyrme.py solves the same differential equation as the first one but using spectral methods. In this concrete case, the shooting method is a simpler approach, however this code is an excellent example for solving a multidomain non-linear unidimensional differential equation. This file needs the presence of the Polynomials.py code to work.

Skyrmions_Plotter.py represents any 3-dimensional skyrmion if the energy density and the pion fields are given. It uses the Python library Mayavi to plot the skyrmion and considers the Leeds colouring convention.

Gradient_Flow.py uses a more advanced procedure to solve a generic skyrmion of B >= 1. It constructs an initial configuration just taking the topological degree of the desired skyrmion and uses a gradient flow method to minimize the energy functional and find the solution.

Skyrme.ipynb is a jupyter notebook which performs some useful symbolic calculations of the Skyrme model. In it you may find the expressions of the Skyrme lagrangian as well as the Isospin inertia tensor for the hedgehog ansatz, but also for the generic SU(2) expansion of the Skyrme field.

RatMap_constructor.ipynb is a very useful jupyter notebook which specifically combined with the Gradient_Flow.py code might be of extreme interest to construct B >= 1 skyrmions. It takes the rational map ansätze and computes the symbolic expressions.
