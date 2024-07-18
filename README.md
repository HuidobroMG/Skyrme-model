# README.md
In this repository you will find multiple codes related to the (3-dimensional) nuclear Skyrme model. We first include some Jupyter notebooks in which a systematic framework of symbolic computations is developed in the Skyrme model for isolated and periodic solutions. These notebooks are helpful to understand some complicated mathematical concepts and get familiar with the model. The C++ codes are used to obtain the solutions due to the computational difficulties that the model presents and the faster performance of the programming language. The Python codes are used to extract the nuclear properties from the solutions previously obtained and the comparison with physical observables, specifically, isolated nuclei and infinite nuclear matter. These codes have been developed and used for the publication of my PhD thesis, which may be found in the following [link](https://minerva.usc.es/xmlui/handle/10347/32925). Additionally, my PhD thesis and the thesis defence are also included in this repository in .pdf format for a complete description of the Skyrme model.

B1.py solves the B = 1 skyrmion in spherical coordinates using a shooting method. This is the simplest code and the first which one should be able to write when working with the Skyrme model. In this code the generalized Skyrme lagrangian may be considered.

Spectral_Skyrme.py solves the same differential equation as the first one but using spectral methods. In this concrete case, the shooting method is a simpler approach, however this code is an excellent example for solving a multidomain non-linear unidimensional differential equation. This file needs the presence of the Polynomials.py code to work.

Skyrmions_Plotter.py represents any 3-dimensional skyrmion if the energy density and the pion fields are given. It uses the Python library Mayavi to plot the skyrmion and considers the Leeds colouring convention.

Gradient_Flow.py uses a more advanced procedure to solve a generic skyrmion of B >= 1. It constructs an initial configuration just taking the topological degree of the desired skyrmion and uses a gradient flow method to minimize the energy functional and find the solution.

Skyrme.ipynb is a jupyter notebook which performs some useful symbolic calculations of the Skyrme model. In it you may find the expressions of the Skyrme lagrangian as well as the Isospin inertia tensor for the hedgehog ansatz, but also for the generic SU(2) expansion of the Skyrme field.

RatMap_constructor.ipynb is a very useful jupyter notebook which specifically combined with the Gradient_Flow.py code might be of extreme interest to construct B >= 1 skyrmions. It takes the rational map ansätze and computes the symbolic expressions.

Skyrmions.cc is a C++ code which constructs (using the rational map approximation) a skyrmions with a given baryon number B = 1-8 and 4N^3 (where N is an integer number) and finds the minimal energy configuration performing an accelerated gradient flow method.

Skyrmions_Parallel.cc is the parallelized version of Skyrmions.cc using the OpenMP library of the C++ language. In order to compile this code, one must run the following command: g++ -fopenmp Skyrmions_Parallel.cc
