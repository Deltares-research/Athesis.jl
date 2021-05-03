# Athesis.jl
*a work in progress*

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Deltares/Athesis.jl/HEAD?filepath=src/groundwater2d3d.jl)

Experimental portable numerical three-dimensional groundwater solver in Julia.
Currently, the model runs on a structured rectangular grid.
The user can specify whether simulations should be run:
- on CPU or GPU
- using single or double precision

The implentation contains two types of sources:
- point sources
- recharge (rainfall)

Additionally, a number of simple boundary condition types have been implemented. This is still under construction.
In the present imeplementation, an insaturated or partially-saturated zone is not yet considered.

Additionally, some general code snippets can be found in src/sandbox/.

- <b>diffusion2d_CPU_GPU</b> contains a 2d implementation of a simple diffusion solver, that can run both on CPU and GPU.
- <b>multiple_dispatch</b> contains a 1d implementation of a computational kernel that employs multiple dispatch of the kernel (UPDATE) function
based on the contents of "parameters"
Depending on the argument, 8 methods can be invoked
(through multiple dispatch of UPDATE function):
  - explicit advection
  - explicit diffusion
  - (diagonally-)implicit advection
  - (diagonally-)implicit diffusion
  - explicit advection, explicit diffusion
  - (diagonally-)implicit advection, explicit diffusion
  - (diagonally-)implicit advection, explicit diffusion
  - (diagonally-)implicit advection, (diagonally-)implicit diffusion
- <b>loopabstraction</b> prototypes the framework to run loops over the grid on both the CPU and GPU efficiently
