# Porteau.jl
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://openearth.github.io/Porteau.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://openearth.github.io/Porteau.jl/dev)
[![Build Status](https://travis-ci.com/evetion/Porteau.jl.svg?branch=master)](https://travis-ci.com/openearth/Porteau.jl)


*a work in progress*

Experimental portable water related numerical solver in Julia.

Some first code snippets can be found in src/sandbox/.

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
