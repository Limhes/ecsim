# ecsim

Fast and general simulation of electrode processes coupled to homogeneous reactions

This is the source code of the algorithm used on [limhes.net/ecsim](http://limhes.net/ecsim), which should also be consulted for a description of the simulator's functionality and inner workings.

## how to use

A full working example of a cyclic voltammetry simulation illustrating all functionality is given in `/general/main.cpp`.

An interface to WebAssembly is implemented in `/webassembly/main.cpp` and can be compiled with the included installer script.

## dependencies

Eigen library version 3

emscripten (for compilation to WebAssembly)
