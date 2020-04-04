# ecsim

Fast and general simulation of electrode processes coupled to homogeneous reactions

This is the source code of the algorithm used on [limhes.net/ecsim](http://limhes.net/ecsim), which should also be consulted for a description of the simulator's functionality and inner workings.

## how to use

A full working example of a cyclic voltammetry simulation illustrating all functionality is given in `/general/main.cpp`.

An interface to WebAssembly is implemented in `/webassembly/main.cpp` and can be compiled with the included installer script.

A Qt5 GUI application can be found in `/qt5/` which can be built using `qmake QESP` and then `make`. This version uses a slightly different source code (e.g. some classes inherit QObject to make use of signal/slot) and this application should thus be viewed as an entity separate from the main code base.

## dependencies

Dependencies (with example Arch Linux installation commands) are given below:

1. Eigen version 3 (`sudo pacman -S eigen`)
2. For the Qt5 version: Qt5 (`sudo pacman -S qt5-base qt5-tools`)
3. For compilation to WebAssembly: emscripten (`sudo pacman -S emscripten`)
