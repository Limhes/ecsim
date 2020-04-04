# ecsim

Fast and general simulation of electrode processes coupled to homogeneous reactions

This is the source code of the algorithm used on http://limhes.net/ecsim

Functionality:
* asdasd
* asdasd

## how to use

A full working example where all functionality is used, is given in `/general/main.cpp`

A layer to WebAssembly is given in `/webassembly/main.cpp` and can be compiled with the included installer script.

em++ main.cpp system.cpp electrodes.cpp experiment.cpp coefs_alpha_beta.cpp simulation.cpp -o libwasmecsim.js -O3 -std=c++11 -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]'

## dependencies

Eigen

emscripten (for compilation to WebAssembly)

