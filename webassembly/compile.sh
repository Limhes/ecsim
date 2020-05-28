#!/bin/bash

em++ ./main.cpp ../source/system.cpp ../source/electrodes.cpp ../source/experiment.cpp ../source/coefs_alpha_beta.cpp ../source/matrixsystem.cpp ../source/simulation.cpp -I ../source/ -I ~/.local/include/emscripten -o ./libwasmecsim.js -O3 -std=c++11 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['ccall', 'cwrap']"
