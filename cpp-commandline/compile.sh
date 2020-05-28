#!/bin/bash

g++ ./main.cpp ../source/system.cpp ../source/electrodes.cpp ../source/experiment.cpp ../source/coefs_alpha_beta.cpp ../source/simulation.cpp ../source/simulationcore.cpp -I ../source/ -O3 -std=c++11
