#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 12JUN2020

@author: limhes

UNIT TEST FOR THE SIMULATION OF A REVERSIBLE (NERSTIAN) REDOX COUPLE

Notes:
This test makes sure that, for a reversible (Nernstian) redox couple, the simulation
gives peak currents that deviate less than 1% from the Randles-Sevcik equation and gives
peak potentials that are within one step potential from the theoretical value.

The test passes for any combination of parameters in the ranges:
100 < T < 800           T: temperature [K]
1.0e-6 < r < 1.0e-2     r: electrode radius [m]
1.0e-2 < nu < 1.0e3     nu: CV scan rate [V/s]
0.3 < alpha < 0.7       alpha: symmetry factor of the redox step [-]
1.0e-3 < C < 1.0e3      C: concentration [mol/m^3 or mM]
With either:            n: #electrons in redox step [-], Dred & Dox: diffusion coefficient [m/s^2]
    1.0e-7 < Dox,Dred < 1.0e-11 AND 1.0e-2 < Dox/Drex < 1.0e2 AND n == 1
    1.0e-7 < Dox,Dred < 1.0e-11 AND Dox/Drex == 1.0 AND 1 < n < 8

The test breaks with unequal diffusion constants and n > 1. The Randles-Sevcik equation is,
as far as I am aware of, not necessarily valid in this situation (I cannot find information
on this particular case). It should be noted that ANY scenario where n > 1 is regarded as
non-existent/non-physical by some electrochemists, and the recommended way to simulate such
systems is by using multiple, successive one-electron redox steps, each with their own
symmetry factor, heterogeneous rate constant and standard potential.

The test also breaks when very low temperatures (T << 100 K) are combined with n > 2. This is
probably an intrinsic limitation of the simulator. Since these use cases seem rare (and not
necessary, see previous note), I will not attempt to fix this.

To run the test, install pytest and run `pytest -v test_reversible.py` which gives on my laptop:
======================= 3402 passed in 74.92s (0:01:14) ========================

"""

import pytest
import numpy as np
import pyecsim as ecs

F = 96485.3
R = 8.31446

test_params = []
for T in [100.0, 293.15, 800.0]:
    for r in [1.0e-6, 1.0e-4, 1.0e-2]:
        for nu in [0.01, 1.0, 1000.0]:
            for alpha in [0.3, 0.5, 0.7]:
                for C in [1.0e-3, 1.0, 1.0e3]:
                    for D in [1.0e-8, 1.0e-10]:
                        for Ddiff in [0.1, 1.0, 10.0]:
                            test_params.append([T, r, nu, 1, alpha, C, D/Ddiff, D*Ddiff])
                    for D in [1.0e-7, 1.0e-11]:
                        for n in [1, 2, 4, 8]:
                            test_params.append([T, r, nu, n, alpha, C, D, D])

@pytest.mark.parametrize('T,r,nu,n,alpha,C,Dred,Dox', test_params)
def test_reversible(T, r, nu, n, alpha, C, Dred, Dox):
    global R, F
    f = float(n)*F/R/T
    
    sim = ecs.Simulation(False)
    sim.setPotentialSizing(0.2/float(n)) # TODO: this should happen inside the algorithm
    sim.env.setTemperature(T)
    sim.el.disk(r)
    
    ox = ecs.Species('ox', 0.0, Dox)
    red = ecs.Species('red', C, Dred)
    rdx1 = ecs.Redox(ox, red, int(n), 0.0, 10.0, alpha).enable()
    sim.sys.addRedox(rdx1)
    
    sim.exper.setScanPotentials(-0.5, [0.5], -0.5)
    sim.exper.setScanRate(nu)

    [potential, current] = sim.run()
    metrics = sim.metrics()
    
    # Bard & Faulkner chapter 6:
    i_p_theoretical = 0.4463*float(n)*F*(np.pi*r*r)*C*np.sqrt(f*nu*Dred) # Randles-Sevcik equation
    E_half = np.log(np.sqrt(Dred/Dox)) / f
    E_p_theoretical = E_half + 1.109 / f

    i_p_err = np.abs(2*(i_p_theoretical - metrics[0])/(i_p_theoretical + metrics[0]))
    E_p_diff = np.abs(E_p_theoretical - metrics[1]) / np.abs(potential[0]-potential[1])
    
    assert i_p_err < 0.01 # less than 1% error in peak current
    assert E_p_diff < 1.0 # peak potential within step potential
    
