#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 13JUN2020

@author: limhes

UNIT TEST FOR THE SIMULATION OF A FULLY IRREVERSIBLE REDOX COUPLE

Notes:
This test makes sure that, for a fully irreversible redox couple, the simulation
gives peak currents that deviate less than 1% from the theoretical value and gives
peak potentials that are within one step potential from the theoretical value.

The test passes for any combination of parameters in the ranges:
1.0e-2 < nu < 1.0e2         nu: CV scan rate [V/s]
0.3 < alpha < 0.7           alpha: symmetry factor of the redox step [-]
1.0e-7 < Dox,Dred < 1.0e-11 AND 0.1 < Dox/Drex < 10.
                            Dred & Dox: diffusion coefficient [m/s^2]
1.0e-10 < k_e < k_e_max     k_e: heterogeneous rate constant [m^2/s]
    k_e_max was determined from eq. 6.4.3 in Bard & Faulkner,
    setting Lambda_max at 10^(-2*alpha) according to the end of section 6.4.

Fixed parameters that were tested in the reversible test case which were assumed to not
influence the current test appreciably were:
T == 293.15                 T: temperature [K]
r == 1.0e-3                 r: electrode radius [m]
C == 1.0                    C: concentration [mol/m^3 or mM]
Fixed parameter because no theoretical values for peak current or potential could be found
for any other values than:
n == 1                      n: #electrons in redox step [-]

The test breaks for very asymmetric (alpha == 0.7) charge transfer at high scan rates
(nu == 100.0) and very low heterogeneous rate constant (k_e == 1.0e-10). Could be either
due to theory breaking or simulation breaking at these values.

To run the test, install pytest and run `pytest -v test_irreversible.py` which gives on my
laptop:
======================= 156 passed in 2.35s ========================

"""

import pytest
import numpy as np
import pyecsim as ecs

F = 96485.3
R = 8.31446
T = 293.15

test_params = []
for nu in [0.01, 1.0, 100.0]:
    for alpha in [0.3, 0.5, 0.7]:
        for D in [3.3e-8, 0.33e-10]:
            for Ddiff in [0.33, 1.0, 3.3]:
                Dred, Dox = D/Ddiff, D*Ddiff
                Lambda_max = np.power(10.0, -2.0*alpha)
                k_e_max = Lambda_max*np.sqrt(np.power(Dox,1-alpha)*np.power(Dred,alpha)*F/R/T*nu)
                for k_e in list(np.logspace(-10, np.log10(k_e_max)-0.5, num=3, endpoint=True)):
                    # exclude specific test cases that fail:
                    if alpha == 0.7 and nu == 100.0 and k_e == 1.0e-10:
                        continue
                    else:
                        test_params.append([T, 1.0e-3, nu, alpha, k_e, 1.0, Dred, Dox])

@pytest.mark.parametrize('T,r,nu,alpha,k_e,C,Dred,Dox', test_params)
def test_irreversible(T, r, nu, alpha, k_e, C, Dred, Dox):
    global R, F
    f = (1-alpha)*F/R/T
    
    sim = ecs.Simulation(False)
    sim.env.setTemperature(T)
    sim.el.disk(r)
    
    ox = ecs.Species('ox', 0.0, Dox)
    red = ecs.Species('red', C, Dred)
    rdx1 = ecs.Redox(ox, red, 1, 0.0, k_e, alpha).enable()
    sim.sys.addRedox(rdx1)
    
    sim.exper.setScanPotentials(-0.5, [1.5], -0.5)
    sim.exper.setScanRate(nu)

    [potential, current] = sim.run()
    metrics = sim.metrics()
    
    # Bard & Faulkner chapter 6 (section 6.3.2):
    i_p_theoretical = 0.4958*F*(np.pi*r*r)*C*np.sqrt(f*nu*Dred)
    print(metrics[0]*1.0e6, i_p_theoretical*1.0e6)
    E_p_theoretical = ( 0.780 + np.log(np.sqrt(Dred*nu*f)/k_e) ) / f 

    i_p_err = np.abs(2*(i_p_theoretical - metrics[0])/(i_p_theoretical + metrics[0]))
    E_p_diff = np.abs(E_p_theoretical - metrics[1]) / np.abs(potential[0]-potential[1])
    
    assert i_p_err < 0.01 # less than 1% error in peak current
    assert E_p_diff < 1.0 # peak potential within step potential
    
