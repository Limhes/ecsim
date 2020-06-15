#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 12JUN2020

@author: limhes

UNIT TEST FOR THE SIMULATION OF A REVERSIBLE (NERSTIAN) REDOX COUPLE FOLLOWED BY A
HOMOGENEOUS CHEMICAL REACTION

Notes:
This test makes sure that, for a reversible (Nernstian) redox couple followed by a chemical
reaction, the simulation gives peak currents that deviate less than 1%/2.5% from the theoretical values
and gives peak potentials that are within one step potential from the theoretical value.

The test passes for any combination of parameters in the ranges:
1.0e-2 < nu < 1.0e3                     nu: CV scan rate [V/s]
1.0e-3 < C < 1.0e3                      C: concentration [mol/m^3 or mM]

Fixed parameters that were tested in the reversible test case which were assumed to not
influence the current test appreciably were:
T == 293.15                 T: temperature [K]
r == 1.0e-3                 r: electrode radius [m]
alpha == 0.5                alpha: symmetry factor of the redox step [-]
D == 1.0e-9                 Dred & Dox & Dprod: diffusion coefficient [m/s^2]
Fixed parameter because no theoretical values for peak current or potential could be found
for any other values than:
n == 1                      n: #electrons in redox step [-]

Parameter values for ErCi (EC in KP):
max(1.0, k_f_min) < k_f < 1.0e10    k_f: forward rate constant [/s]
    k_f_min was determined from Savéant (Elements of ...) section 2.2
    to force the system into the pure kinetic zone
AND k_b == 0                        k_b: backward rate constant [/s]
In the center of the KP zone, the ErCi peak current simulates within 1% error.

Parameter values for ErCr (EC in DE, which seamlessly transitions into DO to become indistinguishable
from Er):
Selected parameter values from Fig. 2.1 (zone diagram for EC) from Savéant (Elements of ...).
(loglambda, logK) in [(4.0, -2.0), (10.0, -2.0), (5.0, 0.0), (10.0, 0.0), (10.0, 2.5)]
In the center of the DE zone, the ErCr peak current simulates within 1% error.

To run the test, install pytest and run `pytest -v test_EC.py` which gives on my laptop:
======================= 105 passed in 2.45s ========================

"""

import pytest
import numpy as np
import pyecsim as ecs

F = 96485.3
R = 8.31446
T = 293.15
f = F/R/T # ~ 40
r = 1.0e-3
D = 1.0e-9
n = 1
alpha = 0.5

test_params_ErCi = []
for nu in [0.01, 1.0, 1000.0]:
    lambda_min = 100.0 # from DO to KP at high K
    k_f_min = lambda_min*nu/f
    for C in [1.0e-3, 1.0, 1.0e3]:
        for k_f in list( np.logspace(np.log(k_f_min)+2.0, 10.0, num=7, endpoint=True) ):
            if k_f < 1.0:
                continue
            else:
                test_params_ErCi.append([nu,C,k_f])

test_params_ErCr = []
for nu in [0.01, 1.0, 1000.0]:
    for C in [1.0e-3, 1.0, 1.0e3]:
        for (loglambda, logK) in [(4.0, -2.0), (10.0, -2.0), (5.0, 0.0), (10.0, 0.0), (10.0, 2.5)]:
            k = (10.0**loglambda)*nu/f
            k_b = k/(1.0+10.0**logK)
            k_f = k-k_b
            print(k_f, k_b)
            test_params_ErCr.append([nu,C,k_f,k_b])



@pytest.mark.parametrize('nu,C,k_f', test_params_ErCi)
def test_ErCi(nu, C, k_f):
    
    sim = ecs.Simulation(True)
    sim.env.setTemperature(T)
    sim.el.disk(r)
    
    red = ecs.Species('red', C, D)
    ox = ecs.Species('ox', 0.0, D)
    prod = ecs.Species('prod', 0.0, D)
    rdx1 = ecs.Redox(ox, red, n, 0.0, 10.0, alpha).enable()
    sim.sys.addRedox(rdx1)
    rxn1 = ecs.Reaction(ox, None, prod, None, k_f, 0.0).enable()
    sim.sys.addReaction(rxn1)
    
    sim.exper.setScanPotentials(-0.5, [0.5], -0.5)
    sim.exper.setScanRate(nu)

    [potential, current] = sim.run()
    metrics = sim.metrics()
    
    # Savéant (Elements of ...) page 83:
    i_p_theoretical = 0.496*F*(np.pi*r*r)*C*np.sqrt(f*nu*D)
    print(metrics[0], i_p_theoretical)
    E_p_theoretical = (0.78 - np.log(k_f/nu/f)/2) / f
    print(metrics[1], E_p_theoretical)

    i_p_err = np.abs(2*(i_p_theoretical - metrics[0])/(i_p_theoretical + metrics[0]))
    E_p_diff = np.abs(E_p_theoretical - metrics[1]) / np.abs(potential[0]-potential[1])
    
    assert i_p_err < 0.025 # less than 2.5% error in peak current
    assert E_p_diff < 1.0 # peak potential within step potential

@pytest.mark.parametrize('nu,C,k_f,k_b', test_params_ErCr)
def test_ErCr(nu, C, k_f, k_b):
    
    sim = ecs.Simulation(True)
    sim.env.setTemperature(T)
    sim.el.disk(r)
    
    red = ecs.Species('red', C, D)
    ox = ecs.Species('ox', 0.0, D)
    prod = ecs.Species('prod', 0.0, D)
    rdx1 = ecs.Redox(ox, red, n, 0.0, 10.0, alpha).enable()
    sim.sys.addRedox(rdx1)
    rxn1 = ecs.Reaction(ox, None, prod, None, k_f, k_b).enable()
    sim.sys.addReaction(rxn1)
    
    sim.exper.setScanPotentials(-1.0, [0.5], -1.0)
    sim.exper.setScanRate(nu)
    #sim.setPotentialSizing(0.1)

    [potential, current] = sim.run()
    metrics = sim.metrics()
    
    # Savéant (Elements of ...) page 85:
    i_p_theoretical = 0.4463*F*(np.pi*r*r)*C*np.sqrt(f*nu*D) # Randles-Sevcik equation
    print(metrics[0], i_p_theoretical)
    E_p_theoretical = (1.109 - np.log(1+k_f/k_b)) / f
    print(metrics[1], E_p_theoretical)

    i_p_err = np.abs(2*(i_p_theoretical - metrics[0])/(i_p_theoretical + metrics[0]))
    E_p_diff = np.abs(E_p_theoretical - metrics[1]) / np.abs(potential[0]-potential[1])
    
    assert i_p_err < 0.015 # less than 2.5% error in peak current
    assert E_p_diff < 1.0 # peak potential within 2 step potential
    
