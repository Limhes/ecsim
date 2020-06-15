#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 12JUN2020

@author: limhes

UNIT TEST FOR THE SIMULATION OF A HOMOGENEOUS CHEMICAL REACTION FOLLOWED BY A
REVERSIBLE (NERSTIAN) REDOX COUPLE 

Notes:
This test makes sure that, for a chemical reaction followed by a reversible (Nernstian)
redox couple, the simulation gives peak currents that deviate less than 1%/1.5% from the
theoretical values and gives peak potentials that are within one step potential from the
theoretical value.

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

Tests pass only for values that are well outside the zone boundaries in Fig. 2.8 from
Savéant (Elements of ...), implying that the zone transitions are wide, as opposed to the
situation for the EC reaction. Visual inspection of simulated CE voltammograms in the
transitional area confirmed this observation.

Selected parameter values for which the tests hold:
DE: (loglambda, logK) in [(7.0, 1.0), (10.0, -1.5)]
KP: (loglambda, logK) in [(3.0, -3.0), (3.0, -6.0), (4.5, -9.0)]

To run the test, install pytest and run `pytest -v test_CE.py` which gives on my laptop:
======================= 45 passed in 1.12s ========================

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

test_params_CE_KP = []
for nu in [0.01, 1.0, 1000.0]:
    for C in [1.0e-3, 1.0, 1.0e3]:
        for (loglambda, logK) in [(3.0, -3.0), (3.0, -6.0), (4.5, -9.0)]:
            k = (10.0**loglambda)*nu/f
            k_b = k/(1.0+10.0**logK)
            k_f = k-k_b
            print(k_f, k_b)
            test_params_CE_KP.append([nu,C,k_f,k_b])

test_params_CE_DE = []
for nu in [0.01, 1.0, 1000.0]:
    for C in [1.0e-3, 1.0, 1.0e3]:
        for (loglambda, logK) in [(7.0, 1.0), (10.0, -1.5)]:
            k = (10.0**loglambda)*nu/f
            k_b = k/(1.0+10.0**logK)
            k_f = k-k_b
            print(k_f, k_b)
            test_params_CE_DE.append([nu,C,k_f,k_b])
            

@pytest.mark.parametrize('nu,C,k_f,k_b', test_params_CE_KP)
def test_CE_KP(nu, C, k_f, k_b):
    
    sim = ecs.Simulation(True)
    sim.env.setTemperature(T)
    sim.el.disk(r)
    
    red = ecs.Species('red', 0.0, D)
    ox = ecs.Species('ox', 0.0, D)
    prod = ecs.Species('prod', C, D)
    rdx1 = ecs.Redox(ox, red, n, 0.0, 10.0, alpha).enable()
    sim.sys.addRedox(rdx1)
    rxn1 = ecs.Reaction(prod, None, red, None, k_f, k_b).enable()
    sim.sys.addReaction(rxn1)
    
    sim.exper.setScanPotentials(-0.5, [0.5], -0.5)
    sim.exper.setScanRate(nu)

    [potential, current] = sim.run()
    metrics = sim.metrics()
    
    # Savéant (Elements of ...) page 93:
    i_p_theoretical = F*(np.pi*r*r)*C*np.sqrt(D*k_b)*(k_f/k_b)
    print(metrics[0], i_p_theoretical)

    i_p_err = np.abs(2*(i_p_theoretical - metrics[0])/(i_p_theoretical + metrics[0]))
    
    assert i_p_err < 0.01 # less than 1% error in plateau current

@pytest.mark.parametrize('nu,C,k_f,k_b', test_params_CE_DE)
def test_CE_DE(nu, C, k_f, k_b):
    
    sim = ecs.Simulation(True)
    sim.env.setTemperature(T)
    sim.el.disk(r)
    
    red = ecs.Species('red', 0.0, D)
    ox = ecs.Species('ox', 0.0, D)
    prod = ecs.Species('prod', C, D)
    rdx1 = ecs.Redox(ox, red, n, 0.0, 10.0, alpha).enable()
    sim.sys.addRedox(rdx1)
    rxn1 = ecs.Reaction(prod, None, red, None, k_f, k_b).enable()
    sim.sys.addReaction(rxn1)
    
    sim.exper.setScanPotentials(-1.0, [1.0], -1.0)
    sim.exper.setScanRate(nu)
    #sim.setPotentialSizing(0.1)

    [potential, current] = sim.run()
    metrics = sim.metrics()
    
    # Savéant (Elements of ...) page 94:
    i_p_theoretical = 0.4463*F*(np.pi*r*r)*C*np.sqrt(f*nu*D) # Randles-Sevcik equation
    print(metrics[0], i_p_theoretical)
    K = k_f/k_b
    E_p_theoretical = (1.109 - np.log(K/(1+K))) / f
    print(metrics[1], E_p_theoretical)

    i_p_err = np.abs(2*(i_p_theoretical - metrics[0])/(i_p_theoretical + metrics[0]))
    E_p_diff = np.abs(E_p_theoretical - metrics[1]) / np.abs(potential[0]-potential[1])
    
    assert i_p_err < 0.015 # less than 1.5% error in peak current
    assert E_p_diff < 1.0 # peak potential within 2 step potential
    
