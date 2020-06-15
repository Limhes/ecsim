#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 13JUN2020

@author: limhes

UNIT TEST FOR THE SIMULATION OF A QUASI-REVERSIBLE REDOX COUPLE

Notes:
Since there are no analytical expressions to test a quasi-reversible redox couple
against, this test is a VISUAL CHECK, and is not run with pytest.

The functions generate_plot_X() will generate plots of parameter X at selected
values of alpha and Lambda. An overlay of each plot (in blue) onto the corresponding
figure from Bard & Faulkner is stored as test_Eq_X.png, where X is:
    X = psi     Fig. 6.4.1
    X = kappa   Fig. 6.4.2
    X = ksi     Fig. 6.4.3 (used DeltaTheta = 0.05 to generate smoother graph)

FYI: I do not have permission from Wiley and Sons to use these figures.

"""

import numpy as np
import pyecsim as ecs
import matplotlib.pyplot as plt

F = 96485.3
R = 8.31446
T = 293.15
f = F/R/T

def simulate_Eq(T, r, nu, alpha, k_e, C, Dred, Dox, deltaTheta=0.2):
    sim = ecs.Simulation(True)
    sim.env.setTemperature(T)
    sim.el.disk(r)
    
    ox = ecs.Species('ox', 0.0, Dox)
    red = ecs.Species('red', C, Dred)
    rdx1 = ecs.Redox(ox, red, 1, 0.0, k_e, 1-alpha).enable()
    sim.sys.addRedox(rdx1)
    
    sim.exper.setScanPotentials(-0.2, [], 1.0)
    sim.exper.setScanRate(nu)
    
    sim.setPotentialSizing(deltaTheta)
    [potential, current] = sim.run()
    return [potential, current, sim.metrics()]
    
def generate_plot_psi():
    fig, ax = plt.subplots(3,1,figsize=(10,18))
    for plot_position, alpha in enumerate([0.7, 0.5, 0.3]):
        for Lambda in [0.01, 0.1, 1.0, 10.0]:
            Dox = 1.0e-9
            Dred = 1.0e-9
            nu = 1.0
            r = 1.0e-3
            C = 1.0
            
            k_e = Lambda*np.sqrt(np.power(Dox,1-alpha)*np.power(Dred,alpha)*F/R/T*nu)
        
            [potential, current, metrics] = simulate_quasireversible(T, r, nu, alpha, k_e, C, Dred, Dox)
            
            i_to_psi = F*(np.pi*r*r)*C*np.sqrt(f*nu*Dred)
            psi = [i/i_to_psi for i in current]
            
            ax[plot_position].plot(potential, psi, 'k-')
            
        ax[plot_position].set_xlim(-0.128, 0.46)
        ax[plot_position].set_xticks([-0.128, 0.0, 0.128, 0.257, 0.385])
        ax[plot_position].set_ylim(0.0, 0.5)
        ax[plot_position].grid()
        
def generate_plot_kappa():
    fig, ax = plt.subplots(1,1,figsize=(6,10))
    Lambda_list = list(np.logspace(2.5, -4.5, num=50))
    for alpha in [0.7, 0.6, 0.5, 0.4, 0.3]:
        
        kappa = []
        for Lambda in Lambda_list:
            Dox = 1.0e-9
            Dred = 1.0e-9
            nu = 1.0
            r = 1.0e-3
            C = 1.0
            
            k_e = Lambda*np.sqrt(np.power(Dox,1-alpha)*np.power(Dred,alpha)*f*nu)
        
            [potential, current, metrics] = simulate_quasireversible(T, r, nu, alpha, k_e, C, Dred, Dox)
            
            i_to_irev = 0.4463*F*(np.pi*r*r)*C*np.sqrt(f*nu*Dred)
            kappa.append(metrics[0]/i_to_irev)
            
        ax.semilogx(Lambda_list, kappa, 'k-')
        ax.set_xlim(10**2.5, 10**-4.5)
        ax.set_ylim(0.0, 1.15)
                
def generate_plot_ksi():
    fig, ax = plt.subplots(1,1,figsize=(6,8))
    Lambda_list = list(np.logspace(2.5, -4.5, num=50))
    for alpha in [0.7, 0.6, 0.5, 0.4, 0.3]:
        
        ksi = []
        for Lambda in Lambda_list:
            Dox = 1.0e-9
            Dred = 1.0e-9
            nu = 1.0
            r = 1.0e-3
            C = 1.0
            
            k_e = Lambda*np.sqrt(np.power(Dox,1-alpha)*np.power(Dred,alpha)*F/R/T*nu)
        
            [potential, current, metrics] = simulate_quasireversible(T, r, nu, alpha, k_e, C, Dred, Dox, 0.05)
            E_half = np.log(np.sqrt(Dred/Dox))/f
            ksi.append( (metrics[1]-E_half)*f )
            
        ax.semilogx(Lambda_list, ksi, 'k-')
        ax.set_xlim(10**2.5, 10**-4.5)
        ax.set_ylim(0.0, 20.0)

if __name__ == '__main__':
    #generate_plot_psi()
    #generate_plot_kappa()
    #generate_plot_ksi()
    pass
