#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 11:06:45 2020

@author: limhes
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyecsim as ecs

if __name__ == '__main__':

    sim = ecs.Simulation(True)  # bool verbosity

    A = ecs.Species('A', 1.0, 1.0e-9)
    B = ecs.Species('B', 0.0, 1.0e-9)
    C = ecs.Species('C', 0.0, 1.0e-9)

    rxn1 = ecs.Reaction(B, None, C, None, 10.0, 0.0).enable()
    sim.sys.addReaction(rxn1)

    rdx1 = ecs.Redox(A, B, 1, -0.5, 10.0, 0.5).enable()
    sim.sys.addRedox(rdx1)

    sim.el.disk(1.0e-3)

    sim.exper.setScanPotentials(0.0, [-0.7], 0.0)
    sim.exper.setScanRate(1.0)

    k_fs = np.logspace(-5.0, 9.0, num=10)
    E_pc = []

    cmap_conc = mpl.cm.get_cmap('nipy_spectral')

    plt.figure()
    for index, k_f in enumerate(list(k_fs)):
        cval = 1.0-index/10
        rxn1.setKf(k_f)
        [potential, current] = sim.run()
        plt.plot(potential, [i*1.0e6 for i in current], color=cmap_conc(cval))
        E_pc.append(sim.metrics()[3])  # metrics() returns [i_pa, E_pa, i_pc, E_pc]
    plt.xlabel('Potential [V]')
    plt.ylabel(r'Current [$\mu$A]')
    plt.title(r'CVs of the $E_{r}C_{i}$ scheme at various $k_f$ @ $\nu$ = 1.0 V/s')

    # print variation of E_peak,cathodic with forward rate constant:
    E_pc = np.array(E_pc)
    plt.figure()
    plt.semilogx(k_fs, E_pc, '.')
    plt.xlabel('Forward rate constant [1/s]')
    plt.ylabel('Cathodic peak potential [V]')
    plt.title(r'$E_{peak}$ vs. $k_f$ for the $E_{r}C_{i}$ scheme @ $\nu$ = 1.0 V/s')
