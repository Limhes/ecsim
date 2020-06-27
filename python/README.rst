PyECsim is a fast and general simulator of voltammetry experiments. Using the latest state-of-the-art algorithms, it simulates

- any number of electrode processes, using Butler-Volmer charge transfer kinetics, coupled to
- any number of homogeneous reactions (first or second order), using
- user-defined voltammetry waveforms.

A typical simulation, even with very large reaction rate constants (10e9 /s), takes around 10 ms to run. Because the simulations are fast, you can use them inside the target function of an optimization routine.

In-depth documentation can be found on `read the docs <https://pyecsim.readthedocs.io/>`_.

# Install
---------

Installing PyECsim is easy:

.. code-block::

    pip install --user pyecsim

# Quick-start
-------------

The below example is a minimal piece of code that simulates a charge transfer step (ox + e <-> red) coupled to a chemical reaction (ox <-> prod).

Import PyECsim into your Python code (and matplotlib to plot the voltammogram):

.. code-block:: python

    import pyecsim as ecs
    import matplotlib.pyplot as plt

Then create a simulator and species, a charge transfer step and a chemical reaction:

.. code-block:: python

    sim = ecs.Simulation(True) # or False to suppress output

    red = ecs.Species('reduced species', 1.0, 1.0e-9) # name, concentration [mol/m3], diffusion coefficient [m2/s]
    ox = ecs.Species('oxidized species', 0.0, 1.0e-9)
    prod = ecs.Species('reaction product', 0.0, 1.0e-9)
    
    rdx1 = ecs.Redox(ox, red, 1, 0.0, 10.0, 0.5).enable() # ox, red, n, E_0 [V], k_e [m/s], alpha
    sim.sys.addRedox( rdx1 )
    
    rxn1 = ecs.Reaction(ox, None, prod, None, 5.0, 0.0).enable() # reactant1, reactant2, product1, product2, k_f, k_b
    sim.sys.addReaction( rxn1 )

Lastly, we set the electrode type and radius, the voltammograms initial/vertex/final potentials and the scan rate:

.. code-block:: python

    sim.el.disk(1.0e-3) # radius [m]
    sim.exper.setScanPotentials(-0.5, [0.5], -0.5) # potentials [V]: initial, [0 or more vertices], final 
    sim.exper.setScanRate(1.0) # scan rate [V/s]

And then we can run the simulation and plot the results:

.. code-block:: python

    [potential, current] = sim.run()
    plt.plot(potential, current)

Which gives:

.. image:: http://limhes.net/upload/images_ecsim/cv_readme.png
