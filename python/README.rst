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

Import PyECsim into your Python code (and matplotlib to plot the voltammogram):

.. code-block:: python

    import pyecsim as ecs
    import matplotlib.pyplot as plt

Then create a simulator and species, a charge transfer step and a chemical reaction:

.. code-block:: python

    sim = ecs.Simulation(True)

    red = ecs.Species('reduced species', 1.0, 1.0e-9)
    ox = ecs.Species('oxidized species', 0.0, 1.0e-9)
    prod = ecs.Species('reaction product', 0.0, 1.0e-9)

    sim.sys.addRedox( rdx1 := ecs.Redox(ox, red, 1, 0.0, 10.0, 0.5).enable() ) # Python >= 3.8
    sim.sys.addReaction( rxn1 := ecs.Reaction(ox, None, prod, None, 5.0, 0.0).enable() ) # Python >= 3.8

Lastly, we set the electrode type and radius, the voltammograms initial/vertex/final potentials and the scan rate:

.. code-block:: python

    sim.el.disk(1.0e-3)
    sim.exper.setScanPotentials(-0.5, [0.5], -0.5)
    sim.exper.setScanRate(1.0)

And then we can run the simulation and plot the results:

.. code-block:: python

    [potential, current] = sim.run()
    plt.plot(potential, current)

Which gives:

.. image:: http://limhes.net/upload/images_ecsim/cv_readme.png
