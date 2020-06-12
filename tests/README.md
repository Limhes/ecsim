# Unit tests of ecsim using pyecsim and pytest

## Testing of a single redox step

Of the system Red <-> Ox + n*e

### Reversible (Nernstian) charge transfer

Test script: `test_reversible.py`

Results: `test_reversible_results.txt`

This test makes sure that, for a reversible (Nernstian) redox couple, the simulation gives peak currents that deviate less than 1% from the Randles-Sevcik equation and gives peak potentials that are within one step potential from the theoretical value.

The test passes for any combination of parameters in the ranges:

```
100 < T < 800           T: temperature [K]
1.0e-6 < r < 1.0e-2     r: electrode radius [m]
1.0e-2 < nu < 1.0e3     nu: CV scan rate [V/s]
0.3 < alpha < 0.7       alpha: symmetry factor of the redox step [-]
1.0e-3 < C < 1.0e3      C: concentration [mol/m^3 or mM]
With either:            n: #electrons in redox step [-], Dred & Dox: diffusion coefficient [m/s^2]
    1.0e-7 < Dox,Dred < 1.0e-11 AND 1.0e-2 < Dox/Drex < 1.0e2 AND n == 1
    1.0e-7 < Dox,Dred < 1.0e-11 AND Dox/Drex == 1.0 AND 1 < n < 8
```

The test breaks with unequal diffusion constants and n > 1. The Randles-Sevcik equation is, as far as I am aware of, not necessarily valid in this situation (I cannot find information on this particular case). It should be noted that ANY scenario where n > 1 is regarded as non-existent/non-physical by some electrochemists, and the recommended way to simulate such systems is by using multiple, successive one-electron redox steps, each with their own symmetry factor, heterogeneous rate constant and standard potential.

The test also breaks when very low temperatures (T << 100 K) are combined with n > 2. This is probably an intrinsic limitation of the simulator. Since these use cases seem rare (and not necessary, see previous note), I will not attempt to fix this.

To run the test, install pytest and run `pytest -v test_reversible.py` which gives on my laptop: `3402 passed in 74.92s (0:01:14)`

### Quasi-reversible charge transfer

Since there are no analytical expressions to test a quasi-reversible redox couple
against, this test is a VISUAL CHECK, and is not run with pytest.

The functions generate_plot_X() will generate plots of parameter X at selected
values of alpha and Lambda. An overlay of each plot (in blue) onto the corresponding
figure from Bard & Faulkner is stored as test_quasireversible_X.png, where X is:
    X = psi     Fig. 6.4.1
    X = kappa   Fig. 6.4.2
    X = ksi     Fig. 6.4.3 (used DeltaTheta = 0.05 to generate smoother graph)

![Psi vs Lambda](/tests/test_quasireversible_psi.png)

![Kappa vs Lambda](/tests/test_quasireversible_kappa.png)

![Ksi vs Lambda](/tests/test_quasireversible_ksi.png)

FYI: I do not have permission from Wiley and Sons to use these figures.

### Fully irreversible charge transfer

Test script: `test_irreversible.py`

Results: `test_irreversible_results.txt`

This test makes sure that, for a fully irreversible redox couple, the simulation gives peak currents that deviate less than 1% from the theoretical value and gives peak potentials that are within one step potential from the theoretical value.

The test passes for any combination of parameters in the ranges:

```
1.0e-2 < nu < 1.0e2         nu: CV scan rate [V/s]
0.3 < alpha < 0.7           alpha: symmetry factor of the redox step [-]
1.0e-7 < Dox,Dred < 1.0e-11 AND 0.1 < Dox/Drex < 10.
                            Dred & Dox: diffusion coefficient [m/s^2]
1.0e-10 < k_e < k_e_max     k_e: heterogeneous rate constant [m^2/s]
    k_e_max was determined from eq. 6.4.3 in Bard & Faulkner,
    setting Lambda_max at 10^(-2*alpha) according to the end of section 6.4.
```

Fixed parameters that were tested in the reversible test case which were assumed to not
influence the current test appreciably were:

```
T == 293.15                 T: temperature [K]
r == 1.0e-3                 r: electrode radius [m]
C == 1.0                    C: concentration [mol/m^3 or mM]
```

Fixed parameter because no theoretical values for peak current or potential could be found
for any other values than:

```
n == 1                      n: #electrons in redox step [-]
```

The test breaks for very asymmetric (alpha == 0.7) charge transfer at high scan rates (nu == 100.0) and very low heterogeneous rate constant (k_e == 1.0e-10). Could be either due to theory breaking or simulation breaking at these values.

To run the test, install pytest and run `pytest -v test_irreversible.py` which gives on my laptop: `156 passed in 2.35s`

