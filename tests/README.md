# Unit tests of ecsim using pyecsim and pytest

## Testing of a single redox step

Of the system Ox + n*e <-> Red

### Reversible (Nernstian) charge transfer: Er

Test script: `test_Er.py`

Results: `test_Er_results.txt`

This test makes sure that, for a reversible (Nernstian) redox couple, the simulation gives peak currents that deviate less than 1% from the Randles-Sevcik equation and gives peak potentials that are within one step potential from the theoretical value.

![ip for Er](/tests/formulae/formula_Er_ip.png)

![Ep for Er](/tests/formulae/formula_Er_Ep.png)

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

To run the test, install pytest and run `pytest -v test_Er.py` which gives on my laptop: `3402 passed in 74.92s (0:01:14)`

### Quasi-reversible charge transfer: Eq

Since there are no analytical expressions to test a quasi-reversible redox couple
against, this test is a VISUAL CHECK, and is not run with pytest.

The functions `generate_plot_X()` will generate plots of parameter X at selected
values of alpha and Lambda. An overlay of each plot (in blue) onto the corresponding
figure from Bard & Faulkner is stored as `test_Eq_X.png`, where X is:

```
X = psi     Fig. 6.4.1
X = kappa   Fig. 6.4.2
X = ksi     Fig. 6.4.3 (used DeltaTheta = 0.05 to generate smoother graph)
```

![Psi vs Lambda](/tests/test_Eq_psi.png)

![Kappa vs Lambda](/tests/test_Eq_kappa.png)

![Ksi vs Lambda](/tests/test_Eq_ksi.png)

FYI: I do not have permission from Wiley and Sons to use these figures.

### Fully irreversible charge transfer: Ei

Test script: `test_Ei.py`

Results: `test_Ei_results.txt`

This test makes sure that, for a fully irreversible redox couple, the simulation gives peak currents that deviate less than 1% from the theoretical value and gives peak potentials that are within one step potential from the theoretical value.

![ip for Ei](/tests/formulae/formula_Ei_ip.png)

![Ep for Ei](/tests/formulae/formula_Ei_Ep.png)

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

To run the test, install pytest and run `pytest -v test_Ei.py` which gives on my laptop: `162 passed in 2.87s`

## Charge transfer coupled to homogeneous reactions

### Charge transfer preceding chemical reactions: ErC

Test script: `test_EC.py`

Results: `test_EC_results.txt`

Of the system Ox + n*e <-> Red <-> Prod

This test makes sure that, for a reversible (Nernstian) redox couple followed by a chemical reaction, the simulation gives peak currents that deviate less than 1%/2.5% from the theoretical values and gives peak potentials that are within one step potential from the theoretical value.

Following the nomenclature of Savéant, we have a reaction with k_f and k_b: equilibrium constant K = k_f/k_b and k = k_f+k_b. For the overall system we have the kinetic parameter lambda = R*T/F * k/nu. There are two well-defined regions: KP for ErCi and DE for ErCr with large K and large k. These are distinct cases, and will be treated separately.

The tests pass for any combination of parameters in the ranges:

```
1.0e-2 < nu < 1.0e3                     nu: CV scan rate [V/s]
1.0e-3 < C < 1.0e3                      C: concentration [mol/m^3 or mM]
```

Fixed parameters that were tested in the reversible test case which were assumed to not
influence the current test appreciably were:

```
T == 293.15                 T: temperature [K]
r == 1.0e-3                 r: electrode radius [m]
alpha == 0.5                alpha: symmetry factor of the redox step [-]
D == 1.0e-9                 Dred & Dox & Dprod: diffusion coefficient [m/s^2]
```

Fixed parameter because no theoretical values for peak current or potential could be found
for any other values than:

```
n == 1                      n: #electrons in redox step [-]
```

#### ErCi

![ip for ErCi](/tests/formulae/formula_ErCi_ip.png)

![Ep for ErCi](/tests/formulae/formula_ErCi_Ep.png)

Parameter values for ErCi (EC in KP):

```
max(1.0, k_f_min) < k_f < 1.0e10    k_f: forward rate constant [/s]
    k_f_min was determined from Savéant (Elements of ...) section 2.2
    to force the system into the pure kinetic zone
AND k_b == 0                        k_b: backward rate constant [/s]
```

For all values in the KP zone, the ErCi peak current simulates within 2.5% error.
In the center of the KP zone, the ErCi peak current simulates within 1% error.

#### ErCr

![ip for ErCr](/tests/formulae/formula_ErCr_ip.png)

![Ep for ErCr](/tests/formulae/formula_ErCr_Ep.png)

Parameter values for ErCr (EC in DE, which seamlessly transitions into DO to become indistinguishable from Er):
Selected parameter values from Fig. 2.1 (zone diagram for EC) from Savéant (Elements of ...).

```
(loglambda, logK) in [(4.0, -2.0), (10.0, -2.0), (5.0, 0.0), (10.0, 0.0), (10.0, 2.5)]
```

For all values in the DE zone, the ErCr peak current simulates within 1.5% error.
In the center of the DE zone, the ErCr peak current simulates within 1% error.
