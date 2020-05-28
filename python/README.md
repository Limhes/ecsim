# pyecsim

A Python wrapper around ecsim. Documentation can be found [here](https://limhes.github.io/pyecsim)

## how to use

Import the module using e.g. `import pyecsim as ecs` and start simulating. Have a look at the following examples:

* See the [Jupyter Notebook](tutorial.ipynb) for an interactive example
* See `main.py` for an example script

## installation

To build and install the module, you'll need to install `python3`, the `eigen3` library, `git`, `pybind11` and Python's `setuptools`.

```
git clone https://github.com/Limhes/ecsim
cd ecsim/python/
python3 setup.py install --user
```
