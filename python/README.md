# pyecsim

A Python wrapper around ecsim. Documentation can be found [here](https://limhes.github.io/pyecsim)

## how to use

Import the module using e.g. `import pyecsim as ecs` and start simulating. Have a look at the following examples:

* See the [Jupyter Notebook](tutorial.ipynb) for an interactive example
* See `main.py` for an example script

## installation

### option 1 - download and install

This is the easiest option and recommended for Windows users. Download the `.msi` file (32-bit or 64-bit) [here](dist) and run the installer.

### option 2 - build on your PC

This is a more advanced option and recommended for Linux users. To build and install the module, you'll need `python3`, the `eigen3` library, `git`, `pybind11` and Python's `setuptools`.

```
git clone https://github.com/Limhes/ecsim
cd ecsim/python/
python3 setup.py install --user
```
