# ecsim

Fast and general simulation of electrode processes coupled to homogeneous reactions

This is the C++ source code (`/source`) of the algorithm used on [limhes.net/ecsim](http://limhes.net/ecsim), which should also be consulted for a description of the simulator's functionality and inner workings.

## how to use

* **C++ / commandline**: A full working example of a cyclic voltammetry simulation illustrating all functionality is given in `/cpp-commandline/main.cpp`.
* **C++ / Qt GUI**: A Qt5 GUI application can be found in `/cpp-qt5gui/` which can be built using `qmake QESP` and then `make`. This version uses a slightly different source code (e.g. some classes inherit QObject to make use of signal/slot) and this application should thus be viewed as an entity separate from the main code base in `/source`.
* **JavaScript**: An interface to WebAssembly is implemented in `/webassembly/main.cpp` and can be compiled with the included installer script.
* **Python**: A Python module (using pybind11) can be found in `/python/`.

## dependencies

Dependencies (with example Arch Linux installation commands) are given below:

1. Eigen version 3 (`sudo pacman -S eigen`)
2. Optional (for the Qt5 version): Qt5 (`sudo pacman -S qt5-base qt5-tools`)
3. Optional (for compilation to WebAssembly): emscripten (`sudo pacman -S emscripten`)

## CI/CD Workflows

The code in this repository is built automatically using GitHub Actions workflows.
On every pull-request or commit(s) to the main branch a [build workflow](.github/workflows/build.yaml) is triggered to run the Python tests.

On every git tag created (and pushed to the remote repository) a [release workflow](.github/workflows/release.yaml) is triggered to create a new release.
```bash
git tag -m "my-version"
git push --tags 
```
The release workflow will first build both Linux and Windows artifacts, and then will collect those artifacts in a [GitHub release on this repository](https://github.com/Limhes/ecsim/releases).

The building process for Linux leverages a docker image that contains multiple python versions.
The build itself is described in the [Dockerfile](./Dockerfile) committed to the repository.

The building process for Windows leverages GitHub Actions matrix strategy to create different build jobs for each version of Python that the artefacts are built against.
