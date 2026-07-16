# Getting started {#getting_started}

This page covers a clean source installation and the shortest route to the C++,
Python, CLI and graphical interfaces.

## Requirements

A standard source build requires:

- Linux with GCC or Clang supporting C++20;
- CMake 3.20 or newer and Ninja or Make;
- GNU Scientific Library (GSL);
- Eigen 3;
- Python 3.10--3.12 and `pybind11` for the Python package.

MARTY and SOFTSUSY are optional. The bundled 2HDMC and MinuitCpp components are
built only where required by the selected configuration.

## Build the C++ core and CLI

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential cmake ninja-build pkg-config \
  libgsl-dev libeigen3-dev

cmake -S Hyperiso/Hyperiso/core -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_WITH_CLI=ON

cmake --build build --parallel
cmake --install build --prefix "$HOME/.local"
```

Add the installation prefix when needed:

```bash
export PATH="$HOME/.local/bin:$PATH"
export CMAKE_PREFIX_PATH="$HOME/.local:$CMAKE_PREFIX_PATH"
```

Verify the CLI:

```bash
hyperiso-ui --version
hyperiso-ui --help
```

## Install the Python package from source

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ./Hyperiso/Hyperiso
```

For development and tests:

```bash
python -m pip install -e "./Hyperiso/Hyperiso[test,dev]"
python -m pytest Hyperiso/Hyperiso/pyhyperiso/test
```

The unrelated package named `pybind` is not required; HyperIso uses
`pybind11`.

## Verify the Python installation

```bash
python - <<'PY'
import pyhyperiso
from pyhyperiso.phyperiso import pyhyperiso as native

assert pyhyperiso.__version__ == native.__version__
print("pyhyperiso", pyhyperiso.__version__, "import OK")
PY
```

## Run the command-line interface

Wilson coefficients:

```bash
hyperiso-ui wilson summary \
  --model SM \
  --lha Assets/lha/si_input.flha \
  --groups BCoefficients \
  --coeffs C7,C9,C10 \
  --qmatch 81 \
  --q 4.8 \
  --order NNLO
```

Observables:

```bash
hyperiso-ui observable summary \
  --model SM \
  --lha Assets/lha/si_input.flha \
  --observables BR_Bs__mu_mu,BR_B__Xs_gamma \
  --order NNLO
```

## Start the Dash interface

```bash
cd GHyperiso/HyperisoDashGUI
python -m pip install -r requirements.txt
python app.py
```

Open `http://127.0.0.1:8050` in a browser.

## Run tests

```bash
cmake -S Hyperiso/Hyperiso/core -B build-test \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_TESTS=ON
cmake --build build-test --parallel
ctest --test-dir build-test --output-on-failure
```

Heavy numerical comparison tests are opt-in through
`-DHYPERISO_BUILD_COMPARISON_TESTS=ON`.

## Next steps

- Read @ref inputs before changing numerical inputs.
- Read @ref interfaces to choose the most appropriate frontend.
- Read @ref statistics before running fits or confidence contours.
- Use the installed examples under `Hyperiso/Hyperiso/examples_cpp` and
  `Hyperiso/Hyperiso/examples_python` as executable reference workflows.
