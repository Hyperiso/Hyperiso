# Installation guide

## System packages

On Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential cmake ninja-build pkg-config \
  libgsl-dev libeigen3-dev \
  python3-dev python3-pip python3-venv
```

Optional documentation tools:

```bash
sudo apt-get install -y doxygen graphviz
```

## C++ core

```bash
cmake -S Hyperiso/Hyperiso/core -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_WITH_CLI=ON

cmake --build build -j
cmake --install build --prefix "$HOME/.local"
```

## Python package

```bash
python -m pip install --upgrade pip build
python -m pip install ./Hyperiso/Hyperiso
```

The build backend installs `pybind11` automatically in an isolated build
environment. The unrelated package named `pybind` must not be installed.

For an editable development installation:

```bash
python -m pip install -e "./Hyperiso/Hyperiso[test,dev]"
```

## Optional backends

Optional backends should be enabled only when their dependencies are installed and the related tests are needed.

```bash
cmake -S Hyperiso/Hyperiso/core -B build-marty \
  -DBUILD_WITH_MARTY=ON
```

```bash
cmake -S Hyperiso/Hyperiso/core -B build-softsusy \
  -DBUILD_WITH_SOFTSUSY=ON
```

## Verify the installation

```bash
hyperiso-ui --version
python - <<'PY'
import pyhyperiso
from pyhyperiso.phyperiso import pyhyperiso as native

assert pyhyperiso.__version__ == native.__version__ == "1.0.0"
print(f"pyhyperiso {pyhyperiso.__version__} import OK")
PY
```
