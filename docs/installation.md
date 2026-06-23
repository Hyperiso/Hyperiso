# Installation guide

## System packages

On Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake ninja-build libgsl-dev python3-dev python3-pip
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
  -DBUILD_WITH_APP=ON

cmake --build build -j
cmake --install build --prefix "$HOME/.local"
```

## Python package

```bash
python -m pip install --upgrade pip build
python -m pip install ./Hyperiso/Hyperiso
```

## Optional backends

Optional backends should be enabled only when their dependencies are installed and the related tests are needed.

```bash
cmake -S Hyperiso/Hyperiso/core -B build-marty \
  -DBUILD_WITH_MARTY=ON
```

```bash
cmake -S Hyperiso/Hyperiso/core -B build-2hdmc \
  -DBUILD_WITH_2HDMC=ON
```

```bash
cmake -S Hyperiso/Hyperiso/core -B build-softsusy \
  -DBUILD_WITH_SOFTSUSY=ON
```

## Verify the installation

```bash
hyperiso-ui --help
python -c "import pyhyperiso; print('pyhyperiso import OK')"
```
