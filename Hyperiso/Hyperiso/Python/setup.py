# setup.py
from skbuild import setup

setup(
    name="phyperiso",
    version="0.1.0",
    description="Hyperiso Python binding (C++ backend via pybind11)",
    author="Théo Reymermier, Niels Fardeau",
    packages=["pyhyperiso", "pyhyperiso.phyperiso"],
    package_dir={"": "."},
    cmake_install_dir="pyhyperiso/phyperiso",
    python_requires=">=3.7",
)