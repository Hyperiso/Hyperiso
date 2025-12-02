# Hyperiso

## Introduction
Hyperiso is a modern redesign of the SuperIso software, used in particle physics for flavour observables calculations. Its greatest benefit is the uses of the Marty software, allowing calculations in 
any BSM models (with BSM masses above 100 GeV). Hyperiso is capable of performing wilson coefficients calculations at matching scale (at Leading Order using Marty, or NNLO for SM, SUSY and THDM), running of these coefficients up TO NNLO, calculation of a lot of flavour observables (mainly from B decays).

## Features
- Calculation of Wilson coefficients
- Retrieval of SLHA/FLHA data
- Marty usages

## Installation

### Prerequisites
- Operating System: Windows, Linux (not MacOS, go buy a real computer)
- gcc
- gsl

### Installation Steps (C++)
1. Clone the repository: `git clone https://github.com/Hyperiso/Hyperiso.git`
2. Go to the repository: `cd hyperiso`
3. Run cmake: `mkdir build && cd build && cmake ../Hyperiso/Hyperiso/core` (with options like -DBUILD_WITH_MARTY=ON, -DBUILD_WITH_2HDMC=ON or -DBUILD_WITH_SOFTSUSY=ON)
3. Compile: `cmake --build .`
4. Install using `cmake --install build --prefix "$HOME/.local"` (the prefix option is used to avoid permission error)

### Installation Steps (Python)
1. Clone the repository: `git clone https://github.com/Hyperiso/Hyperiso.git`
2. Go to the repository: `cd Hyperiso`
3. Run pip install: `pip install Hyperiso/Hyperiso/Python`

## Usage
Depending on the interface, create your main file, link it to the Hyperiso library and do what you want !

### Examples
- Example 1: Reading SLHA file

## Testing
Tests can be performed by running the command `ctest` in the build folder. Multiple targets, like `testWilson` and `testDatabase`, are available.

Feel free to adjust the formatting or provide additional details as needed!
