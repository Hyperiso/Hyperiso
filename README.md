# Hyperiso

## Introduction
HyperIso is a modern redesign of the SuperIso software, used in particle physics for comparing theory and experiment. This new version aims to enhance the user interface, optimize calculations and data generation, as well as enable the use of new features and integration with other software.

## Features
- Calculation of Wilson coefficients
- Retrieval of SLHA/FLHA data
- Marty usages

## Installation

### Prerequisites
- Operating System: Windows, Linux (not MacOS)
- gcc

### Installation Steps
1. Clone the repository: `git clone https://github.com/Hyperiso/Hyperiso.git`
2. Run cmake: `cmake .`
3. Compile: `make`

## Usage
Execute the command `./hyperiso` and enjoy.

### Examples
- Example 1: Reading SLHA file

## Testing
Tests can be performed by running the command `make test` in the test folder. Two targets, `testWilson` and `testDatabase`, are available.

Feel free to adjust the formatting or provide additional details as needed!
