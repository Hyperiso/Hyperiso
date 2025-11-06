#ifndef RNGHelper
#define RNGHelper

#include <vector>
#include <iostream>
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

static constexpr double EPS_SYM = 1e-10;
static constexpr double EPS_DIAG = 1e-8;

Matrix readMatrixFromStdin();

void printVector(const Vector& v);

void printUsage(const char* prog);

#endif