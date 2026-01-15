#ifndef RNGHelper
#define RNGHelper

#include <vector>
#include <iostream>
#include <iomanip>
#include "Matrix.h"

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

static constexpr double EPS_SYM = 1e-10;
static constexpr double EPS_DIAG = 1e-8;

Matrix readMatrixFromStdin();

void printVector(const Vector& v);
void printRealMatrix(const RealMatrix& A);

void printUsage(const char* prog);

#endif