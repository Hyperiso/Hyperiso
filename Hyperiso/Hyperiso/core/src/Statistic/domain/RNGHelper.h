#ifndef RNG_HELPER_H
#define RNG_HELPER_H

#include <vector>
#include <iostream>
#include <iomanip>

#include "Matrix.h"

/**
 * @file RNGHelper.h
 * @brief Small helper utilities for matrix-based random generation tools.
 *
 * This header provides lightweight I/O and display helpers commonly used by
 * random-number-generation or covariance/correlation-matrix utilities.
 *
 * It includes:
 * - a simple dense matrix alias for raw stdin parsing,
 * - numerical tolerances used when validating input matrices,
 * - helpers to read a square matrix from standard input,
 * - pretty-printers for vectors and @ref RealMatrix objects,
 * - a usage printer for command-line tools.
 *
 * Expected stdin matrix format:
 * @code
 * n
 * a11 a12 ... a1n
 * ...
 * an1 an2 ... ann
 * @endcode
 *
 * @see RealMatrix
 * @see Vector
 */

/// Simple nested-vector matrix type used for raw input parsing.
using Matrix = std::vector<std::vector<double>>;

/**
 * @brief Tolerance used when checking matrix symmetry.
 *
 * This is typically used to decide whether `A(i,j)` and `A(j,i)` should be
 * considered equal up to numerical precision.
 */
static constexpr double EPS_SYM = 1e-10;

/**
 * @brief Tolerance used when checking diagonal normalization.
 *
 * This is typically used for correlation matrices, where diagonal entries are
 * expected to be close to `1.0`.
 */
static constexpr double EPS_DIAG = 1e-8;

/**
 * @brief Reads a square matrix from standard input.
 *
 * The expected format is:
 * @code
 * n
 * r11 r12 ... r1n
 * ...
 * rn1 rn2 ... rnn
 * @endcode
 *
 * The function first reads the matrix size `n`, then reads `n*n` values in
 * row-major order into a nested `std::vector<std::vector<double>>`.
 *
 * @return Parsed square matrix.
 *
 * @throws std::runtime_error if:
 * - the size `n` cannot be read,
 * - `n <= 0`,
 * - one of the matrix coefficients cannot be read.
 */
Matrix readMatrixFromStdin();


/**
 * @brief Prints a vector to standard output on a single line.
 *
 * Output format:
 * - fixed notation,
 * - precision set to 15 digits,
 * - values separated by spaces.
 *
 * @param v Vector to print.
 */
void printVector(const Vector& v);

/**
 * @brief Prints a @ref RealMatrix to standard output.
 *
 * Output format:
 * - scientific notation,
 * - precision set to 5 digits,
 * - one matrix row per line,
 * - each row enclosed in square brackets.
 *
 * Example:
 * @code
 * [1.00000e+00, 2.00000e+00]
 * [3.00000e+00, 4.00000e+00]
 * @endcode
 *
 * @param A Matrix to print.
 */
void printRealMatrix(const RealMatrix& A);

/**
 * @brief Prints command-line usage information for a matrix-based RNG tool.
 *
 * The message describes:
 * - expected command-line arguments,
 * - supported distributions,
 * - expected stdin matrix format,
 * - a usage example.
 *
 * @param prog Program name, typically `argv[0]`.
 */
void printUsage(const char* prog);

#endif  // RNG_HELPER_H