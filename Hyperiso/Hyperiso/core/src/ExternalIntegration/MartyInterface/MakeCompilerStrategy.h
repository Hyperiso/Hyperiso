#ifndef MAKE_COMPILER_STRATEGY_H
#define MAKE_COMPILER_STRATEGY_H

#include "CompilerStrategy.h"

#include <filesystem>
#include <utility>

/**
 * @file MakeCompilerStrategy.h
 * @brief Declares a compilation strategy based on `make`.
 *
 * This header defines ::MakeCompilerStrategy, a ::CompilerStrategy
 * that delegates the build process to `make` in a given source directory,
 * and then runs the resulting binary with an optional matching scale.
 */

 /**
 * @class MakeCompilerStrategy
 * @ingroup CodeGenerationModule
 * @brief Compilation strategy using a make-based build system.
 *
 * MakeCompilerStrategy assumes that:
 *  - @p sourceFile passed to ::compile() is a directory containing a Makefile,
 *  - running `make` in that directory builds the desired binary,
 *  - ::compile_run() can then execute the binary with a `-Q` option specifying
 *    a matching scale.
 */
class MakeCompilerStrategy : public CompilerStrategy {
public:
    /**
     * @brief Constructs a MakeCompilerStrategy for a given (model, wilson) pair.
     *
     * @param model   Model name (metadata only).
     * @param wilson  Wilson basis / label (metadata only).
     */
    MakeCompilerStrategy(std::string model, std::string wilson) : CompilerStrategy(model, wilson) {}

    /**
     * @brief Compiles if necessary and then runs the resulting binary.
     *
     * If ::check_if_compile() returns false, ::compile() is called first.
     * The resulting binary is then run with:
     *  `outputBinary -Q <Q_match>`
     *
     * @param sourceFile    Directory containing the Makefile.
     * @param outputBinary  Path/name of the binary to run.
     */
    void compile_run(const std::string& sourceFile, const std::string& outputBinary) override;

    /**
     * @brief Invokes `make` in the given source directory.
     *
     * @param sourceFile    Directory in which to run `make`.
     * @param outputBinary  Unused here except for consistency with interface;
     *                      the Makefile determines the real output name.
     */
    void compile(const std::string& sourceFile, const std::string& outputBinary) override;

    /**
     * @brief Sets the matching scale @f$Q_{\text{match}}@f$ passed to the binary.
     *
     * The scale is used by ::compile_run() when executing:
     *  `outputBinary -Q <Q_match>`.
     *
     * @param Q_match  Matching scale in GeV.
     */
    void set_Q_match(double Q_match) {this->Q_match = Q_match;}

    /** Use an invocation-local parameter CSV instead of the shared default. */
    void set_param_file(std::filesystem::path path) { param_file_ = std::move(path); }

    /** Use an invocation-local Wilson CSV instead of the shared default. */
    void set_output_file(std::filesystem::path path) { output_file_ = std::move(path); }

    /// Matching scale (default: 80.379 GeV, the physical W mass).
    double Q_match = 80.379;

private:
    std::filesystem::path param_file_{};
    std::filesystem::path output_file_{};

};

#endif