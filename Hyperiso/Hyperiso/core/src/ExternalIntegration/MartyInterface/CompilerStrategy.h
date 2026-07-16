#ifndef COMPILER_STRATEGY_H
#define COMPILER_STRATEGY_H

#include <string>
#include <fstream>
#include <sys/stat.h>
#include <iostream>
#include <cstdio>
#include <memory>
#include <array>
#include <stdexcept>

#include "Include.h"

/**
 * @file CompilerStrategy.h
 * @brief Declares compilation strategies and utilities for executing shell commands.
 *
 * This header defines:
 *  - ::CompilerStrategy: an abstract base class for different compilation backends
 *    (e.g. g++, make-based).
 *  - ::executeCommand(): a small helper to run shell commands and capture errors.
 *
 * These are mainly used in the MARTY interface layer to build and run generated
 * C++ code for Wilson coefficients and matching routines.
 */

/**
 * @brief Executes a shell command and captures its output.
 *
 * The @p command string is executed via `popen()`, with stderr redirected to
 * stdout (using `2>&1`). The complete output is read and, on failure, printed
 * to `std::cerr`.
 *
 * @param command Shell command to execute.
 * @return `true` if the command exited with status 0.
 *
 * @throws std::runtime_error if the command cannot be started or exits with
 *         a non-zero status. This intentionally stops the MARTY pipeline at
 *         the first failing command instead of cascading into missing-file
 *         errors.
 */
bool executeCommand(const std::string& command);

/**
 * @class CompilerStrategy
 * @ingroup CodeGenerationModule
 * @brief Abstract base class for compilation and run strategies.
 *
 * CompilerStrategy encapsulates how generated C++ code is compiled and/or
 * executed. Concrete strategies can:
 *
 *  - call `g++` directly (see ::GppCompilerStrategy),
 *  - invoke a `make`-based build system (see ::MakeCompilerStrategy),
 *  - implement custom logic for different toolchains.
 *
 * Typical usage:
 *  - construct a strategy with a given (model, wilson) pair,
 *  - call ::compile_run() on the generated source and target binary.
 */
class CompilerStrategy {
public:
    /**
     * @brief Constructs a CompilerStrategy with a given model and Wilson label.
     *
     * @param model   Name of the physics model (e.g. "SM", "SUSY").
     * @param wilson  Name of the Wilson basis or identifier.
     */
    CompilerStrategy(std::string model, std::string wilson) : model(model), wilson(wilson) {}

    /// Virtual destructor.
    virtual ~CompilerStrategy() = default;

    /**
     * @brief Compiles (if needed) and then runs the resulting binary.
     *
     * The typical implementation:
     *  - calls ::check_if_compile(),
     *  - if not already compiled, invokes ::compile(),
     *  - finally executes the resulting binary.
     *
     * @param sourceFile    Path to the generated C++ source.
     * @param outputBinary  Path/name of the binary to produce and run.
     */
    virtual void compile_run(const std::string& sourceFile, const std::string& outputBinary) = 0;

    /**
     * @brief Compiles the given source file into an output binary.
     *
     * Each concrete strategy is responsible for building the appropriate
     * compiler command (e.g. `g++`, `make`, etc.) and calling ::executeCommand().
     *
     * @param sourceFile    Path to the generated C++ source.
     * @param outputBinary  Path/name of the binary to produce.
     */
    virtual void compile(const std::string& sourceFile, const std::string& outputBinary) = 0;

    /**
     * @brief Checks whether a compiled binary already exists and is non-empty.
     *
     * This is used to avoid recompiling unchanged binaries on repeated calls.
     *
     * @param outputBinary  Path to the binary.
     * @return `true` if the file exists and has non-zero size, `false` otherwise.
     */
    virtual bool check_if_compile(const std::string& outputBinary);

protected:
    /// Model name.
    std::string model;
    /// Wilson basis / label name.
    std::string wilson;
};

#endif
