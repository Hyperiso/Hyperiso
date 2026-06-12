#ifndef GPP_COMPILER_STRATEGY_H
#define GPP_COMPILER_STRATEGY_H

#include <cstdlib>

#include "FileNameManager.h"
#include "CompilerStrategy.h"
#include "config.hpp"


/**
 * @file GppCompilerStrategy.h
 * @brief Declares a compilation strategy using g++ directly.
 *
 * This header defines ::GppCompilerStrategy, a concrete implementation of
 * ::CompilerStrategy that:
 *  - invokes `g++` with the appropriate MARTY libraries / rpath,
 *  - then executes the resulting binary in the correct output directory.
 */

 /**
 * @class GppCompilerStrategy
 * @ingroup CodeGenerationModule
 * @brief Compilation strategy using the g++ compiler.
 *
 * GppCompilerStrategy:
 *  - checks whether a binary already exists (via ::check_if_compile()),
 *  - if needed, compiles the provided source file into an executable,
 *  - finally runs the executable from the output directory managed by
 *    ::FileNameManager.
 */
class GppCompilerStrategy : public CompilerStrategy {
public:
    /**
     * @brief Constructs a GppCompilerStrategy for a given (model, wilson) pair.
     *
     * @param model   Model name (used only as metadata; compilation flags
     *                are global).
     * @param wilson  Wilson basis / label (used to locate output directories).
     */
    GppCompilerStrategy(std::string model, std::string wilson) : CompilerStrategy(model, wilson) {}

    /**
     * @brief Compiles if necessary and then runs the resulting binary.
     *
     * If ::check_if_compile() returns false, ::compile() is called first.
     * The binary is then executed from the directory returned by
     * ::FileNameManager::getOutputDir().
     *
     * @param sourceFile    Path to the generated C++ source file.
     * @param outputBinary  Name/path of the executable to create and run.
     */
    void compile_run(const std::string& sourceFile, const std::string& outputBinary) override;

    /**
     * @brief Compiles the given source file using g++ and links against MARTY.
     *
     * The command line typically includes:
     *  - `-L<MARTY_INSTALL>/lib`
     *  - `-Wl,-rpath,<MARTY_INSTALL>/lib`
     *  - `-lmarty`
     *  - the absolute `libgfortran.so` path discovered from g++/gfortran
     *
     * @param sourceFile    Path to the generated C++ source.
     * @param outputBinary  Path/name of the executable to produce.
     */
    void compile(const std::string& sourceFile, const std::string& outputBinary) override;
};

#endif
