#ifndef SMPARAMSETTER_H
#define SMPARAMSETTER_H

#include <cmath>
#include <set>
#include <string>
#include <iostream>
#include <vector>
#include <optional>

#include "config.hpp"
#include "LhaID.h"
#include "Interpreter.h"
#include "IMartyParameterProxy.h"
#include "CinematicExtractor.h"

/**
 * @class SMParamSetter
 * @brief Converts Hyperiso parameters into MARTY input parameters.
 */
class SMParamSetter {
public:
    /**
     * @brief Constructs a parameter setter for a given model.
     *
     * @param model                  Model name, e.g. "SM" or "ZPrime".
     * @param special_blocks         Blocks treated with special formulas.
     * @param sm_proxy               Proxy used to access SM parameters.
     * @param bsm_proxy              Optional proxy used to access BSM parameters.
     * @param cinematic_template     Optional generated/template C++ file containing
     *                               computeWilsonCoefficients(...). When provided,
     *                               KIN invariants are inferred from the process legs.
     */
    SMParamSetter(const std::string& model,
                  std::set<std::string> special_blocks,
                  std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> sm_proxy,
                  std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> bsm_proxy = nullptr,
                  std::string cinematic_template = "");

    std::unordered_map<std::string, double> setParam(const std::string& name, const InterpretedParam& interpretedParam);

private:
    scalar_t calculateValue(const InterpretedParam& interpretedParam);

    scalar_t calculateKinematicInvariant(const LhaID& code) const;
    scalar_t calculateOneToThreeInvariant(const LhaID& code, const std::vector<scalar_t>& masses) const;
    scalar_t calculateOneToTwoInvariant(const LhaID& code, const std::vector<scalar_t>& masses) const;
    scalar_t legacyKinematicInvariant(const LhaID& code) const;

    std::vector<scalar_t> extractMassesForCurrentProcess() const;
    scalar_t massValueForParticle(const std::string& particle_name) const;

    Model model_type;
    std::set<std::string> special_blocks;

    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> sm_proxy;
    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> bsm_proxy;

    std::string cinematic_template;
    std::optional<CinematicProcess> cinematic_process;
};

#endif // SMPARAMSETTER_H
