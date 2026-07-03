#include "SMParamSetter.h"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <utility>

namespace {

double sqr(double x) { return x * x; }
double fourth(double x) { const double x2 = x * x; return x2 * x2; }

} // namespace

SMParamSetter::SMParamSetter(const std::string& model,
                             std::set<std::string> special_blocks,
                             std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> sm_proxy,
                             std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> bsm_proxy,
                             std::string cinematic_template)
    : special_blocks(std::move(special_blocks)),
      sm_proxy(std::move(sm_proxy)),
      bsm_proxy(std::move(bsm_proxy)),
      cinematic_template(std::move(cinematic_template))
{
    if (model == "SM") {
        this->model_type = Model::SM;
    } else if (model == "THDM") {
        this->model_type = Model::THDM;
    } else if (model == "MSSM" || model == "NMSSM") {
        this->model_type = Model::SUSY;
    } else {
        this->model_type = Model::MARTY;
    }

    if (!this->cinematic_template.empty()) {
        const auto process = CinematicExtractor().extract_process(this->cinematic_template);
        if (!process.empty()) {
            this->cinematic_process = process;
        }
    }
}

std::unordered_map<std::string, double> SMParamSetter::setParam(const std::string& name, const InterpretedParam& interpretedParam) {

    std::unordered_map<std::string, double> params {};

    LOG_DEBUG("setting parameter", name, interpretedParam.block, interpretedParam.code);
    std::set<std::string> special = this->special_blocks;
    if (special.find(interpretedParam.block) != special.end()) {
        params[name] = calculateValue(interpretedParam);
    } else if (interpretedParam.block == "MASS" && (interpretedParam.code == LhaID(5) || interpretedParam.code == LhaID(6))) {
        if (interpretedParam.code == LhaID(5)) {
            params[name] = (*sm_proxy)("MASS_EW_SCALE", LhaID(5, 1));
        } else {
            params[name] = (*sm_proxy)("MASS_EW_SCALE", 6);
        }
    } else if (interpretedParam.block == "GAUGE" && interpretedParam.code == LhaID(4)) {
        params[name] = calculateValue(interpretedParam);
    } else {
        if (interpretedParam.is_bsm) {
            if (interpretedParam.is_complex) {
                params[name+ "_rel"] = (*bsm_proxy)(interpretedParam.block, interpretedParam.code).real();
                params[name + "_img"] = (*bsm_proxy)(interpretedParam.block, interpretedParam.code).imag();
            } else {
                params[name] = (*bsm_proxy)(interpretedParam.block, interpretedParam.code);
            }
        } else {
            if (interpretedParam.is_complex) {
                params[name+ "_rel"] = (*sm_proxy)(interpretedParam.block, interpretedParam.code).real();
                params[name + "_img"] = (*sm_proxy)(interpretedParam.block, interpretedParam.code).imag();
            } else {
                params[name] = (*sm_proxy)(interpretedParam.block, interpretedParam.code);
            }
        }
    }
    return params;
}

scalar_t SMParamSetter::calculateValue(const InterpretedParam& interpretedParam) {
    if (interpretedParam.block == "KIN") {
        return calculateKinematicInvariant(interpretedParam.code);
    }
    if (interpretedParam.block == "WEIN") {
        return asin(sqrt((*sm_proxy)("SMINPUTS", LhaID(7, 1))));
    }
    if (interpretedParam.block == "REGPROP") {
        return 1e-6;
    }
    if (interpretedParam.block == "BETA") {
        return atan((*bsm_proxy)("MINPAR", 3));
    }
    if(interpretedParam.block == "GAUGE") {
        if (interpretedParam.code == LhaID(4)) {
            return std::sqrt((*sm_proxy)("SMINPUTS", 2) * std::sqrt(2))
                 * (*sm_proxy)("SMINPUTS", 4)
                 * std::sin(2 * asin(sqrt((*sm_proxy)("SMINPUTS", LhaID(7, 1)))));
        }
    }
    return 1.0;
}

scalar_t SMParamSetter::calculateKinematicInvariant(const LhaID& code) const {
    if (!cinematic_process.has_value()) {
        return legacyKinematicInvariant(code);
    }

    const auto masses = extractMassesForCurrentProcess();

    if (cinematic_process->incoming_count() == 1 && cinematic_process->outgoing_count() == 3) {
        return calculateOneToThreeInvariant(code, masses);
    }

    if (cinematic_process->incoming_count() == 1 && cinematic_process->outgoing_count() == 2) {
        return calculateOneToTwoInvariant(code, masses);
    }

    LOG_WARN("SMParamSetter", "Unsupported MARTY kinematics. Falling back to legacy KIN rule.",
             cinematic_process->incoming_count(), "incoming and", cinematic_process->outgoing_count(), "outgoing particles.");
    return legacyKinematicInvariant(code);
}

scalar_t SMParamSetter::calculateOneToThreeInvariant(const LhaID& code, const std::vector<scalar_t>& masses) const {
    if (masses.size() != 4) {
        return legacyKinematicInvariant(code);
    }

    const double m1 = masses[0];
    const double m2 = masses[1];
    const double m3 = masses[2];
    const double m4 = masses[3];
    const double denom = m1 - m2;
    const double denom2 = sqr(denom);

    if (std::abs(denom) < 1e-15) {
        LOG_WARN("SMParamSetter", "Singular 1->3 KIN denominator m1-m2. Falling back to legacy KIN rule.");
        return legacyKinematicInvariant(code);
    }

    if (code == LhaID(12)) {
        return m1 * m2;
    }
    if (code == LhaID(13)) {
        return m1 * (sqr(m1) - 2*m1*m2 + sqr(m2) + sqr(m3) - sqr(m4)) / (2 * denom);
    }
    if (code == LhaID(14)) {
        return m1 * (sqr(m1) - 2*m1*m2 + sqr(m2) - sqr(m3) + sqr(m4)) / (2 * denom);
    }
    if (code == LhaID(23)) {
        return m2 * (sqr(m1) - 2*m1*m2 + sqr(m2) + sqr(m3) - sqr(m4)) / (2 * denom);
    }
    if (code == LhaID(24)) {
        return m2 * (sqr(m1) - 2*m1*m2 + sqr(m2) - sqr(m3) + sqr(m4)) / (2 * denom);
    }
    if (code == LhaID(34)) {
        return (
            sqr(m1)*sqr(m3) + sqr(m1)*sqr(m4)
          - 2*m1*m2*sqr(m3) - 2*m1*m2*sqr(m4)
          + sqr(m2)*sqr(m3) + sqr(m2)*sqr(m4)
          - fourth(m3) + 2*sqr(m3)*sqr(m4) - fourth(m4)
        ) / (2 * denom2);
    }

    return legacyKinematicInvariant(code);
}

scalar_t SMParamSetter::calculateOneToTwoInvariant(const LhaID& code, const std::vector<scalar_t>& masses) const {
    if (masses.size() != 3) {
        return legacyKinematicInvariant(code);
    }

    const double m1 = masses[0];
    const double m2 = masses[1];
    const double m3 = masses[2];

    if (std::abs(m1) < 1e-15) {
        LOG_WARN("SMParamSetter", "Singular 1->2 KIN denominator m1. Falling back to legacy KIN rule.");
        return legacyKinematicInvariant(code);
    }

    const double delta23 = sqr(m2) - sqr(m3);
    const double root12_arg = fourth(m1) + 2*sqr(m1)*sqr(m2) - 2*sqr(m1)*sqr(m3) + sqr(delta23);
    const double root13_arg = fourth(m1) - 2*sqr(m1)*sqr(m2) + 2*sqr(m1)*sqr(m3) + sqr(delta23);

    const double root12 = std::sqrt(std::max(0.0, root12_arg));
    const double root13 = std::sqrt(std::max(0.0, root13_arg));

    if (code == LhaID(12)) {
        return root12 / 2.0;
    }
    if (code == LhaID(13)) {
        return root13 / 2.0;
    }
    if (code == LhaID(23)) {
        return sqr(m1)/4.0 - sqr(m2)/2.0 - sqr(m3)/2.0
             + sqr(delta23)/(4.0*sqr(m1))
             + root13*root12/(4.0*sqr(m1));
    }

    return legacyKinematicInvariant(code);
}

scalar_t SMParamSetter::legacyKinematicInvariant(const LhaID& code) const {
    if (code == LhaID(34)) {
        return -pow((*sm_proxy)("MASS", 13), 2.);
    }
    return (pow((*sm_proxy)("MASS_EW_SCALE", LhaID(5, 1)) ,2.) + std::pow((*sm_proxy)("MASS", 3), 2.))/2.;
}

std::vector<scalar_t> SMParamSetter::extractMassesForCurrentProcess() const {
    std::vector<scalar_t> masses;
    if (!cinematic_process.has_value()) {
        return masses;
    }

    const auto particles = cinematic_process->ordered_particles();
    masses.reserve(particles.size());

    for (const auto& particle : particles) {
        masses.push_back(massValueForParticle(particle));
    }
    return masses;
}

scalar_t SMParamSetter::massValueForParticle(const std::string& particle_name) const {
    const auto p = CinematicExtractor::normalize_particle_name(particle_name);

    if (p == "a" || p == "g" || p == "ve" || p == "vmu" || p == "vtau") {
        return 0.0;
    }

    if (p == "b") {
        return (*sm_proxy)("MASS_EW_SCALE", LhaID(5, 1));
    }
    if (p == "t") {
        return (*sm_proxy)("MASS_EW_SCALE", 6);
    }

    if (p == "d") return (*sm_proxy)("MASS", 1);
    if (p == "u") return (*sm_proxy)("MASS", 2);
    if (p == "s") return (*sm_proxy)("MASS", 3);
    if (p == "c") return (*sm_proxy)("MASS", 4);

    if (p == "e")   return (*sm_proxy)("MASS", 11);
    if (p == "mu")  return (*sm_proxy)("MASS", 13);
    if (p == "tau") return (*sm_proxy)("MASS", 15);

    if (p == "z") return (*sm_proxy)("MASS", 23);
    if (p == "w" || p == "w+" || p == "w-") return (*sm_proxy)("MASS", 24);
    if (p == "h") return (*sm_proxy)("MASS", 25);

    LOG_WARN("SMParamSetter", "Unknown particle mass for MARTY kinematics:", particle_name,
             ". Treating it as massless. Add it to SMParamSetter::massValueForParticle if needed.");
    return 0.0;
}
