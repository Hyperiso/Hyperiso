#include "HyperisoMaster.h"
#include "ParameterProvider.h"
#include "Include.h"
#include "Logger.h"
#include "CompositeParamCreator.h"
#include "QCDProvider.h"
#include "ParameterSetter.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto hyp = HyperisoMaster(); // Create the interface for hyperiso.

    HyperisoConfig config_hyp; // Config struct where we can put all the options we want for Hyperiso (general options)

    config_hyp.flags[ExternalFlag::IS_LHA_SPECTRUM] = true; // For SUSY and THDM model, tells Hyperiso that the lha already contain all the spectrum and does not need SoftSusy or 2HDMC to calculate it.
    config_hyp.flags[ExternalFlag::HAS_WILSON_INPUT] = false; // Tell Hyperiso that the lha contains the wilson coefficients as input. If only BSM coefficients are provided, SM will be calculated within Hyperiso.
    config_hyp.flags[ExternalFlag::HYP_AS_SM_MARTY] = false; // If true, Hyperiso will calculate all the wilson coefficients for the SM (up to NNLO) instead of Marty. To be used for better precision.
    config_hyp.flags[ExternalFlag::HAS_TH_OBSERVABLE_INPUT] = false; // Tell Hyperiso that the lha contains the theoretical observables as input (not used in hyperiso for the moment)

    config_hyp.model = Model::SM; // The model we want to use, SM by default. If not THDM or SUSY, MARTY is needed.

    config_hyp.mty_model_name = "ZPrime"; // Only if Config.model = Model::MARTY, name of the bsm model.
    config_hyp.mty_model_path = "/home/theo/hyperiso/Assets/input_files/marty_model/ZPrime.h"; // Only if Config.model = Model::MARTY, path of the bsm model.

    hyp.init("lha/si_input.flha", config_hyp); // Initialize program manager with LHA file and the config. Search in the Assets directory if relative path.

    ParameterProvider sm {ParameterType::SM};
    ParameterProvider wil {ParameterType::WILSON};
    LOG_INFO(sm("SMINPUTS", 6));

    CompositeParamAdapter cpc;
    std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"SMINPUTS", "MASS"}}};

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        QCDProvider qcd;
        auto xt = std::pow(qcd(MassConfig{6, 80, MassType::MSBAR, MassType::POLE}) / src.at("MASS")->retrieve(24)->get_val(), 2);
        dep_block->store_or_assign(1, std::make_shared<Parameter>(Parameter({ParameterType::WILSON, "WPARAM", 1}, xt, 0., 0.)));
    };
    cpc.add_block_dependency("WPARAM", src, ParameterType::WILSON, func);

    ParameterSetter ps;
    LOG_INFO("Before: m_W =", sm("MASS", 24), ", x_t =", wil("WPARAM", 1));
    ps.mutate({ParameterType::SM, "MASS", 24}, 100);
    LOG_INFO("After: m_W =", sm("MASS", 24), ", x_t =", wil("WPARAM", 1));

    return 0;
}