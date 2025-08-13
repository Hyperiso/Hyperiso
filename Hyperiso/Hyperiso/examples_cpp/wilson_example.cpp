#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"

int main() {
    auto hyp = HyperisoMaster(); // Create the interface for hyperiso (reading/writing in the lha and all the options we want to use)

    Config config_hyp; // Config struct where we can put all the options we want for Hyperiso (general options)

    config_hyp.flags[ExternalFlag::USE_MARTY] = false; // Rather we want to use MARTY (mandatory for new models) or hardcoded function (for THDM, SUSY and SM)

    config_hyp.model = Model::SM; // The model we want to use, SM by default. If not THDM or SUSY, MARTY is needed.

    config_hyp.mty_model_name = "ZPrime"; // Only if Config.model = Model::Custom, name of the bsm model.
    config_hyp.mty_model_path = "/home/theo/hyperiso/Assets/input_files/marty_model/ZPrime.h"; // Only if Config.model = Model::Custom, path of the bsm model.

    hyp.init("lha/camilia.flha", config_hyp); // Initialize program manager with LHA file and the config. Search in the Assets directory if relative path.

    auto wi = WilsonInterface(); // Initialize interface for wilson coefficient calculation.

    WilsonBuildConfig config({WGroup::B}, 81, 42, QCDOrder::LO); // Configuration for the wilsons coefficients. First argument is the name of the group, second is the matching scale, third is the running scale and finally the last is the order for the coefficient.

    wi.build(config); // build all the coefficients using the given config.

    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD"); // Log a block from the HyperisoMaster (lha or created by the wilsons).

    LOG_INFO("Parameters created");
    LOG_INFO("C9(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::LO, ContributionType::TOTAL)); //getR -> get running coefficient (getM -> get matching coefficient)
    // ContributionType can be SM, BSM or TOTAL (SM+BSM)
    // If you want to add all QCDOrder, use getFR (or get FM for matching)
    
    LOG_INFO("AGAIN");
    LOG_INFO("C9(mu_h) at NLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NLO, ContributionType::TOTAL));
    LOG_INFO("C9(mu_h) at NNLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    
    LOG_INFO("C9(mu_h) full =", wi.getFR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    return 0;
}