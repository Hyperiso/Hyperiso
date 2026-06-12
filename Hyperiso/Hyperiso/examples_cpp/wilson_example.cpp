#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"
#include "BlockProxy.h"

int main() {
    auto hyp = HyperisoMaster(); // Create the interface for hyperiso.

    HyperisoConfig config_hyp; // Config struct where we can put all the options we want for Hyperiso (general options)

    config_hyp.model = Model::SM; // The model we want to use, SM by default. If not THDM or SUSY, MARTY is needed.

    hyp.init("lha/si_input.flha", config_hyp); // Initialize program manager with LHA file and the config. Search in the Assets directory if relative path.

    auto wi = WilsonInterface(); // Initialize interface for wilson coefficient calculation.

    WilsonBuildConfig config({WGroup::B}, 81, 42, QCDOrder::LO); // Configuration for the wilsons coefficients. First argument is the name of the group, second is the matching scale, third is the running scale and finally the last is the order for the coefficient.

    wi.build(config); // build all the coefficients using the given config.

    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD"); // Log a block from the HyperisoMaster (lha or created by the wilsons).

    LOG_INFO("C9(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::LO, ContributionType::TOTAL)); //getR -> get running coefficient (getM -> get matching coefficient)
    // ContributionType can be SM, BSM or TOTAL (SM+BSM)
    // If you want to add all QCDOrder, use getFR (or get FM for matching)
    
    LOG_INFO("C9(mu_h) at NLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NLO, ContributionType::TOTAL));
    LOG_INFO("C9(mu_h) at NNLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    
    LOG_INFO("C9(mu_h) full =", wi.getFR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    return 0;
}