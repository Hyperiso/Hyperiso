#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"
#include "BlockProxy.h"

int main() {
    auto hyp = HyperisoMaster();

    HyperisoConfig config_hyp; 

    config_hyp.model = Model::MARTY; // We specify here that we want to use MARTY to calculate the wilson coefficients
    //Be sure to have MARTY installed, if you installed it within Hyperiso, then you don't have anything to do
    //If you have MARTY on your machine please use the following line:
    // hyp.pre_init_set_marty_path("/path/to/marty");

    //If you have in your LHA new blocks which are not in the flha convention, please add this line :
    // hyp.pre_init_add_block("NAME_OF_THE_BLOCK");
    //By default block are id + value. if you have another type of block (with 2 or more id) please use the options of the pre_init_add_block API.

    config_hyp.flags[ExternalFlag::HYP_AS_SM_MARTY] = true; //If true, then Hyperiso will use the hardcoded version of wilson coefficients for their SM part. MARTY will only be used for the BSM part.
    //This is helpful to get more precision, with the SM part up to NNLO within Hyperiso.

    config_hyp.mty_model_name = "ZPrime"; //Use the name of the class you used in your model file (here Zprime.h) within MARTY.
    config_hyp.mty_model_path = "/home/theo/hyperiso/Assets/input_files/marty_model/ZPrime.h"; //Path to the MARTY model file. If relative path is used, the path start from the Assets folder.
    
    hyp.init("lha/zprime_input.flha", config_hyp);

    //The rest of the code is the same than the base_wilson_example.cpp example.
    auto wi = WilsonInterface();

    WilsonBuildConfig config({WGroup::B}, 81, 42, QCDOrder::NLO);

    wi.build(config);

    LOG_INFO("C7(mu_h) at LO =", wi.getM(WGroup::B, WCoef::C7, QCDOrder::LO, ContributionType::TOTAL));
    LOG_INFO("C9(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::LO, ContributionType::TOTAL)); //getR -> get running coefficient (getM -> get matching coefficient)

    LOG_INFO("C9(mu_h) at NLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NLO, ContributionType::TOTAL));
    LOG_INFO("C9(mu_h) at NNLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    
    LOG_INFO("C9(mu_h) full =", wi.getFR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    return 0;
}