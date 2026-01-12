#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"
#include "BlockProxy.h"

int main() {
    auto hyp = HyperisoMaster();
    HyperisoConfig config_hyp;
    // config_hyp.flags[ExternalFlag::USE_MARTY] = false;
    config_hyp.model = Model::THDM;
    config_hyp.mty_model_name = "ZPrime";
    config_hyp.mty_model_path = "/home/theo/hyperiso/Assets/input_files/marty_model/ZPrime.h";

    hyp.init("lha/testinput_thdm.lha", config_hyp); // Initialize program manager with LHA file

    auto wi = WilsonInterface(); // Initialize interface and build the required groups

    WilsonBuildConfig config({WGroup::B}, 81, 42, QCDOrder::LO);

    wi.build(config);

    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD");

    LOG_INFO("Parameters created");
    LOG_INFO("C9(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::LO, ContributionType::TOTAL));
    LOG_INFO("AGAIN");
    LOG_INFO("C9(mu_h) at NLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NLO, ContributionType::TOTAL));
    LOG_INFO("C9(mu_h) at NNLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    
    LOG_INFO("C9(mu_h) full =", wi.getFR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));

    HyperisoConfig config_v2;

    BlockProxy().log_all_blocks(ParameterType::WILSON);
    // LOG_INFO("changing config");
    hyp.switch_lha("default/lha/testInput.flha", config_v2);

    BlockProxy().log_all_blocks(ParameterType::WILSON);

    // WilsonBuildConfig config2({WGroup::B}, 81, 5, QCDOrder::LO);

    // wi.build(config2);

    // std::shared_ptr<Block> a = std::make_shared<Block>();
    // std::shared_ptr<Block> b = std::make_shared<Block>();

    // a->blockname = "A";
    // b->blockname = "B";

    // auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
    //     // src.at("MASS")->retrieve(25)->get_val();
    //     dep_block->store_or_assign(0, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", 1}, 42, 0., 0.));
    // };
    // std::unordered_map<std::string, std::shared_ptr<Block>> truc = {};

    // truc["A"] = a;
    // truc["B"] = b;

    // a->store(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", 1}, 5, 0., 0.));
    // a->store(3, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", 3}, 10, 0., 0.));

    // b->store(11, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "BBB", 11}, 2929, 0., 0.));

    // std::shared_ptr<Block> c = std::make_shared<DependentBlock>(truc, func);
    // c->blockname = "C";
    // c->init();
    // c->update();

    // auto func_d = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
    //     // src.at("MASS")->retrieve(25)->get_val();
    //     dep_block->store_or_assign(5, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "AAAAAH", 5}, 666, 0., 0.));
    // };
    // std::unordered_map<std::string, std::shared_ptr<Block>> truc_d = {};

    // truc_d["B"] = b;
    // truc_d["C"] = c;

    // std::shared_ptr<Block> d = std::make_shared<DependentBlock>(truc_d, func_d);
    // d->blockname = "D";
    // d->init();
    // d->update();


    // auto func_e = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
    //     // src.at("MASS")->retrieve(25)->get_val();
    //     dep_block->store_or_assign(55, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EAEAER", 55}, 10000, 0., 0.));
    // };
    // std::unordered_map<std::string, std::shared_ptr<Block>> truc_e = {};

    // // truc_e["B"] = b;
    // truc_e["D"] = d;

    // std::shared_ptr<Block> e = std::make_shared<DependentBlock>(truc_e, func_e);
    // e->blockname = "E";
    // e->init();
    // e->update();


    // auto obs = d->getObservers();
    // auto obs_b = b->getObservers();
    // for (auto ob : obs) {
    //     std::cout << "truc : " << ob->blockname << std::endl;
    // }

    // for (auto ob : obs_b) {
    //     std::cout << "truc_b : " << ob->blockname << std::endl;
    // }

    // std::cout << "a :" << std::endl;
    // std::cout << a << std::endl;

    // std::cout << "b :" << std::endl;
    // std::cout << b << std::endl;

    // std::cout << "c :" << std::endl;
    // std::cout << c << std::endl;

    // std::cout << "d :" << std::endl;
    // std::cout << d << std::endl;

    // std::cout << "e :" << std::endl;
    // std::cout << e << std::endl;
    // c->destroy();

    // std::cout << "a :" << std::endl;
    // std::cout << a << std::endl;

    // std::cout << "b :" << std::endl;
    // std::cout << b << std::endl;

    // std::cout << "c :" << std::endl;
    // std::cout << c << std::endl;

    // std::cout << "d :" << std::endl;
    // std::cout << d << std::endl;

    // std::cout << "e :" << std::endl;
    // std::cout << e << std::endl;
    return 0;
}