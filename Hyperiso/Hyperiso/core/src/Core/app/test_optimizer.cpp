#include "Parameters.h"
#include "Include.h"
#include "ParamOptimizer.h"
#include "HyperisoMaster.h"

int main() {

    Config config;
    config.model = Model::THDM;

    HyperisoMaster hi;
    hi.init("lha/testinput_thdm.lha", config);
    auto sm  = Parameters::GetInstance(ParameterType::SM)->get_block_accessor();
    auto bsm = Parameters::GetInstance(ParameterType::BSM)->get_block_accessor();
    std::cout << "here ! " << std::endl;

    ParamOptimizer opt({sm, bsm});

    opt.set_value("YU", LhaID(6), 1.15);
    opt.set_value("YU", LhaID(6), 1.18);
    opt.set_param("GAUGE", LhaID(1), std::make_shared<Parameter>(ParamId("GAUGE", LhaID(1)), 0.357, 0., 0.));
    opt.remove("MASS", LhaID(24));

    opt.commit();

}