#include "HyperisoMaster.h"
#include "ParameterProvider.h"
#include "Include.h"
#include "Logger.h"
#include "CompositeParamCreator.h"
#include "QCDProvider.h"
#include "ParameterSetter.h"
#include "DependencyPruner.h"
#include "DependantBlockInfoProvider.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    auto hyp = HyperisoMaster();

    HyperisoConfig config_hyp;

    config_hyp.model = Model::SM;

    hyp.init("lha/si_input.flha", config_hyp);

    ParameterProvider sm;
    DependencyPruner dp;
    DependantBlockInfoProvider dbip;

    std::cout << dbip.is_dependent_block(ParameterType::SM, "VCKM") << std::endl;

    dp.detach_block(ParameterType::SM, "VCKM");
    
    std::cout << dbip.is_dependent_block(ParameterType::SM, "VCKM") << std::endl;


    return 0;
}