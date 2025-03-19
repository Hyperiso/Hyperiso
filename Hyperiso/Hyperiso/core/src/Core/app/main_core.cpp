#include "HyperisoMaster.h"
#include "ParameterProvider.h"
#include "General.h"
#include "Logger.h"
#include "CompositeParamCreator.h"
#include "QCDProvider.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    HyperisoMaster hi;
    hi.init("lha/testInput.flha");

    ParameterProvider provider {ParameterType::SM};
    ParameterProvider wparam_provider;
    LOG_INFO(provider("SMINPUTS", 6));
    
    ParamId xt_pid {ParameterType::WILSON, "SMINPUTS", 6};

    CompositeParamCreator cpc;
    std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"SMINPUTS", "MASS"}}};
    auto func = [xt_pid] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        QCDProvider qcd;
        auto xt = std::pow(qcd(MassConfig{6, 80, MassType::POLE, MassType::POLE}) / src.at("MASS")->retrieve(24).get_val(), 2);
        dep_block->store(1, Parameter(xt_pid, xt, 0., 0.));
    };

    cpc.add_dependency("WPARAM", src, ParameterType::WILSON, func);
    LOG_INFO(wparam_provider(xt_pid));

    return 0;
}