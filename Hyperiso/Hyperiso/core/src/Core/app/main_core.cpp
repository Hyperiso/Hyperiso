#include "HyperisoMaster.h"
#include "ParameterProvider.h"
#include "Include.h"
#include "Logger.h"
#include "CompositeParamCreator.h"
#include "QCDProvider.h"
#include "ParameterSetter.h"
#include "BlockProvider.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    Config config;
    config.model = Model::SM;

    config.flags[ExternalFlag::HAS_WILSON_INPUT] = true;

    HyperisoMaster hi;
    hi.init("lha/testInput.flha", config);

    BlockProvider().log_all_blocks(ParameterType::WILSON);
    
    ParameterProvider sm {ParameterType::SM};
    ParameterProvider wil {ParameterType::WILSON};
    ParameterProvider obs {ParameterType::OBSERVABLE};

    auto obs_p = obs.get_parameter({ParameterType::OBSERVABLE, "FOBS", {521, 2, 3, 321, 13}});
    LOG_INFO(*obs_p);
    auto bin = obs_p->get_bin();
    LOG_INFO("A_FB(B > K mu mu) = [", bin.first, ",", bin.second, "] =", obs_p->get_val());
    

    // 521_11_3_423_-15_16

    auto obs_p2 = obs.get_parameter({ParameterType::OBSERVABLE, "FOBS", {521, 11, 3, 423, -15, 16}});
    LOG_INFO(*obs_p2);
    auto bin2 = obs_p2->get_bin();
    LOG_INFO("A_FB(B > K mu mu) = [", bin2.first, ",", bin2.second, "] =", obs_p2->get_val());

    LOG_INFO(sm("SMINPUTS", 6));

    CompositeParamAdapter cpc;
    std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"SMINPUTS", "MASS"}}};

    auto func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
        QCDProvider qcd;
        auto xt = std::pow(qcd(MassConfig{6, 80, MassType::MSBAR, MassType::POLE}) / src.get_val("MASS",24), 2);
        dep_block->store_or_assign(1, std::make_shared<Parameter>(Parameter({ParameterType::WILSON, "WPARAM", 1}, xt, 0., 0.)));
    };
    cpc.add_block_dependency("WPARAM", src, ParameterType::WILSON, func);

    ParameterSetter ps;
    LOG_INFO("Before: m_W =", sm("MASS", 24), ", x_t =", wil("WPARAM", 1));
    ps.mutate({ParameterType::SM, "MASS", 24}, 100);
    LOG_INFO("After: m_W =", sm("MASS", 24), ", x_t =", wil("WPARAM", 1));

    

    return 0;
}