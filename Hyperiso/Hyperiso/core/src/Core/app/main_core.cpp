#include "HyperisoMaster.h"
#include "ParameterProvider.h"
#include "General.h"
#include "Logger.h"
#include "CompositeParamCreator.h"
#include "QCDProvider.h"
#include "ParameterSetter.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    HyperisoMaster hi;
    hi.init("lha/testInput.flha");

    ParameterProvider sm {ParameterType::SM};
    ParameterProvider wil {ParameterType::WILSON};
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