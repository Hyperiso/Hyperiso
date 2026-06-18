#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

#include "CompositeParamCreator.h"
#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ParameterProvider.h"
#include "ParameterSetter.h"
#include "QCDProvider.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // Dev example: create a dependent block from existing source blocks.
    // This is useful for derived quantities which should update lazily when an input changes.
    HyperisoConfig config;
    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    ParameterProvider sm(ParameterType::SM);
    ParameterProvider wil(ParameterType::WILSON);

    CompositeParamAdapter composite_creator;

    // The new WILSON::WPARAM block depends on SM::SMINPUTS and SM::MASS.
    std::unordered_map<ParameterType, std::vector<std::string>> sources = {
        {ParameterType::SM, {"SMINPUTS", "MASS"}}
    };

    auto update_func = [](const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
        QCDProvider qcd;

        // Example derived quantity: x_t = (m_t(mu) / m_W)^2.
        const double m_top = qcd(MassConfig(6, 80.0, MassType::MSBAR, MassType::POLE));
        const double m_w = src.get_val("MASS", LhaID(24));
        const double x_t = std::pow(m_top / m_w, 2.0);

        dep_block->store_or_assign(
            LhaID(1),
            std::make_shared<Parameter>(Parameter(ParamId(ParameterType::WILSON, "WPARAM", LhaID(1)), x_t, 0.0, 0.0))
        );
    };

    composite_creator.add_block_dependency("WPARAM", sources, ParameterType::WILSON, update_func);

    ParameterSetter setter;

    std::cout << "Before: m_W = " << sm("MASS", LhaID(24))
              << ", x_t = " << wil("WPARAM", LhaID(1)) << "\n";

    // Because WPARAM is a dependent block, the value changes when its source changes.
    setter.mutate(ParamId(ParameterType::SM, "MASS", LhaID(24)), 100.0);

    std::cout << "After: m_W = " << sm("MASS", LhaID(24))
              << ", x_t = " << wil("WPARAM", LhaID(1)) << "\n";

    return 0;
}
