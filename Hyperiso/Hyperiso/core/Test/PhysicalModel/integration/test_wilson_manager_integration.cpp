#include <cassert>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "WilsonManager.h"
#include "CoefficientGroupBuilder.h"
#include "WilsonCoefficientRegistry.h"
#include "GroupDefinition.h"
#include "GroupMapper.h"
#include "wcoef_ids.hpp"
#include "WilsonGroupAdapterConfig.h"
#include "IMartyWilsonProxy.h"
#include "InterpretedParam.h"
#include "Block.h"
#include "SourcesView.h"

namespace {


struct RecordingBlockComposer : IBlockComposer {
    struct ParamCall {
        ParamId dest;
        std::unordered_set<ParamId> sources;
    };
    struct BlockCall {
        std::string block_name;
        std::unordered_map<ParameterType, std::vector<std::string>> src;
    };

    std::vector<ParamCall>  param_calls;
    std::vector<BlockCall>  block_calls;
    bool cleared_all = false;

    void compose_block(const std::string& name,
                       const std::unordered_map<ParameterType, std::vector<std::string>>& src,
                       const DepUpdateFunc&) override
    {
        block_calls.push_back({name, src});
        composed_blocks.insert(name);
    }

    void compose_parameter(const ParamId& pid,
                           const std::unordered_set<ParamId>& deps,
                           const DepParamUpdateFunc&) override
    {
        param_calls.push_back({pid, deps});
    }

    void remove_block(const std::string& name) override {
        composed_blocks.erase(name);
    }

    void update(const std::string& /*name*/) override {}

    void remove_all_composed_blocks() override {
        cleared_all = true;
        composed_blocks.clear();
    }
};

struct DummyParamProxy : IParameterProxy<std::string, LhaID> {
    scalar_t operator()(const std::string&, const LhaID&) const override {
        return 4.0 * PI;
    }
    bool exist(const std::string&, const LhaID&) const override {
        return true;
    }
    double get_scale(const std::string&) const override {
        return 1.0;
    }
};

template<typename T>
struct DummyCoreAPI : ICoreAPI<T> {
    T value;
    explicit DummyCoreAPI(const T& v) : value(v) {}
    T get() override { return value; }
};

struct DummyScaleSetter : IParamSetter<ScaleType> {
    std::vector<ScaleType> switches;
    std::vector<double>    values;

    void set(double v) override { values.push_back(v); }
    void switch_param(ScaleType s) override { switches.push_back(s); }
};

struct DummyMartyProxy : IMartyWilsonProxy<InterpretedParam> {
    void calculate(std::string, std::string, double, std::string) override {}
    std::set<std::string>  get_special_blocks() override { return {}; }
    std::unordered_set<InterpretedParam> get_dependencies(std::string) override { return {}; }
};

WilsonGroupAdapterConfig make_adapters(
    const std::shared_ptr<IParameterProxy<std::string, LhaID>>& wilson_proxy,
    const std::shared_ptr<IBlockComposer>& iblock_c,
    const std::shared_ptr<DummyCoreAPI<bool>>& use_marty)
{
    auto marty_model_name = std::make_shared<DummyCoreAPI<std::string>>(std::string{""});
    auto marty_model_path = std::make_shared<DummyCoreAPI<fs::path>>(fs::path{});
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr;

    return WilsonGroupAdapterConfig{
        wilson_proxy,
        iblock_c,
        use_marty,
        marty_model_name,
        marty_model_path,
        marty_proxy
    };
}

}

int main() {
    std::cout << "== CoefficientManager INTEGRATION ==\n";


    CoefficientRegistry reg;
    register_B(reg);

    auto iblock_c     = std::make_shared<RecordingBlockComposer>();
    auto wilson_proxy = std::make_shared<DummyParamProxy>();
    auto use_marty    = std::make_shared<DummyCoreAPI<bool>>(false);

    WilsonGroupAdapterConfig adapters =
        make_adapters(wilson_proxy, iblock_c, use_marty);

    BuildContext ctx{adapters, Model::SM, Backend::Builtin, ContributionType::SM, GroupMapper::to_id(WGroup::B)};

    CoefficientGroupBuilder builder{reg};
    auto grp = builder.build(ctx);
    assert(grp && "B group must be buildable");


    auto model_api    = std::make_shared<DummyCoreAPI<Model>>(Model::SM);
    auto scale_setter = std::make_shared<DummyScaleSetter>();
    auto has_wilson_api    = std::make_shared<DummyCoreAPI<bool>>(false);
    auto hyp_as_sm    = std::make_shared<DummyCoreAPI<bool>>(false);

    WilsonPortsConfig ports{iblock_c, wilson_proxy, use_marty, has_wilson_api, model_api, scale_setter, hyp_as_sm};
    ports.build_group = nullptr;

    std::map<std::string, std::shared_ptr<CoefficientGroup>> groups;
    std::string group_name = GroupMapper::str(ctx.group_id);
    groups.emplace(group_name, grp);

    double mu_W = 80.0;
    double mu_h = 5.0;
    std::string order = "LO";

    auto manager = CoefficientManager::Builder(
        groups,
        mu_W,
        mu_h,
        order,
        ports
    );
    assert(manager && "Builder must return a valid CoefficientManager");


    assert(manager->getModel() == ModelMapper::str(Model::SM));

    assert(scale_setter->switches.size() == 2);
    assert(scale_setter->values.size()   == 2);
    assert(scale_setter->switches[0] == ScaleType::MATCHING);
    assert(scale_setter->values[0]   == mu_W);
    assert(scale_setter->switches[1] == ScaleType::HADRONIC);
    assert(scale_setter->values[1]   == mu_h);

    auto all = manager->getGroups();
    assert(all.size() == 1);
    assert(all.contains(group_name));

    auto bases = manager->getGroupBases(ctx.group_id);
    assert(bases.contains(WilsonBasis::B_STANDARD));

    assert(!iblock_c->param_calls.empty()
           && "Builder should have defined dependent parameters for matching / combination");
    assert(!iblock_c->block_calls.empty()
           && "Builder should have composed at least one hadronic block");

    std::string expected_block =
        GroupMapper::str(GroupMapper::enum_elt(group_name),
                         ScaleType::HADRONIC,
                         WilsonBasis::B_STANDARD);

    bool found_hadronic = false;
    for (const auto& bc : iblock_c->block_calls) {
        if (bc.block_name == expected_block) {
            found_hadronic = true;
            break;
        }
    }
    assert(found_hadronic && "Hadronic block for B_STANDARD should have been composed");

    manager->update(100.0, 2.0);
    assert(scale_setter->switches.size() == 4);
    assert(scale_setter->values.size()   == 4);

    std::cout << " CoefficientManager integration suite passed.\n";
    return 0;
}
