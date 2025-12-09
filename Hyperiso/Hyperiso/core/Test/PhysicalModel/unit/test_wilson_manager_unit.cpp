#include <cassert>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "WilsonManager.h"
#include "CustomWilson.h"
#include "CustomWilsonGroup.h"
#include "GroupMapper.h"
#include "wcoef_ids.hpp"
#include "SourcesView.h"
#include "Block.h"
#include "WilsonGroupAdapterConfig.h"
#include "IMartyWilsonProxy.h"
#include "InterpretedParam.h"

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
    std::unordered_set<std::string> removed_blocks;
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
        removed_blocks.insert(name);
        composed_blocks.erase(name);
    }

    void update(const std::string& /*name*/) override {
    }

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
    void calculate(std::string, std::string, double, std::string, bool = false) override {}
    std::set<std::string>  get_special_blocks() override { return {}; }
    std::unordered_set<InterpretedParam> get_dependencies(std::string) override { return {}; }
};

WilsonGroupAdapterConfig make_adapters(
    const std::shared_ptr<IParameterProxy<std::string, LhaID>>& wilson_proxy,
    const std::shared_ptr<IBlockComposer>& iblock_c)
{
    auto use_marty        = std::make_shared<DummyCoreAPI<bool>>(false);
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
    std::cout << "== CoefficientManager UNIT ==\n";

    auto iblock_c     = std::make_shared<RecordingBlockComposer>();
    auto wilson_proxy = std::make_shared<DummyParamProxy>();
    auto use_marty    = std::make_shared<DummyCoreAPI<bool>>(false);
    auto model_api    = std::make_shared<DummyCoreAPI<Model>>(Model::SM);
    auto scale_setter = std::make_shared<DummyScaleSetter>();

    PortsConfig ports{iblock_c, wilson_proxy, use_marty, model_api, scale_setter};
    ports.build_group = nullptr;

    CoefficientManager manager{ports};

    {
        assert(manager.getModel() == ModelMapper::str(Model::SM));
        model_api->value = Model::SUSY;
        assert(manager.getModel() == ModelMapper::str(Model::SUSY));
        model_api->value = Model::SM;
    }

    {
        manager.set_matching_scale(80.0);
        manager.set_hadronic_scale(5.0);

        assert(scale_setter->switches.size() == 2);
        assert(scale_setter->values.size()   == 2);

        assert(scale_setter->switches[0] == ScaleType::MATCHING);
        assert(scale_setter->values[0]   == 80.0);

        assert(scale_setter->switches[1] == ScaleType::HADRONIC);
        assert(scale_setter->values[1]   == 5.0);
    }

    {
        auto adapters = make_adapters(wilson_proxy, iblock_c);
        WGroupId gid  = GroupMapper::to_id(WGroup::BScalar);
        std::string group_name = GroupMapper::str(gid);

        auto grp = std::make_shared<CustomCoefficientGroup>(
            adapters,
            gid,
            "MY_CUSTOM_GROUP",
            ContributionType::SM
        );

        std::unordered_map<ParameterType, std::vector<std::string>> src_names{
            {ParameterType::SM, {"MASS"}},
        };
        grp->set_basis_order_sources_and_running(
            WilsonBasis::B_STANDARD,
            QCDOrder::LO,
            src_names,
            &CustomCoefficientGroup::identity_running
        );

        manager.registerCoefficientGroup(group_name, grp);

        auto got = manager.getCoefficientGroup(group_name);
        assert(got.get() == grp.get());

        auto all = manager.getGroups();
        assert(all.size() == 1);
        assert(all.contains(group_name));

        auto bases = manager.getGroupBases(gid);
        assert(bases.size() == 1);
        assert(bases.contains(WilsonBasis::B_STANDARD));

    }

    {
        auto adapters = make_adapters(wilson_proxy, iblock_c);
        WGroupId gid  = GroupMapper::to_id(WGroup::BScalar);
        std::string group_name = GroupMapper::str(gid);

        auto grp = std::make_shared<CustomCoefficientGroup>(
            adapters,
            gid,
            "MY_CUSTOM_GROUP2",
            ContributionType::SM
        );

        std::unordered_map<ParameterType, std::vector<std::string>> src_LO{
            {ParameterType::SM,     {"MASS"}},
            {ParameterType::WILSON, {"W_BLOCK_LO"}},
        };
        grp->set_basis_order_sources_and_running(
            WilsonBasis::B_STANDARD,
            QCDOrder::LO,
            src_LO,
            &CustomCoefficientGroup::identity_running
        );

        std::unordered_map<ParameterType, std::vector<std::string>> src_NLO{
            {ParameterType::SM,     {"GAUGE"}},
            {ParameterType::WILSON, {"W_BLOCK_NLO"}},
        };
        grp->set_basis_order_sources_and_running(
            WilsonBasis::B_STANDARD,
            QCDOrder::NLO,
            src_NLO,
            &CustomCoefficientGroup::identity_running
        );

        CoefficientManager mgr2{ports};
        mgr2.registerCoefficientGroup(group_name, grp);

        std::unordered_map<ParameterType, std::vector<std::string>> collected;

        mgr2.fill_sources_for_group(group_name, "NLO", collected, WilsonBasis::B_STANDARD);

        assert(collected.size() == 2);
        auto itSM = collected.find(ParameterType::SM);
        auto itW  = collected.find(ParameterType::WILSON);
        assert(itSM != collected.end());
        assert(itW  != collected.end());

        const auto& vSM = itSM->second;
        bool has_mass  = std::find(vSM.begin(), vSM.end(), "MASS")  != vSM.end();
        bool has_gauge = std::find(vSM.begin(), vSM.end(), "GAUGE") != vSM.end();
        assert(has_mass && has_gauge);

        const auto& vW = itW->second;
        bool has_lo  = std::find(vW.begin(), vW.end(), "W_BLOCK_LO")  != vW.end();
        bool has_nlo = std::find(vW.begin(), vW.end(), "W_BLOCK_NLO") != vW.end();
        assert(has_lo && has_nlo);
    }

    std::cout << " CoefficientManager unit suite passed.\n";
    return 0;
}
