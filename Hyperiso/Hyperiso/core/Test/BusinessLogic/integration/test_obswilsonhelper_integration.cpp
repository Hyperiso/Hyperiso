#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <vector>

#include "Include.h"

#ifdef LOG_DEBUG
#undef LOG_DEBUG
#endif
#define LOG_DEBUG(...) do {} while (0)

#include "ObsWilsonHelper.h"
#include "IObsWilsonBuilder.h"
#include "IWilsonFreezer.h"
#include "Configs.h"

class SpyObsWilsonBuilder : public IObsWilsonBuilder {
public:
    int build_calls = 0;
    std::vector<std::unordered_set<WGroupId>> history_groups;

    void build(std::shared_ptr<AbstractConfig> cfg) override {
        build_calls++;
        auto wil_cfg = std::dynamic_pointer_cast<WilsonBuildConfig>(cfg);
        if (wil_cfg) history_groups.push_back(wil_cfg->groups);
        else         history_groups.push_back({});
    }

    void add_custom_group(const CustomWilsonGroupConfig&) override {}
    std::shared_ptr<IObsWilsonProxy> get_proxy() override { return nullptr; }
};

class SpyWilsonFreezer : public IWilsonFreezer<WGroupId> {
public:
    std::vector<WGroupId> froze;
    std::vector<WGroupId> unfroze;

    void freeze(WGroupId g) override   { froze.push_back(g); }
    void unfreeze(WGroupId g) override { unfroze.push_back(g); }

    void clear() { froze.clear(); unfroze.clear(); }
};

static std::unordered_set<WGroupId> set_from(std::initializer_list<WGroupId> xs) {
    return std::unordered_set<WGroupId>(xs.begin(), xs.end());
}
static std::unordered_set<WGroupId> set_from(const std::vector<WGroupId>& xs) {
    return std::unordered_set<WGroupId>(xs.begin(), xs.end());
}
static WilsonBuildConfig make_cfg(std::initializer_list<WGroupId> groups) {
    WilsonBuildConfig cfg{};
    cfg.groups = set_from(groups);
    return cfg;
}

int main() {
    std::cout << "== ObsWilsonHelper INTEGRATION (ports) ==\n";

    ObsWilsonHelper(true);

    auto builder_spy = std::make_shared<SpyObsWilsonBuilder>();
    std::shared_ptr<IObsWilsonBuilder> builder = builder_spy;
    auto freezer = std::make_shared<SpyWilsonFreezer>();

    const WGroupId A = GroupMapper::to_id(WGroup::B);
    const WGroupId B = GroupMapper::to_id(WGroup::BPrime);
    const WGroupId C = GroupMapper::to_id(WGroup::BScalar);

    // {A,B} => build({A,B})
    {
        freezer->clear();
        auto cfg = make_cfg({A, B});
        ObsWilsonHelper::build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 1);
        assert(builder_spy->history_groups.back() == set_from({A, B}));
        assert(freezer->froze.empty());
        assert(freezer->unfroze.empty());
    }

    // {B} => A not needed => freeze(A)
    {
        freezer->clear();
        auto cfg = make_cfg({B});
        ObsWilsonHelper::build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 1); 
        assert(set_from(freezer->froze) == set_from({A}));
        assert(freezer->unfroze.empty());
    }

    // {A,B} => A était frozen => unfreeze(A)
    {
        freezer->clear();
        auto cfg = make_cfg({A, B});
        ObsWilsonHelper::build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 1);
        assert(freezer->froze.empty());
        assert(set_from(freezer->unfroze) == set_from({A}));
    }

    // {} => freeze (A,B)
    {
        freezer->clear();
        auto cfg = make_cfg({});
        ObsWilsonHelper::build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 1); 
        assert(set_from(freezer->froze) == set_from({A, B}));
        assert(freezer->unfroze.empty());
    }

    // {C,B} => B frozen => unfreeze(B), C  => build({C})
    //          A => freeze(A) (A est connu)
    {
        freezer->clear();
        auto cfg = make_cfg({C, B});
        ObsWilsonHelper::build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 2);
        assert(builder_spy->history_groups.back() == set_from({C}));

        assert(set_from(freezer->unfroze) == set_from({B}));
        assert(set_from(freezer->froze)   == set_from({A}));
    }

    std::cout << "INTEGRATION OK\n";
    return 0;
}
