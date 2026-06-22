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
    std::unordered_set<WGroupId> last_groups;

    void build(std::shared_ptr<AbstractConfig> cfg) override {
        build_calls++;
        last_groups.clear();

        auto wil_cfg = std::dynamic_pointer_cast<WilsonBuildConfig>(cfg);
        if (wil_cfg) last_groups = wil_cfg->groups;
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
    std::cout << "== ObsWilsonHelper UNIT ==\n";

    ObsWilsonHelper helper;

    auto builder_spy = std::make_shared<SpyObsWilsonBuilder>();
    std::shared_ptr<IObsWilsonBuilder> builder = builder_spy;

    auto freezer_spy = std::make_shared<SpyWilsonFreezer>();
    std::shared_ptr<IWilsonFreezer<WGroupId>> freezer = freezer_spy;

    const WGroupId g1 = GroupMapper::to_id(WGroup::B);
    const WGroupId g2 = GroupMapper::to_id(WGroup::BPrime);
    const WGroupId g3 = GroupMapper::to_id(WGroup::BScalar);

    {
        freezer_spy->clear();
        auto cfg = make_cfg({g1, g2});

        helper.build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 1);
        assert(builder_spy->last_groups == set_from({g1, g2}));
        assert(freezer_spy->froze.empty());
        assert(freezer_spy->unfroze.empty());
    }

    {
        freezer_spy->clear();
        auto cfg = make_cfg({g1, g2});

        helper.build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 1);
        assert(freezer_spy->froze.empty());
        assert(freezer_spy->unfroze.empty());
    }

    {
        freezer_spy->clear();
        auto cfg = make_cfg({g1});

        helper.build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 1);
        assert(set_from(freezer_spy->froze) == set_from({g2}));
        assert(freezer_spy->unfroze.empty());
    }

    {
        freezer_spy->clear();
        auto cfg = make_cfg({g1, g2});

        helper.build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 1);
        assert(freezer_spy->froze.empty());
        assert(set_from(freezer_spy->unfroze) == set_from({g2}));
    }

    {
        freezer_spy->clear();
        auto cfg = make_cfg({g3});

        helper.build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 2);
        assert(builder_spy->last_groups == set_from({g3}));
        assert(set_from(freezer_spy->froze) == set_from({g1, g2}));
        assert(freezer_spy->unfroze.empty());
    }

    {
        helper.clear();
        freezer_spy->clear();
        auto cfg = make_cfg({g1});

        helper.build(cfg, builder, freezer);

        assert(builder_spy->build_calls == 3);
        assert(builder_spy->last_groups == set_from({g1}));
        assert(freezer_spy->froze.empty());
        assert(freezer_spy->unfroze.empty());
    }

    std::cout << "UNIT OK\n";
    return 0;
}
