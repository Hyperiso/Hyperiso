#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>

#include "mapper_hub.hpp"
#include "registry_init.hpp"

int main() {
    std::cout << "== Dynamic mapper UNIT ==\n";

    init_all_builtins();

    GroupMapper::register_custom("UNIT_DYNAMIC_GROUP", {"unit-dyn-group"});
    WCoefMapper::register_custom("C_UNIT_DYNAMIC", {"c-unit-dyn"}, std::pair<int, int>{930001, 1});

    WGroupId group = GroupMapper::id_of("unit-dyn-group");
    WCoefId coeff = WCoefMapper::id_of("c-unit-dyn");

    assert(GroupMapper::str(group) == "UNIT_DYNAMIC_GROUP");
    assert(WCoefMapper::str(coeff) == "C_UNIT_DYNAMIC");
    assert(!GroupMapper::enum_of(group).has_value());
    assert(!WCoefMapper::enum_of(coeff).has_value());
    const auto expected_flha = std::make_pair(930001, 1);
    assert(WCoefMapper::flha_base(coeff) == expected_flha);

    CustomObservableSpec obs;
    obs.canonical = "UNIT_DYNAMIC_OBSERVABLE";
    obs.aliases = {"unit-dyn-observable"};
    obs.ext = LhaID(930002, 1);

    DecayMapper::register_custom_with_observables(
        "UNIT_DYNAMIC_DECAY",
        {"unit-dyn-decay"},
        {obs}
    );

    DecayId decay = DecayMapper::id_of("unit-dyn-decay");
    ObservableId observable = ObservableMapper::id_of("unit-dyn-observable");

    auto parent = DecayMapper::get_decay_id(observable);
    assert(parent.has_value());
    assert(parent.value() == decay);
    assert(!ObservableMapper::enum_of(observable).has_value());

    auto members = DecayMapper::get_observables(decay);
    assert(std::find(members.begin(), members.end(), observable) != members.end());

    std::cout << "UNIT OK\n";
    return 0;
}
