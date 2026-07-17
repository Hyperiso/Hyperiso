#include <iostream>
#include <set>
#include <stdexcept>
#include <string_view>

#include "BinnedObservableId.h"
#include "GeneralEnum.h"
#include "LhaID.h"
#include "Map.h"
#include "Mapper/observable_ids.hpp"
#include "ObservableTypeId.h"

namespace {

void require(bool condition, std::string_view message) {
    if (!condition) {
        throw std::runtime_error(std::string(message));
    }
}

template <typename Callable>
void require_invalid_argument(Callable&& callable, std::string_view message) {
    try {
        callable();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(std::string(message));
}

} // namespace

int main() {
    using enum Observables;

    static_assert(make_polarization_type_id(PolarizationKind::Fermion, 2) == 9212);
    static_assert(make_polarization_type_id(PolarizationKind::LongitudinalVector, 1) == 9221);
    static_assert(make_polarization_type_id(PolarizationKind::TransverseVector, 1) == 9231);

    require_invalid_argument(
        [] { (void)make_polarization_type_id(PolarizationKind::Fermion, 0); },
        "daughter position zero must be rejected"
    );
    require_invalid_argument(
        [] { (void)make_polarization_type_id(PolarizationKind::Fermion, 10); },
        "daughter positions above nine must be rejected"
    );
    require_invalid_argument(
        [] {
            (void)make_polarization_type_id(
                static_cast<PolarizationKind>(0),
                1
            );
        },
        "unknown polarization kinds must be rejected"
    );

    const auto& mapping = observable_flha_mapping();

    require(
        mapping.at(P_TAU_B0__D_TAU_NU) == LhaID(511, 9212, 3, 411, -15, 16),
        "B0 -> D tau nu must use fermion polarization 9212"
    );
    require(
        mapping.at(P_TAU_B__D0_TAU_NU) == LhaID(521, 9212, 3, 421, -15, 16),
        "B -> D0 tau nu must use fermion polarization 9212"
    );
    require(
        mapping.at(P_D_B0__DSTAR_TAU_NU) == LhaID(511, 9221, 3, 413, -15, 16),
        "B0 -> D* tau nu must use longitudinal-vector polarization 9221"
    );
    require(
        mapping.at(P_D_B__DSTAR0_TAU_NU) == LhaID(521, 9221, 3, 423, -15, 16),
        "B -> D*0 tau nu must use longitudinal-vector polarization 9221"
    );

    require(
        mapping.at(F_T_B0__KSTAR0_MU_MU) ==
            LhaID(511, kTransverseFractionTypeId, 3, 313, 13, -13),
        "F_T must use type 932"
    );
    require(
        mapping.at(ALPHA_K_B0__KSTAR0_MU_MU) ==
            LhaID(511, kAlphaKTypeId, 3, 313, 13, -13),
        "alpha_K must use type 933"
    );
    require(
        mapping.at(F_T_B0__KSTAR0_MU_MU) !=
            mapping.at(ALPHA_K_B0__KSTAR0_MU_MU),
        "F_T and alpha_K must have distinct FLHA ids"
    );

    const auto type_932 = ObservableMapper::from_flha(
        LhaID(511, kTransverseFractionTypeId, 3, 313, 13, -13)
    );
    require(type_932.has_value(), "type 932 must resolve");
    require(
        type_932.value() == ObservableMapper::to_id(F_T_B0__KSTAR0_MU_MU),
        "type 932 must remain the canonical F_T identifier"
    );

    const auto legacy_tau = ObservableMapper::from_flha(
        LhaID(511, kLegacyFermionPolarizationTypeId, 3, 411, -15, 16)
    );
    require(legacy_tau.has_value(), "legacy fermion polarization must resolve");
    require(
        legacy_tau.value() == ObservableMapper::to_id(P_TAU_B0__D_TAU_NU),
        "legacy fermion polarization must resolve to P_TAU"
    );

    const auto legacy_dstar = ObservableMapper::from_flha(
        LhaID(511, kLegacyLongitudinalVectorPolarizationTypeId, 3, 413, -15, 16)
    );
    require(
        legacy_dstar.has_value(),
        "legacy longitudinal-vector polarization must resolve"
    );
    require(
        legacy_dstar.value() == ObservableMapper::to_id(P_D_B0__DSTAR_TAU_NU),
        "legacy longitudinal-vector polarization must resolve to P_D"
    );

    require(
        ObservableMapper::enum_elt(
            LhaID(511, kLegacyFermionPolarizationTypeId, 3, 411, -15, 16)
        ) == P_TAU_B0__D_TAU_NU,
        "enum_elt must accept the legacy fermion polarization id"
    );

    const BinnedObservableId expected_binned(
        P_TAU_B0__D_TAU_NU,
        {1.0, 2.0}
    );
    auto legacy_binned_parts = expected_binned.flha().get_parts();
    legacy_binned_parts[1] = kLegacyFermionPolarizationTypeId;
    const auto legacy_binned = ObservableMapper::from_binned_flha(
        LhaID(std::move(legacy_binned_parts))
    );
    require(legacy_binned.has_value(), "legacy binned polarization must resolve");
    require(
        legacy_binned.value() == expected_binned,
        "legacy binned polarization must preserve the observable and bin"
    );

    std::set<LhaID> unique_ids;
    for (const auto& [observable, flha] : mapping) {
        (void)observable;
        require(
            unique_ids.insert(flha).second,
            "builtin observables must have unique FLHA ids"
        );
    }

    std::cout << "Observable type-id tests passed.\n";
    return 0;
}
