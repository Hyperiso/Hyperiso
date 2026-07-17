#include "observable_ids.hpp"

#include "BinnedObservableId.h"
#include "decay_ids.hpp"
#include "ObservableTypeId.h"

#include <optional>
#include <stdexcept>
#include <utility>

namespace {

std::optional<LhaID> canonicalize_legacy_polarization_id(const LhaID& ext) {
    auto parts = ext.get_parts();
    if (parts.size() < 2) {
        return std::nullopt;
    }

    if (parts[1] == kLegacyFermionPolarizationTypeId) {
        parts[1] = make_polarization_type_id(PolarizationKind::Fermion, 2);
    } else if (parts[1] == kLegacyLongitudinalVectorPolarizationTypeId) {
        parts[1] = make_polarization_type_id(
            PolarizationKind::LongitudinalVector,
            1
        );
    } else {
        return std::nullopt;
    }

    return LhaID(std::move(parts));
}

} // namespace

Observables ObservableMapper::enum_elt_legacy(std::string_view s) {
    const auto key = nk(s);

    for (const auto& [e, name] : observable_mapping()) {
        if (nk(name) == key) {
            return e;
        }
    }

    throw std::out_of_range("Unknown builtin observable: " + std::string(s));
}

Observables ObservableMapper::enum_elt(const LhaID& ext) {
    auto maybeId = from_flha(ext);
    if (!maybeId) {
        throw std::out_of_range("Unknown observable (FLHA): " + ext.to_string());
    }

    auto maybeEnum = enum_of(*maybeId);
    if (!maybeEnum) {
        throw std::out_of_range("ObservableId has no builtin enum: " + maybeId->str());
    }

    return *maybeEnum;
}

std::optional<ObservableId> ObservableMapper::from_flha(const LhaID& ext) {
    if (auto direct = from_external(ext)) {
        return direct;
    }

    const auto canonical = canonicalize_legacy_polarization_id(ext);
    if (!canonical) {
        return std::nullopt;
    }
    return from_external(*canonical);
}

std::optional<LhaID> ObservableMapper::flha_of(const ObservableId& id) {
    return external_of(id);
}

std::optional<BinnedObservableId> ObservableMapper::from_binned_flha(const LhaID& ext) {
    try {
        return BinnedObservableId::from_flha(ext);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

std::optional<LhaID> ObservableMapper::binned_flha_of(const BinnedObservableId& id) {
    try {
        return id.flha();
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

BinnedObservableId ObservableMapper::binned_from_flha(const LhaID& ext) {
    auto out = from_binned_flha(ext);
    if (!out) {
        throw std::runtime_error("Unknown binned observable (FLHA): " + ext.to_string());
    }
    return *out;
}

LhaID ObservableMapper::flha(const Observables& id) {
    return observable_flha_mapping().at(id);
}

bool ObservableMapper::register_custom(const std::string& canonical,
                                       std::vector<std::string> aliases,
                                       std::optional<LhaID> ext,
                                       const DecayId& parent_decay)
{
    // Resolve through DecayMapper so aliases are canonicalized and unknown
    // parent decays are rejected before observable registration.
    const DecayId canonical_decay = DecayMapper::id_of(parent_decay.str());

    const bool ok = Base::register_custom(
        canonical,
        std::move(aliases),
        std::move(ext)
    );

    if (!ok) {
        return false;
    }

    const ObservableId obs_id = Base::id_of(canonical);
    DecayGraph::instance().link(canonical_decay, obs_id);
    return true;
}

bool ObservableMapper::register_custom(const std::string& canonical,
                                       std::vector<std::string> aliases,
                                       std::optional<LhaID> ext,
                                       Decays parent_decay_enum)
{
    return register_custom(
        canonical,
        std::move(aliases),
        std::move(ext),
        DecayMapper::to_id(parent_decay_enum)
    );
}

bool ObservableMapper::register_custom(const std::string& canonical,
                                       std::vector<std::string> aliases,
                                       std::optional<LhaID> ext,
                                       std::string_view parent_decay_str)
{
    return register_custom(
        canonical,
        std::move(aliases),
        std::move(ext),
        DecayMapper::id_of(parent_decay_str)
    );
}

LhaID ObservableMapper::flha(const ObservableId& id) {
    auto k = flha_of(id);
    if (!k) {
        throw std::runtime_error("Observable has no FLHA: " + id.str());
    }
    return *k;
}

LhaID ObservableMapper::flha(const BinnedObservableId& id) {
    auto k = binned_flha_of(id);
    if (!k) {
        throw std::runtime_error("BinnedObservableId has no FLHA: " + id.str());
    }
    return *k;
}
