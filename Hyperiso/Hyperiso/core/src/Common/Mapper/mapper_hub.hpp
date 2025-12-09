#ifndef MAPPER_HUB_H
#define MAPPER_HUB_H

#include <memory>
#include <vector>
#include <string>
#include <string_view>
#include "observable_ids.hpp"
#include "wcoef_ids.hpp"
#include "wgroup_ids.hpp"
#include "parametertype_ids.hpp"
#include "model_ids.hpp"
#include "wilsonbasis_ids.hpp"
#include "contributiontype_ids.hpp"
#include "masstype_ids.hpp"
#include "scaletype_ids.hpp"
#include "uncertaintytype_ids.hpp"
#include "decay_ids.hpp"
#include "qcdorder_ids.hpp"

/**
 * @file mapper_hub.h
 * @brief Unified view and minimal interface for all mapper types.
 *
 * This header defines:
 *   - MapperKind: an enum selecting which mapper family to use,
 *   - IMapperView: a minimal common interface (canonicalization, listing, registration),
 *   - MapperViewNoExt / MapperViewWithExt: adapters from concrete mappers to IMapperView,
 *   - mapper_view(): factory returning a reference to the appropriate adapter,
 *   - helper functions (canonical_of, list_all, register_custom, register_custom_wcoef).
 *
 * It is intended for code that needs to interact generically with
 * different mapper types without templating on each of them.
 */

/**
 * @brief Enumeration of all supported mapper kinds.
 *
 * Each entry corresponds to a concrete mapper type, e.g. ObservableMapper,
 * WCoefMapper, GroupMapper, etc.
 */
enum class MapperKind {
    Observable,
    WCoef,
    Group,
    ParameterType,
    Model,
    WilsonBasis,
    ContributionType,
    MassType,
    ScaleType,
    UncertaintyType,
    Decay,
    QCDOrder
};

/**
 * @brief Minimal non-templated interface for mapper operations.
 *
 * This interface exposes only:
 *   - canonical: normalize/resolve a string to its canonical name,
 *   - list_all: list all known canonical names,
 *   - register_custom: register a user-defined symbol with aliases.
 */
struct IMapperView {
    virtual ~IMapperView() = default;

    /// Resolves a name (case-insensitive) to the canonical string.
    virtual std::string canonical(std::string_view s) const = 0;

    /// Returns all canonical names known to this mapper.
    virtual std::vector<std::string> list_all() const = 0;

    /// Registers a custom symbol with aliases (no external key).
    virtual bool register_custom(std::string canonical, std::vector<std::string> aliases) = 0;
};

/**
 * @brief Adapter for mappers without external keys.
 *
 * The Concrete type is expected to expose:
 *   - enum_elt(std::string_view)
 *   - str(id)
 *   - list_all()
 *   - register_custom(std::string, std::vector<std::string>)
 */
template<class Concrete> struct MapperViewNoExt : IMapperView {
    std::string canonical(std::string_view s) const override {
        auto id = Concrete::enum_elt(s);
        return Concrete::str(id);
    }
    std::vector<std::string> list_all() const override {
        std::vector<std::string> out;
        for (auto& id : Concrete::list_all()) out.push_back(Concrete::str(id));
        return out;
    }
    bool register_custom(std::string canonical, std::vector<std::string> aliases) override {
        return Concrete::register_custom(std::move(canonical), std::move(aliases));
    }
};

/**
 * @brief Adapter for mappers with external keys.
 *
 * The Concrete type is expected to expose:
 *   - enum_elt(std::string_view)
 *   - str(id)
 *   - list_all()
 *
 * By default, register_custom returns false here; a dedicated overload
 * (e.g. register_custom_wcoef) is provided for specific types that
 * require an external key.
 */
template<class Concrete> struct MapperViewWithExt : IMapperView {
    std::string canonical(std::string_view s) const override {
        auto id = Concrete::enum_elt(s);
        return Concrete::str(id);
    }
    std::vector<std::string> list_all() const override {
        std::vector<std::string> out;
        for (auto& id : Concrete::list_all()) out.push_back(Concrete::str(id));
        return out;
    }
    bool register_custom(std::string canonical, std::vector<std::string> aliases) override {
        // Registration requiring an external key is not supported through this interface.
        return false;
    }
};

/**
 * @brief Returns a reference to the IMapperView corresponding to the given kind.
 *
 * The returned object is a static adapter (no ownership transfer) and
 * must not be deleted by the caller.
 */
inline IMapperView& mapper_view(MapperKind k){
    static MapperViewNoExt<ObservableMapper>      obs;
    static MapperViewWithExt<WCoefMapper>         wcoef;
    static MapperViewNoExt<GroupMapper>           group;
    static MapperViewNoExt<ParameterTypeMapper>   ptype;
    static MapperViewNoExt<ModelMapper>           model;
    static MapperViewNoExt<WilsonBasisMapper>     basis;
    static MapperViewNoExt<ContributionTypeMapper>contri;
    static MapperViewNoExt<MassTypeMapper>        mass;
    static MapperViewNoExt<ScaleTypeMapper>       scale;
    static MapperViewNoExt<UncertaintyTypeMapper> unc;
    static MapperViewNoExt<DecayMapper>           decay;
    static MapperViewNoExt<OrderMapper>           qcd;

    switch (k){
        case MapperKind::Observable:     return obs;
        case MapperKind::WCoef:          return wcoef;
        case MapperKind::Group:          return group;
        case MapperKind::ParameterType:  return ptype;
        case MapperKind::Model:          return model;
        case MapperKind::WilsonBasis:    return basis;
        case MapperKind::ContributionType:return contri;
        case MapperKind::MassType:       return mass;
        case MapperKind::ScaleType:      return scale;
        case MapperKind::UncertaintyType:return unc;
        case MapperKind::Decay:          return decay;
        case MapperKind::QCDOrder:       return qcd;
    }
    // Fallback: should not happen, default to observable mapper.
    return obs;
}

// ------------------------------------------------------------------
// Helper free functions that delegate to the selected mapper.
// ------------------------------------------------------------------

/**
 * @brief Resolves a name to its canonical representation for a given mapper kind.
 */
inline std::string canonical_of(MapperKind k, std::string_view s){
    return mapper_view(k).canonical(s);
}

/**
 * @brief Lists all canonical names for a given mapper kind.
 */
inline std::vector<std::string> list_all(MapperKind k){
    return mapper_view(k).list_all();
}

/**
 * @brief Registers a custom symbol (no external key) for a given mapper kind.
 */
inline bool register_custom(MapperKind k, std::string canonical, std::vector<std::string> aliases){
    return mapper_view(k).register_custom(std::move(canonical), std::move(aliases));
}

/**
 * @brief Specialized registration helper for WCoef with a FLHA base key.
 *
 * This bypasses the generic IMapperView interface to allow passing the
 * external FLHA pair (a,b) required by WCoefMapper.
 *
 * @param canonical Canonical name of the Wilson coefficient.
 * @param aliases   Alias names.
 * @param flha      FLHA base indices (a,b).
 * @return true if registration succeeded, false otherwise.
 */
inline bool register_custom_wcoef(std::string canonical, std::vector<std::string> aliases, std::pair<int,int> flha){
    return WCoefMapper::register_custom(std::move(canonical), std::move(aliases), flha);
}

#endif