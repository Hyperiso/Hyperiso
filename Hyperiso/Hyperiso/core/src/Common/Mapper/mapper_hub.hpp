#ifndef MAPPER_HUB_H
#define MAPPER_HUB_H

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

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
 * @file mapper_hub.hpp
 * @brief Non-templated façade over the project mapper families.
 *
 * The individual mappers are template-based and expose slightly different
 * registration APIs.  `mapper_hub.hpp` provides a small common interface for
 * generic tooling such as CLI completion, validation or configuration parsing.
 *
 * @note Generic registration is intentionally conservative.  Mappers requiring
 * extra semantic information, such as a parent decay for observables, expose
 * dedicated helper functions instead of using `IMapperView::register_custom()`.
 */

/**
 * @brief Enumeration of mapper families supported by the hub.
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
 * @brief Minimal runtime-polymorphic view over a mapper.
 */
struct IMapperView {
    virtual ~IMapperView() = default;

    /**
     * @brief Resolve a name or alias to its canonical string.
     *
     * @param s Input name or alias.
     * @return Canonical string stored by the mapper.
     */
    virtual std::string canonical(std::string_view s) const = 0;

    /**
     * @brief List all known canonical strings.
     *
     * @return Builtin and custom canonical strings known to the mapper.
     */
    virtual std::vector<std::string> list_all() const = 0;

    /**
     * @brief Generic custom registration hook.
     *
     * This is supported only for mappers whose custom symbols need no additional
     * semantic information.  For observables and decays, use the dedicated helper
     * functions below so invariants are preserved.
     *
     * @param canonical Canonical name.
     * @param aliases Aliases.
     * @return true on success, false if unsupported or rejected.
     */
    virtual bool register_custom(std::string canonical,
                                 std::vector<std::string> aliases) = 0;
};

/**
 * @brief Adapter for mappers with simple no-extra-data registration.
 *
 * The concrete mapper must expose:
 *   - `id_of(std::string_view)`,
 *   - `str(id)`,
 *   - `list_all()`,
 *   - `register_custom(std::string, std::vector<std::string>)`.
 */
template<class Concrete>
struct MapperViewSimple : IMapperView {
    std::string canonical(std::string_view s) const override {
        auto id = Concrete::id_of(s);
        return Concrete::str(id);
    }

    std::vector<std::string> list_all() const override {
        std::vector<std::string> out;
        for (auto& id : Concrete::list_all()) {
            out.push_back(Concrete::str(id));
        }
        return out;
    }

    bool register_custom(std::string canonical,
                         std::vector<std::string> aliases) override
    {
        return Concrete::register_custom(std::move(canonical), std::move(aliases));
    }
};

/**
 * @brief Adapter for mappers that can be queried generically but not registered.
 *
 * Use this for mapper families that need extra metadata during registration,
 * for example observables requiring a parent decay.
 */
template<class Concrete>
struct MapperViewReadOnly : IMapperView {
    std::string canonical(std::string_view s) const override {
        auto id = Concrete::id_of(s);
        return Concrete::str(id);
    }

    std::vector<std::string> list_all() const override {
        std::vector<std::string> out;
        for (auto& id : Concrete::list_all()) {
            out.push_back(Concrete::str(id));
        }
        return out;
    }

    bool register_custom(std::string, std::vector<std::string>) override {
        return false;
    }
};

/**
 * @brief Return the mapper view associated with a mapper kind.
 *
 * @param k Mapper family.
 * @return Reference to a static adapter instance.
 */
inline IMapperView& mapper_view(MapperKind k) {
    static MapperViewReadOnly<ObservableMapper>      obs;
    static MapperViewReadOnly<WCoefMapper>           wcoef;
    static MapperViewSimple<GroupMapper>             group;
    static MapperViewSimple<ParameterTypeMapper>     ptype;
    static MapperViewSimple<ModelMapper>             model;
    static MapperViewSimple<WilsonBasisMapper>       basis;
    static MapperViewSimple<ContributionTypeMapper>  contri;
    static MapperViewSimple<MassTypeMapper>          mass;
    static MapperViewSimple<ScaleTypeMapper>         scale;
    static MapperViewSimple<UncertaintyTypeMapper>   unc;
    static MapperViewReadOnly<DecayMapper>           decay;
    static MapperViewSimple<OrderMapper>             qcd;

    switch (k) {
        case MapperKind::Observable:      return obs;
        case MapperKind::WCoef:           return wcoef;
        case MapperKind::Group:           return group;
        case MapperKind::ParameterType:   return ptype;
        case MapperKind::Model:           return model;
        case MapperKind::WilsonBasis:     return basis;
        case MapperKind::ContributionType:return contri;
        case MapperKind::MassType:        return mass;
        case MapperKind::ScaleType:       return scale;
        case MapperKind::UncertaintyType: return unc;
        case MapperKind::Decay:           return decay;
        case MapperKind::QCDOrder:        return qcd;
    }

    return obs;
}

/**
 * @brief Resolve a name to its canonical representation for a mapper kind.
 */
inline std::string canonical_of(MapperKind k, std::string_view s) {
    return mapper_view(k).canonical(s);
}

/**
 * @brief List all canonical names known to a mapper kind.
 */
inline std::vector<std::string> list_all(MapperKind k) {
    return mapper_view(k).list_all();
}

/**
 * @brief Generic registration helper for simple mapper kinds.
 *
 * @return false for mapper kinds that require dedicated registration metadata.
 */
inline bool register_custom(MapperKind k,
                            std::string canonical,
                            std::vector<std::string> aliases)
{
    return mapper_view(k).register_custom(std::move(canonical), std::move(aliases));
}

/**
 * @brief Specialized registration helper for Wilson coefficients.
 *
 * @param canonical Canonical coefficient name.
 * @param aliases Optional aliases.
 * @param flha FLHA base pair.
 * @return true on success.
 */
inline bool register_custom_wcoef(std::string canonical,
                                  std::vector<std::string> aliases,
                                  std::pair<int,int> flha)
{
    return WCoefMapper::register_custom(
        std::move(canonical),
        std::move(aliases),
        flha
    );
}

/**
 * @brief Specialized registration helper for custom observables.
 *
 * A custom observable must always be attached to a parent decay.
 *
 * @param canonical Canonical observable name.
 * @param aliases Optional aliases.
 * @param ext Optional FLHA id.
 * @param parent_decay Parent decay name or alias.
 * @return true on success.
 */
inline bool register_custom_observable(std::string canonical,
                                       std::vector<std::string> aliases,
                                       std::optional<LhaID> ext,
                                       std::string_view parent_decay)
{
    return ObservableMapper::register_custom(
        canonical,
        std::move(aliases),
        std::move(ext),
        parent_decay
    );
}

/**
 * @brief Specialized registration helper for a custom decay with observables.
 *
 * This enforces the project invariant that a decay must have at least one
 * observable.
 *
 * @param canonical Canonical decay name.
 * @param aliases Optional decay aliases.
 * @param observables Custom observable definitions to attach.
 * @return true on success.
 */
inline bool register_custom_decay_with_observables(
    std::string canonical,
    std::vector<std::string> aliases,
    std::vector<CustomObservableSpec> observables)
{
    return DecayMapper::register_custom_with_observables(
        std::move(canonical),
        std::move(aliases),
        std::move(observables)
    );
}

#endif // MAPPER_HUB_H
