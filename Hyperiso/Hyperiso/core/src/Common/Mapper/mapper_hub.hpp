// mapper_hub.hpp
#pragma once
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

// Minimal Interface
struct IMapperView {
    virtual ~IMapperView() = default;
    virtual std::string canonical(std::string_view s) const = 0;
    virtual std::vector<std::string> list_all() const = 0;
    virtual bool register_custom(std::string canonical, std::vector<std::string> aliases) = 0;
};

// Concrete adapters (no ext)
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

// Concrete adapters (with ext)
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
        // For overload Hub::register_custom(kind, name, aliases, ext)
        // False by default
        return false;
    }
};

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
    return obs;
}

// helpers
inline std::string canonical_of(MapperKind k, std::string_view s){
    return mapper_view(k).canonical(s);
}
inline std::vector<std::string> list_all(MapperKind k){
    return mapper_view(k).list_all();
}
inline bool register_custom(MapperKind k, std::string canonical, std::vector<std::string> aliases){
    return mapper_view(k).register_custom(std::move(canonical), std::move(aliases));
}

// Overload with specific registering, for wcoef.
inline bool register_custom_wcoef(std::string canonical, std::vector<std::string> aliases, std::pair<int,int> flha){
    return WCoefMapper::register_custom(std::move(canonical), std::move(aliases), flha);
}
