#pragma once
#include "generic_mapper.hpp"
#include "Map.h"
#include "General.h"
#include "decay_graph.hpp"

struct ObservableTag {};
using ObservableId = IdOf<ObservableTag>;

class ObservableMapper
: public GenericMapperWithExt<
      ObservableTag,
      Observables,
      LhaID,
      std::hash<LhaID>,
      observable_mapping,
      observable_flha_mapping
  >
{
public:
    using Base = GenericMapperWithExt<
        ObservableTag, Observables, LhaID, std::hash<LhaID>,
        observable_mapping, observable_flha_mapping
    >;

    using Base::init_builtins;
    using Base::id_of;           // string -> ObservableId
    using Base::str;             // str(ObservableId), str(Observables)
    using Base::to_id;           // Observables -> ObservableId
    using Base::enum_of;         // ObservableId -> optional<Observables>
    using Base::list_all;
    using Base::from_external;   // LhaID -> optional<ObservableId>
    using Base::external_of;     // ObservableId -> optional<LhaID>
    using Base::set_external;

    static Observables enum_elt(std::string_view s){
        const auto key = nk(s);
        for (const auto& [e, name] : observable_mapping()){
            if (nk(name) == key) return e;
        }
        throw std::out_of_range("Unknown observable: " + std::string(s));
    }

    static Observables enum_elt(const LhaID& ext){
        auto maybeId = from_external(ext);
        if (!maybeId) throw std::out_of_range("Unknown observable (FLHA): " + ext.to_string());
        auto maybeEnum = enum_of(*maybeId);
        if (!maybeEnum) throw std::out_of_range("ObservableId has no builtin enum: " + maybeId->str());
        return *maybeEnum;
    }

    static std::optional<ObservableId> from_flha(const LhaID& ext){ return from_external(ext); }
    static std::optional<LhaID>        flha_of(const ObservableId& id){ return external_of(id); }
    static LhaID        flha(const Observables& id){ return observable_mapping().at(id); }
    static bool register_custom(const std::string& canonical,
                                std::vector<std::string> aliases = {},
                                std::optional<LhaID> ext = std::nullopt,
                                std::optional<DecayId> parent_decay = std::nullopt)
    {
        const bool ok = Base::register_custom(canonical, std::move(aliases), std::move(ext));
        if (ok && parent_decay){
            // link to decay graph
            DecayGraph::instance().link(*parent_decay, ObservableId(canonical));
        }
        return ok;
    }

    static LhaID flha(const ObservableId& id){
        auto k = flha_of(id);
        if (!k) throw std::runtime_error("Observable has no FLHA: " + id.str());
        return *k;
    }

    // static LhaID flha(Observables e){
    //     return observable_flha_mapping().at(e);
    // }

    static bool register_custom(const std::string& canonical,
                                std::vector<std::string> aliases,
                                std::optional<LhaID> ext,
                                Decays parent_decay_enum);

    static bool register_custom(const std::string& canonical,
                                std::vector<std::string> aliases,
                                std::optional<LhaID> ext,
                                std::string_view parent_decay_str);
};
