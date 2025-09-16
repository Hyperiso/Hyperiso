// #pragma once
// #include "generic_mapper.hpp"
// #include "Map.h"

// struct DecayTag {};
// using DecayId = IdOf<DecayTag>;

// class DecayMapper
// : public GenericMapperNoExt<DecayTag, Decays, decays_mapping>
// {
// public:
//     using Base = GenericMapperNoExt<DecayTag, Decays, decays_mapping>;
//     using Base::init_builtins;
//     using Base::id_of;
//     using Base::str;
//     using Base::to_id;
//     using Base::enum_of;
//     using Base::list_all;

//     static std::vector<Observables> get_observables(Decays d) {
//         return decay_observable_mapping().at(d);
//     }

//     static Decays get_decay(Observables obs) {
//         for (const auto& [d, vec] : decay_observable_mapping()) {
//             if (std::find(vec.begin(), vec.end(), obs) != vec.end()) return d;
//         }
//         throw std::runtime_error("Observable not mapped to any decay");
//     }
// };

// decay_ids.hpp
#pragma once
#include "generic_mapper.hpp"
#include "Map.h"
#include "General.h"   // on choisit LhaID comme clé externe des decays (modifiable)
#include "ObservableMapper.h"

struct DecayTag {};
using DecayId = IdOf<DecayTag>;

// TODO : need builtin for decay
inline const std::map<Decays, LhaID>& decay_external_mapping(){
    static const std::map<Decays, LhaID> empty{};
    return empty;
}

class DecayMapper
: public GenericMapperWithExt<
      DecayTag,
      Decays,
      LhaID,
      std::hash<LhaID>,
      decays_mapping,
      decay_external_mapping
  >
{
public:
    using Base = GenericMapperWithExt<
        DecayTag, Decays, LhaID, std::hash<LhaID>,
        decays_mapping, decay_external_mapping
    >;

    using Base::init_builtins;
    using Base::id_of;
    using Base::str;         // str(DecayId), str(Decays)
    using Base::to_id;       // Decays -> DecayId
    using Base::enum_of;
    using Base::list_all;
    using Base::from_external;   // LhaID -> optional<DecayId>
    using Base::external_of;     // DecayId -> optional<LhaID>
    using Base::set_external; 

    static std::vector<Observables> get_observables(Decays d) {
        return decay_observable_mapping().at(d);
    }
    static Decays get_decay(Observables obs) {
        for (const auto& [d, vec] : decay_observable_mapping())
            if (std::find(vec.begin(), vec.end(), obs) != vec.end()) return d;
        throw std::runtime_error("Observable not mapped to any decay");
    }


    inline static const std::unordered_map<ObservableId, DecayId>& observable_to_decay() {
        static std::unordered_map<ObservableId, DecayId> cache;
        static std::once_flag init;
        std::call_once(init, []{
            for (const auto& [decayEnum, vec] : decay_observable_mapping()) {
                DecayId did = DecayMapper::to_id(decayEnum);
                for (auto o : vec)
                    cache.emplace(ObservableMapper::to_id(o), did);
            }
        });
        return cache;
    }

    inline static std::optional<DecayId> get_decay_id(ObservableId obs) {
        const auto& m = observable_to_decay();
        if (auto it = m.find(obs); it != m.end()) return it->second;
        return std::nullopt;
    }

};
