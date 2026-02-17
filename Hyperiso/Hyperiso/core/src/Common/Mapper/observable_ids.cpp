#include "observable_ids.hpp"
#include "decay_ids.hpp"

Observables ObservableMapper::enum_elt(std::string_view s){
    const auto key = nk(s);
    for (const auto& [e, name] : observable_mapping()){
        if (nk(name) == key) return e;
    }
    throw std::out_of_range("Unknown observable: " + std::string(s));
}

Observables ObservableMapper::enum_elt(const LhaID& ext){
    auto maybeId = from_external(ext);
    if (!maybeId) throw std::out_of_range("Unknown observable (FLHA): " + ext.to_string());
    auto maybeEnum = enum_of(*maybeId);
    if (!maybeEnum) throw std::out_of_range("ObservableId has no builtin enum: " + maybeId->str());
    return *maybeEnum;
}

std::optional<ObservableId> ObservableMapper::from_flha(const LhaID& ext){ return from_external(ext); }

std::optional<LhaID>        ObservableMapper::flha_of(const ObservableId& id){ return external_of(id); }

LhaID        ObservableMapper::flha(const Observables& id){ return observable_flha_mapping().at(id); }

bool ObservableMapper::register_custom(const std::string& canonical,
                        std::vector<std::string> aliases,
                        std::optional<LhaID> ext,
                        std::optional<DecayId> parent_decay)
{
    const bool ok = Base::register_custom(canonical, aliases, ext);
    if (ok && parent_decay) {
        auto obs_id = Base::id_of(canonical);
        DecayGraph::instance().link(*parent_decay, obs_id);
    }
    return ok;
}

LhaID ObservableMapper::flha(const ObservableId& id){
    auto k = flha_of(id);
    if (!k) throw std::runtime_error("Observable has no FLHA: " + id.str());
    return *k;
}

bool ObservableMapper::register_custom(const std::string& canonical,
                            std::vector<std::string> aliases,
                            std::optional<LhaID> ext,
                            Decays parent_decay_enum)
{
    return register_custom(canonical, std::move(aliases), std::move(ext),
                            DecayMapper::to_id(parent_decay_enum));
}

bool ObservableMapper::register_custom(const std::string& canonical,
                            std::vector<std::string> aliases,
                            std::optional<LhaID> ext,
                            std::string_view parent_decay_str)
{
    return register_custom(canonical, std::move(aliases), std::move(ext),
                            DecayMapper::id_of(parent_decay_str));
}