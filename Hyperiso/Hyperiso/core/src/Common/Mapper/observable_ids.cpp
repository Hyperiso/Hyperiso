#include "observable_ids.hpp"
#include "decay_ids.hpp"
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