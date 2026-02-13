#include "ObservableExtractor.h"

std::vector<std::string> ObservableExtractor::schema(const OutputSpec& /*spec*/) const {
    // On veut un ordre stable, donc on trie par nom (string)
    // On DOIT connaître les bins à l’avance pour les inclure dans schema.
    // => Solution simple: on calcule une fois les observables pour découvrir les bins.
    // (c’est OK car schema() est appelé une fois au début du scan)

    std::vector<std::string> cols;

    auto obs_ids = oi_.get_current_observables();

    std::sort(obs_ids.begin(), obs_ids.end(),
              [](const ObservableId& a, const ObservableId& b) {
                  return ObservableMapper::str(a) < ObservableMapper::str(b);
              });

    for (const auto& oid : obs_ids) {
        const std::string name = ObservableMapper::str(oid);

        auto vals = oi_.compute_observable(oid); // renvoie vector<...> (binné ou non)
        if (vals.empty()) {
            // option: rien, ou colonne “NaN”
            continue;
        }

        for (const auto& v : vals) {
            cols.push_back(obs_key(name, ord_, v.bin));
        }
    }

    return cols;
}

void ObservableExtractor::extract(std::unordered_map<std::string, Value>& outY,
                                  const OutputSpec& /*spec*/) const
{
    auto obs_ids = oi_.get_current_observables();


    std::sort(obs_ids.begin(), obs_ids.end(),
              [](const ObservableId& a, const ObservableId& b) {
                  return ObservableMapper::str(a) < ObservableMapper::str(b);
              });

    for (const auto& oid : obs_ids) {
        const std::string name = ObservableMapper::str(oid);

        auto vals = oi_.compute_observable(oid);
        for (const auto& v : vals) {
            outY[obs_key(name, ord_, v.bin)] = v.value;
        }
    }
}
