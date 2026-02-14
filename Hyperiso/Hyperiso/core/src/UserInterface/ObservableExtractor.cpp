#include "ObservableExtractor.h"

static std::string obs_name(const ObservableId& obs) {
    // adapte si ton type a un autre accès au nom
    // dans ton code tu fais elem.str()
    return obs.str();
}

std::string obs_key(const ObservableId& obs, QCDOrder ord) {
    std::string k = "OBS:";
    k += obs_name(obs);
    k += ":";
    k += OrderMapper::str(ord);
    return k;
}

std::string obs_key_bin(const ObservableId& obs, QCDOrder ord, double bmin, double bmax) {
    // clé stable + lisible
    // (tu peux formater avec plus/moins de digits si tu veux)
    std::string k = obs_key(obs, ord);
    k += ":BIN[";
    k += std::to_string(bmin);
    k += ",";
    k += std::to_string(bmax);
    k += "]";
    return k;
}

std::vector<std::string> ObservableExtractor::schema(const OutputSpec& /*spec*/) const {
    std::vector<std::string> cols;

    // IMPORTANT: oi_.get_current_observables() renvoie très probablement un unordered_set
    // => on copie dans un vector pour trier
    auto current = oi_.get_current_observables();
    oi_.enable_obs();
    std::vector<ObservableId> sorted(current.begin(), current.end());
    std::sort(sorted.begin(), sorted.end(),
              [](const ObservableId& a, const ObservableId& b) {
                  return a.str() < b.str();
              });

    // Ici il y a 2 choix:
    // 1) Schema "minimal" : 1 colonne par observable, et si binned -> colonnes créées au moment du 1er extract
    //    (mais CSV header fixe => pas top)
    // 2) Schema "stable" : on appelle compute_observable 1 fois par observable pour connaître les bins
    //    (c’est ce que je fais, comme ça ton header CSV est complet)
    for (const auto& obs : sorted) {
        auto vals = oi_.compute_observable(obs); // renvoie vector de valeurs (avec bin optionnel)
        if (vals.empty()) {
            cols.push_back(obs_key(obs, ord_));
            continue;
        }

        // Si la première valeur n'a pas de bin => on suppose non-binned (1 colonne)
        if (!vals[0].bin.has_value()) {
            cols.push_back(obs_key(obs, ord_));
            continue;
        }

        // Sinon binned => 1 colonne par bin
        for (const auto& v : vals) {
            if (!v.bin.has_value()) {
                // mélange bizarre: fallback sur non-binned
                cols.push_back(obs_key(obs, ord_));
                break;
            }
            const auto& br = v.bin.value();
            cols.push_back(obs_key_bin(obs, ord_, br.first, br.second));
        }
    }

    return cols;
}

void ObservableExtractor::extract(std::unordered_map<std::string, Value>& outY,
                                 const OutputSpec& /*spec*/) const
{
    auto current = oi_.get_current_observables();
    oi_.enable_obs();
    std::vector<ObservableId> sorted(current.begin(), current.end());
    std::sort(sorted.begin(), sorted.end(),
              [](const ObservableId& a, const ObservableId& b) {
                  return a.str() < b.str();
              });

    for (const auto& obs : sorted) {
        auto vals = oi_.compute_observable(obs);

        if (vals.empty()) {
            // si vide, on écrit NaN ou 0 selon ta convention
            outY[obs_key(obs, ord_)] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        if (!vals[0].bin.has_value()) {
            // non-binned: on prend vals[0]
            outY[obs_key(obs, ord_)] = vals[0].value;
            continue;
        }

        // binned
        for (const auto& v : vals) {
            if (!v.bin.has_value()) {
                // fallback
                outY[obs_key(obs, ord_)] = v.value;
                break;
            }
            const auto& br = v.bin.value();
            outY[obs_key_bin(obs, ord_, br.first, br.second)] = v.value;
        }
    }
}