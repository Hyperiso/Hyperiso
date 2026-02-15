#include "ObservableExtractor.h"

static std::string obs_name(const ObservableId& obs) {
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
            cols.push_back(obs_key(obs, ord_));
            continue;
        }

        if (!vals[0].bin.has_value()) {
            cols.push_back(obs_key(obs, ord_));
            continue;
        }

        for (const auto& v : vals) {
            if (!v.bin.has_value()) {
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
            outY[obs_key(obs, ord_)] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        if (!vals[0].bin.has_value()) {
            outY[obs_key(obs, ord_)] = vals[0].value;
            continue;
        }

        for (const auto& v : vals) {
            if (!v.bin.has_value()) {
                outY[obs_key(obs, ord_)] = v.value;
                break;
            }
            const auto& br = v.bin.value();
            outY[obs_key_bin(obs, ord_, br.first, br.second)] = v.value;
        }
    }
}