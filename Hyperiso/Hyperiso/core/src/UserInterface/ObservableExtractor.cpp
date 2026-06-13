#include "ObservableExtractor.h"

#include <limits>
#include <set>

namespace {

static std::string obs_name(const ObservableId& obs) {
    return obs.str();
}

static bool has_specific_bin(const BinnedObservableId& obs) {
    // Dans le core actuel, {0,0} est utilisé comme convention "unbinned/default".
    return obs.p.first != 0.0 || obs.p.second != 0.0;
}

static void append_schema_for_values(std::vector<std::string>& cols,
                                     std::set<std::string>& seen,
                                     const ObservableId& obs,
                                     QCDOrder ord,
                                     const std::vector<ObservableValue>& vals)
{
    auto add_col = [&](std::string key) {
        if (seen.insert(key).second) {
            cols.push_back(std::move(key));
        }
    };

    if (vals.empty()) {
        add_col(obs_key(obs, ord));
        return;
    }

    if (!vals[0].bin.has_value()) {
        add_col(obs_key(obs, ord));
        return;
    }

    for (const auto& v : vals) {
        if (!v.bin.has_value()) {
            add_col(obs_key(obs, ord));
            continue;
        }

        const auto& br = v.bin.value();
        add_col(obs_key_bin(obs, ord, br.first, br.second));
    }
}

static void write_values(std::unordered_map<std::string, Value>& outY,
                         const ObservableId& obs,
                         QCDOrder ord,
                         const std::vector<ObservableValue>& vals)
{
    if (vals.empty()) {
        outY[obs_key(obs, ord)] = std::numeric_limits<double>::quiet_NaN();
        return;
    }

    if (!vals[0].bin.has_value()) {
        outY[obs_key(obs, ord)] = vals[0].value;
        return;
    }

    for (const auto& v : vals) {
        if (!v.bin.has_value()) {
            outY[obs_key(obs, ord)] = v.value;
            continue;
        }

        const auto& br = v.bin.value();
        outY[obs_key_bin(obs, ord, br.first, br.second)] = v.value;
    }
}

} // namespace

std::string obs_key(const ObservableId& obs, QCDOrder ord) {
    std::string k = "OBS:";
    k += obs_name(obs);
    k += ":";
    k += OrderMapper::str(ord);
    return k;
}

std::string obs_key_bin(const ObservableId& obs,
                        QCDOrder ord,
                        double bmin,
                        double bmax)
{
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
    std::set<std::string> seen;

    oi_.enable_obs();

    auto current = oi_.get_current_observables();
    std::sort(current.begin(), current.end());

    for (const auto& obs : current) {
        if (has_specific_bin(obs)) {
            const auto key = obs_key_bin(obs.s, ord_, obs.p.first, obs.p.second);
            if (seen.insert(key).second) {
                cols.push_back(key);
            }
            continue;
        }

        auto vals = oi_.compute_observable(obs.s);
        append_schema_for_values(cols, seen, obs.s, ord_, vals);
    }

    return cols;
}

void ObservableExtractor::extract(std::unordered_map<std::string, Value>& outY,
                                  const OutputSpec& /*spec*/) const
{
    oi_.enable_obs();

    auto current = oi_.get_current_observables();
    std::sort(current.begin(), current.end());

    for (const auto& obs : current) {
        if (has_specific_bin(obs)) {
            auto v = oi_.compute_observable(obs);

            const auto br = v.bin.value_or(obs.p);
            outY[obs_key_bin(obs.s, ord_, br.first, br.second)] = v.value;
            continue;
        }

        auto vals = oi_.compute_observable(obs.s);
        write_values(outY, obs.s, ord_, vals);
    }
}