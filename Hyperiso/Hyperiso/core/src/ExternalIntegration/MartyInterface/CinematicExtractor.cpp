#include "CinematicExtractor.h"

#include <algorithm>
#include <cctype>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace {

std::string read_file_to_string(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file) {
        return {};
    }

    std::ostringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

std::optional<std::string> first_braced_argument_list_after_compute(const std::string& content)
{
    const std::string needle = "computeWilsonCoefficients";
    const std::size_t call_pos = content.find(needle);
    if (call_pos == std::string::npos) {
        return std::nullopt;
    }

    const std::size_t first_brace = content.find('{', call_pos);
    if (first_brace == std::string::npos) {
        return std::nullopt;
    }

    int depth = 0;
    bool in_string = false;
    bool escaped = false;

    for (std::size_t i = first_brace; i < content.size(); ++i) {
        const char c = content[i];

        if (in_string) {
            if (escaped) {
                escaped = false;
            } else if (c == '\\') {
                escaped = true;
            } else if (c == '"') {
                in_string = false;
            }
            continue;
        }

        if (c == '"') {
            in_string = true;
            continue;
        }

        if (c == '{') {
            ++depth;
        } else if (c == '}') {
            --depth;
            if (depth == 0) {
                return content.substr(first_brace, i - first_brace + 1);
            }
        }
    }

    return std::nullopt;
}

} // namespace

int countMatchInRegex(std::string s, std::string re) {
    std::regex words_regex(re);
    auto words_begin = std::sregex_iterator(s.begin(), s.end(), words_regex);
    auto words_end = std::sregex_iterator();

    return std::distance(words_begin, words_end);
}

std::pair<int, int> CinematicExtractor::extract(const std::string& filename) {
    const auto process = extract_process(filename);
    return {static_cast<int>(process.incoming_count()), static_cast<int>(process.outgoing_count())};
}

CinematicProcess CinematicExtractor::extract_process(const std::string& filename) const {
    CinematicProcess process;
    const std::string content = read_file_to_string(filename);
    if (content.empty()) {
        return {};
    }

    // MARTY templates may wrap a physical external leg in helpers such as
    // AntiPart("mu") or Part("x").  Accept nested identifier wrappers before
    // the literal particle name instead of silently dropping that leg.
    const std::regex particle_regex(
        R"regex((Incoming|Outgoing)\s*\(\s*(?:[A-Za-z_][A-Za-z0-9_:]*\s*\(\s*)*"([^"]+)")regex"
    );

    auto parse_particle_list = [&](const std::string& particle_list) {
        CinematicProcess process;
        auto begin = std::sregex_iterator(particle_list.begin(), particle_list.end(), particle_regex);
        auto end = std::sregex_iterator();

        for (auto it = begin; it != end; ++it) {
            const std::string direction = (*it)[1].str();
            const std::string particle = normalize_particle_name((*it)[2].str());

            if (direction == "Incoming") {
                process.incoming.push_back(particle);
            } else {
                process.outgoing.push_back(particle);
            }
        }
        return process;
    };

    CinematicProcess best_process;
    auto keep_best = [&](CinematicProcess candidate) {
        const auto candidate_legs = candidate.incoming_count() + candidate.outgoing_count();
        const auto best_legs = best_process.incoming_count() + best_process.outgoing_count();
        if (candidate_legs > best_legs) {
            best_process = std::move(candidate);
        }
    };

    // Preferred candidate: a literal list in computeWilsonCoefficients(...).
    // Do not return a partial match immediately: a wrapped anti-particle used to
    // make a 1->3 process look like 1->2, while a later helper contained all legs.
    if (const auto particle_list = first_braced_argument_list_after_compute(content);
        particle_list.has_value()) {
        keep_best(parse_particle_list(*particle_list));
    }

    // Newer templates factor insertions in a helper and pass a variable to the
    // Wilson call. Scan all balanced blocks and retain the most complete process.
    for (std::size_t brace = content.find('{'); brace != std::string::npos;
         brace = content.find('{', brace + 1)) {
        int depth = 0;
        bool in_string = false;
        bool escaped = false;

        for (std::size_t i = brace; i < content.size(); ++i) {
            const char c = content[i];

            if (in_string) {
                if (escaped) {
                    escaped = false;
                } else if (c == '\\') {
                    escaped = true;
                } else if (c == '"') {
                    in_string = false;
                }
                continue;
            }

            if (c == '"') {
                in_string = true;
                continue;
            }

            if (c == '{') {
                ++depth;
            } else if (c == '}') {
                --depth;
                if (depth == 0) {
                    keep_best(parse_particle_list(content.substr(brace, i - brace + 1)));
                    break;
                }
            }
         }
    }

    return best_process;
}

std::string CinematicExtractor::normalize_particle_name(std::string particle_name) {
    particle_name.erase(std::remove_if(particle_name.begin(), particle_name.end(),
        [](unsigned char c) { return std::isspace(c); }), particle_name.end());

    if (!particle_name.empty() && particle_name.front() == '~') {
        particle_name.erase(particle_name.begin());
    }

    std::transform(particle_name.begin(), particle_name.end(), particle_name.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (particle_name == "gamma" || particle_name == "photon") {
        return "a";
    }

    return particle_name;
}

std::optional<std::string> CinematicExtractor::mass_symbol_for_particle(const std::string& particle_name) {
    const auto normalized = normalize_particle_name(particle_name);
    const auto& symbols = default_mass_symbols();

    const auto it = symbols.find(normalized);
    if (it == symbols.end()) {
        return std::nullopt;
    }
    return it->second;
}

const std::unordered_map<std::string, std::string>& CinematicExtractor::default_mass_symbols() {
    static const std::unordered_map<std::string, std::string> symbols {
        {"u", "m_u"}, {"d", "m_d"}, {"s", "m_s"},
        {"c", "m_c"}, {"b", "m_b"}, {"t", "m_t"},

        {"e", "m_e"}, {"mu", "m_mu"}, {"tau", "m_tau"},
        {"ve", "0"}, {"vmu", "0"}, {"vtau", "0"},

        {"a", "0"}, {"g", "0"},
        {"w", "m_W"}, {"w+", "m_W"}, {"w-", "m_W"},
        {"z", "m_Z"}, {"h", "m_h"}
    };

    return symbols;
}
