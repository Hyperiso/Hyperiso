#ifndef CINEMATIC_EXTRACTOR_H
#define CINEMATIC_EXTRACTOR_H

#include <string>
#include <utility>
#include <regex>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <optional>

/**
 * @brief Counts the number of occurrences of a regex pattern in a given string.
 */
int countMatchInRegex(std::string s, std::string re);

/**
 * @struct CinematicProcess
 * @brief Incoming/outgoing particles extracted from a MARTY computeWilsonCoefficients call.
 *
 * The ordered leg convention used by the kinematic formulas is:
 *   m1 = first incoming particle,
 *   m2, m3, ... = outgoing particles in the order they appear in the template.
 */
struct CinematicProcess {
    std::vector<std::string> incoming;
    std::vector<std::string> outgoing;

    bool empty() const { return incoming.empty() && outgoing.empty(); }
    std::size_t incoming_count() const { return incoming.size(); }
    std::size_t outgoing_count() const { return outgoing.size(); }

    std::vector<std::string> ordered_particles() const {
        std::vector<std::string> out;
        out.reserve(incoming.size() + outgoing.size());
        out.insert(out.end(), incoming.begin(), incoming.end());
        out.insert(out.end(), outgoing.begin(), outgoing.end());
        return out;
    }
};

/**
 * @class CinematicExtractor
 * @brief Extracts process kinematics from MARTY templates.
 *
 * The extractor looks for the first computeWilsonCoefficients(...) call and
 * parses the particle list entries Incoming("...") and Outgoing("...").
 * If the call appears multiple times, only the first one is used.
 */
class CinematicExtractor {
public:
    /**
     * @brief Legacy API: returns {number_of_incoming, number_of_outgoing}.
     */
    std::pair<int, int> extract(const std::string& filename);

    /**
     * @brief Extracts the first MARTY process from a template file.
     */
    CinematicProcess extract_process(const std::string& filename) const;

    /**
     * @brief Normalizes common MARTY particle names for dictionary lookup.
     */
    static std::string normalize_particle_name(std::string particle_name);

    /**
     * @brief Returns the MARTY mass symbol associated with a particle name.
     *
     * Examples: b -> m_b, s -> m_s, A/g/gamma -> 0.
     */
    static std::optional<std::string> mass_symbol_for_particle(const std::string& particle_name);

    /**
     * @brief Default particle -> MARTY mass-symbol dictionary.
     */
    static const std::unordered_map<std::string, std::string>& default_mass_symbols();
};

#endif // CINEMATIC_EXTRACTOR_H
