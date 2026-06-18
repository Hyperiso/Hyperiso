#ifndef HYPERISO_CLI_UTILS_H
#define HYPERISO_CLI_UTILS_H

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "Include.h"
#include "HyperisoMaster.h"

/**
 * @struct CliOptions
 * @brief Small command-line parser for the Hyperiso terminal interface.
 *
 * The interface deliberately keeps parsing simple: options are passed as
 * ``--key value``, ``--key=value`` or boolean flags.  Comma-separated values are
 * handled by @ref list().  This avoids the heavy historical parser and makes CLI
 * handlers easy to read and maintain.
 */
struct CliOptions {
    std::vector<std::string> positionals;
    std::map<std::string, std::vector<std::string>> options;

    static CliOptions parse(int argc, char* argv[], int first = 1);

    bool has(const std::string& key) const;
    std::string get(const std::string& key, const std::string& fallback = "") const;
    int get_int(const std::string& key, int fallback) const;
    double get_double(const std::string& key, double fallback) const;
    bool flag(const std::string& key, bool fallback = false) const;
    std::vector<std::string> list(const std::string& key, const std::vector<std::string>& fallback = {}) const;
};

/** Print a title surrounded by a light separator. */
void print_section(const std::string& title);

/** Convert a CLI string to a Hyperiso model enum. */
Model parse_model(const std::string& model);

/** Convert a CLI string to QCD order. */
QCDOrder parse_qcd_order(const std::string& order);

/** Build and initialize HyperisoMaster from common CLI options. */
HyperisoMaster init_hyperiso_from_cli(const CliOptions& opts);

#endif
