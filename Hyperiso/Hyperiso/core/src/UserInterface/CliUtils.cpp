#include "CliUtils.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

bool is_option(const std::string& s) {
    return s.rfind("--", 0) == 0;
}

} // namespace

CliOptions CliOptions::parse(int argc, char* argv[], int first) {
    CliOptions out;
    for (int i = first; i < argc; ++i) {
        std::string token = argv[i];
        if (!is_option(token)) {
            out.positionals.push_back(token);
            continue;
        }

        token = token.substr(2);
        std::string key;
        std::string value;

        const auto eq = token.find('=');
        if (eq != std::string::npos) {
            key = token.substr(0, eq);
            value = token.substr(eq + 1);
        } else {
            key = token;
            if (i + 1 < argc && !is_option(argv[i + 1])) {
                value = argv[++i];
            } else {
                value = "true";
            }
        }
        out.options[key].push_back(value);
    }
    return out;
}

bool CliOptions::has(const std::string& key) const {
    return options.contains(key);
}

std::string CliOptions::get(const std::string& key, const std::string& fallback) const {
    auto it = options.find(key);
    if (it == options.end() || it->second.empty()) return fallback;
    return it->second.back();
}

int CliOptions::get_int(const std::string& key, int fallback) const {
    if (!has(key)) return fallback;
    return std::stoi(get(key));
}

double CliOptions::get_double(const std::string& key, double fallback) const {
    if (!has(key)) return fallback;
    return std::stod(get(key));
}

bool CliOptions::flag(const std::string& key, bool fallback) const {
    if (!has(key)) return fallback;
    const std::string v = lower(get(key, "true"));
    return v == "true" || v == "1" || v == "yes" || v == "on";
}

std::vector<std::string> CliOptions::list(const std::string& key, const std::vector<std::string>& fallback) const {
    auto it = options.find(key);
    if (it == options.end() || it->second.empty()) return fallback;

    std::vector<std::string> out;
    for (const auto& raw : it->second) {
        auto parts = split_csv(raw);
        out.insert(out.end(), parts.begin(), parts.end());
    }
    return out.empty() ? fallback : out;
}

void print_section(const std::string& title) {
    std::cout << "\n== " << title << " ==\n";
}

Model parse_model(const std::string& model) {
    const std::string m = lower(model);
    if (m == "sm") return Model::SM;
    if (m == "thdm" || m == "2hdm") return Model::THDM;
    if (m == "mssm" || m == "susy" || m == "nmssm") return Model::SUSY;
    if (m == "marty") return Model::MARTY;
    throw std::invalid_argument("Unknown model: " + model);
}

QCDOrder parse_qcd_order(const std::string& order) {
    const std::string o = lower(order);
    if (o == "lo") return QCDOrder::LO;
    if (o == "nlo") return QCDOrder::NLO;
    if (o == "nnlo") return QCDOrder::NNLO;
    throw std::invalid_argument("Unknown QCD order: " + order);
}


ContributionType parse_contribution(const std::string& contribution) {
    const std::string value = lower(contribution);
    if (value == "sm") return ContributionType::SM;
    if (value == "bsm") return ContributionType::BSM;
    if (value == "total") return ContributionType::TOTAL;
    throw std::invalid_argument("Unknown contribution type: " + contribution);
}

HyperisoMaster init_hyperiso_from_cli(const CliOptions& opts) {
    HyperisoConfig cfg;
    cfg.model = parse_model(opts.get("model", "SM"));
    cfg.flags[ExternalFlag::HYP_AS_SM_MARTY] = opts.flag("sm-marty", false);
    cfg.flags[ExternalFlag::IS_LHA_SPECTRUM] = opts.flag("spectrum", false);

    HyperisoMaster hyp;
    hyp.init(opts.get("lha", "lha/si_input.flha"), cfg);
    return hyp;
}
