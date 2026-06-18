#include "WilsonHandler.h"

#include <iostream>
#include <unordered_set>
#include <vector>

#include "CliUtils.h"
#include "WilsonInterface.h"
#include "mapper_hub.hpp"

namespace {

void print_wilson_usage() {
    std::cout
        << "Usage:\n"
        << "  hyperiso-ui wilson summary [options]\n\n"
        << "Options:\n"
        << "  --groups <csv>     Wilson groups, default BCoefficients\n"
        << "  --coeffs <csv>     Coefficients, default C7,C9,C10\n"
        << "  --qmatch <value>   Matching scale, default 81\n"
        << "  --q <value>        Running scale, default 4.8\n"
        << "  --order <order>    LO, NLO or NNLO, default NNLO\n"
        << "  --model <model>    SM, THDM, MSSM or MARTY, default SM\n"
        << "  --lha <path>       Input LHA/FLHA file\n";
}

} // namespace

int handleWilsonOptions(int argc, char* argv[]) {
    CliOptions opts = CliOptions::parse(argc, argv, 1);
    const std::string command = opts.positionals.empty() ? "summary" : opts.positionals[0];

    if (opts.flag("help", false) || command == "help") {
        print_wilson_usage();
        return 0;
    }
    if (command != "summary") {
        throw std::invalid_argument("Unknown wilson command: " + command);
    }

    auto hyp = init_hyperiso_from_cli(opts);
    init_all_builtins();

    const auto group_names = opts.list("groups", {"BCoefficients"});
    const auto coeff_names = opts.list("coeffs", {"C7", "C9", "C10"});
    const double qmatch = opts.get_double("qmatch", 81.0);
    const double q = opts.get_double("q", 4.8);
    const QCDOrder order = parse_qcd_order(opts.get("order", "NNLO"));

    std::unordered_set<WGroupId> groups;
    for (const auto& name : group_names) {
        groups.insert(GroupMapper::id_of(name));
    }

    WilsonInterface wilson;
    wilson.build(WilsonBuildConfig(groups, qmatch, q, order));

    print_section("Wilson summary");
    std::cout << "order=" << OrderMapper::str(order)
              << ", qmatch=" << qmatch
              << ", q=" << q << "\n";

    for (const auto& group_name : group_names) {
        WGroupId group = GroupMapper::id_of(group_name);
        std::cout << "\nGroup: " << GroupMapper::str(group) << "\n";
        for (const auto& coeff_name : coeff_names) {
            WCoefId coeff = WCoefMapper::id_of(coeff_name);
            const auto m = wilson.getFM(group, coeff, order, ContributionType::TOTAL);
            const auto r = wilson.getFR(group, coeff, order, ContributionType::TOTAL);
            std::cout << "  " << WCoefMapper::str(coeff)
                      << "  matching=" << m
                      << "  running=" << r << "\n";
        }
    }

    return 0;
}
