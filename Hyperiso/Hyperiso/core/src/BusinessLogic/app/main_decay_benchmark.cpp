#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#include "HyperisoMaster.h"
#include "Logger.h"
#include "ObservableInterface.h"

namespace {

struct Options {
    std::string input = "lha/si_input.flha";
    std::string output = "benchmark_decays.csv";
    int repeats = 30;
    int warmups = 3;
    unsigned int threads = 0; // 0 -> hardware_concurrency for decays that support it
    QCDOrder order = QCDOrder::NLO;
};

struct DecayCase {
    Decays decay;
    bool supports_threads = false;
};

QCDOrder parse_order(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    if (value == "LO") return QCDOrder::LO;
    if (value == "NLO") return QCDOrder::NLO;
    if (value == "NNLO") return QCDOrder::NNLO;
    if (value == "NONE") return QCDOrder::NONE;
    LOG_ERROR("ArgumentError", "Unknown QCD order:", value, "expected LO, NLO, NNLO or NONE.");
    return QCDOrder::NONE;
}

void print_usage(const char* exe) {
    std::cerr
        << "Usage: " << exe << " [options]\n"
        << "\nOptions:\n"
        << "  --input PATH     Input FLHA/LHA file (default: default/lha/testInput.flha)\n"
        << "  --out PATH       Output CSV file (default: benchmark_decays.csv)\n"
        << "  --repeats N      Timed repetitions per decay (default: 30)\n"
        << "  --warmups N      Untimed warm-up calls per decay (default: 3)\n"
        << "  --threads N      Threads for decays exposing set_*_threads; 0 uses hardware_concurrency (default: 0)\n"
        << "  --order ORDER    QCD order: LO, NLO, NNLO or NONE (default: NLO)\n"
        << "  --help           Show this message\n";
}

Options parse_args(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) {
                LOG_ERROR("ArgumentError", "Missing value after", name);
            }
            return argv[++i];
        };

        if (arg == "--input") opt.input = require_value(arg);
        else if (arg == "--out") opt.output = require_value(arg);
        else if (arg == "--repeats") opt.repeats = std::stoi(require_value(arg));
        else if (arg == "--warmups") opt.warmups = std::stoi(require_value(arg));
        else if (arg == "--threads") opt.threads = static_cast<unsigned int>(std::stoul(require_value(arg)));
        else if (arg == "--order") opt.order = parse_order(require_value(arg));
        else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            LOG_ERROR("ArgumentError", "Unknown option:", arg);
        }
    }

    if (opt.repeats <= 0) LOG_ERROR("ArgumentError", "--repeats must be positive.");
    if (opt.warmups < 0) LOG_ERROR("ArgumentError", "--warmups must be non-negative.");

    if (opt.threads == 0) {
        opt.threads = std::thread::hardware_concurrency();
        if (opt.threads == 0) opt.threads = 1;
    }
    return opt;
}

double checksum(const std::map<ObservableId, std::vector<ObservableValue>>& results) {
    double sum = 0.0;
    for (const auto& [_, values] : results) {
        for (const auto& value : values) {
            sum += value.value;
        }
    }
    return sum;
}

struct Stats {
    double mean_ms = 0.0;
    double stddev_ms = 0.0;
    double min_ms = 0.0;
    double max_ms = 0.0;
};

Stats compute_stats(const std::vector<double>& samples) {
    Stats stats;
    stats.mean_ms = std::accumulate(samples.begin(), samples.end(), 0.0) / static_cast<double>(samples.size());
    stats.min_ms = *std::min_element(samples.begin(), samples.end());
    stats.max_ms = *std::max_element(samples.begin(), samples.end());

    double variance = 0.0;
    for (double x : samples) {
        variance += (x - stats.mean_ms) * (x - stats.mean_ms);
    }
    variance /= static_cast<double>(samples.size());
    stats.stddev_ms = std::sqrt(variance);
    return stats;
}

std::vector<DecayCase> default_decay_cases() {
    return {
        {Decays::B__D_l_nu, false},
        {Decays::B__Dstar_l_nu, false},
        {Decays::B__Kstar_gamma, false},
        {Decays::B__l_l, false},
        {Decays::B__l_nu, false},
        {Decays::B__Xs_gamma, false},
        {Decays::B__Xs_l_l, false},
        {Decays::M0_Mix, false},
        {Decays::B__Kstar_l_l, true},
        {Decays::B__K_l_l, true},
        {Decays::Bs__phi_l_l, true},
        {Decays::Lambda_b__Lambda_l_l, false},
        {Decays::K__l_l, false},
        {Decays::K__pi_nu_nu, false},
        {Decays::K__l_nu, false},
        {Decays::D__l_nu, false},
        {Decays::Ds__l_nu, false}
    };
}

void apply_thread_setting(ObservableInterface& interface, Decays decay, unsigned int threads) {
    switch (decay) {
        case Decays::B__Kstar_l_l:
            interface.set_bkstarll_threads(threads);
            break;
        case Decays::B__K_l_l:
            interface.set_bkll_threads(threads);
            break;
        case Decays::Bs__phi_l_l:
            interface.set_bsphi_threads(threads);
            break;
        default:
            break;
    }
}

} // namespace

int main(int argc, char** argv) {
    Logger::getInstance()->setLevel(Logger::LogLevel::WARN);
    const Options opt = parse_args(argc, argv);

    HyperisoMaster hyperiso;
    HyperisoConfig config;
    config.model = Model::SM;
    hyperiso.init(opt.input, config);

    std::ofstream csv(opt.output);
    if (!csv) {
        LOG_ERROR("IOError", "Cannot open output CSV:", opt.output);
    }

    csv << "decay,order,threads,repeats,warmups,n_observables,mean_ms,stddev_ms,min_ms,max_ms,checksum\n";

    for (const auto& decay_case : default_decay_cases()) {
        ObservableInterface interface;
        interface.add_observables(decay_case.decay, opt.order, false);
        if (decay_case.supports_threads) {
            apply_thread_setting(interface, decay_case.decay, opt.threads);
        }

        for (int i = 0; i < opt.warmups; ++i) {
            volatile double sink = checksum(interface.compute_all());
            (void)sink;
        }

        std::vector<double> samples_ms;
        samples_ms.reserve(static_cast<size_t>(opt.repeats));
        double last_checksum = 0.0;
        size_t n_observables = 0;

        for (int i = 0; i < opt.repeats; ++i) {
            const auto start = std::chrono::steady_clock::now();
            auto results = interface.compute_all();
            const auto stop = std::chrono::steady_clock::now();

            n_observables = results.size();
            last_checksum = checksum(results);
            samples_ms.push_back(std::chrono::duration<double, std::milli>(stop - start).count());
        }

        const Stats stats = compute_stats(samples_ms);
        const std::string decay_name = DecayMapper::str(decay_case.decay);
        csv << decay_name << ','
            << OrderMapper::str(opt.order) << ','
            << (decay_case.supports_threads ? opt.threads : 1) << ','
            << opt.repeats << ','
            << opt.warmups << ','
            << n_observables << ','
            << stats.mean_ms << ','
            << stats.stddev_ms << ','
            << stats.min_ms << ','
            << stats.max_ms << ','
            << last_checksum << '\n';
        csv.flush();

        std::cout << "decay=" << decay_name
                  << " mean_ms=" << stats.mean_ms
                  << " stddev_ms=" << stats.stddev_ms
                  << " n_observables=" << n_observables << '\n';
    }

    std::cout << "Wrote " << opt.output << '\n';
    return 0;
}
