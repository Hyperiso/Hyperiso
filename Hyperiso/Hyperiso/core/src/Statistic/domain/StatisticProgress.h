#ifndef STATISTIC_PROGRESS_H
#define STATISTIC_PROGRESS_H

#include <chrono>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

/**
 * @class StatisticProgressReporter
 * @brief Lightweight terminal progress reporter with live ETA estimates.
 *
 * The reporter is intentionally independent from StatisticManager.  It measures
 * the wall-clock time of accepted Monte Carlo draws and prints compact progress
 * lines only when explicitly enabled.  The first few accepted draws are used as a
 * warm-up window so that the ETA is based on measured runtime instead of a fixed
 * guess.
 */
class StatisticProgressReporter {
public:
    StatisticProgressReporter(bool enabled,
                              std::size_t total,
                              std::size_t probe_draws = 5,
                              std::size_t update_every = 1,
                              std::ostream& os = std::cout);

    /** Notify the reporter that one accepted item was produced. */
    void accepted(std::size_t accepted_count);

    /** Print a final newline/summary when progress was enabled. */
    void finish(std::size_t accepted_count, std::size_t attempts, std::size_t failures);

private:
    std::string eta_string(std::size_t accepted_count) const;
    std::string bar(std::size_t accepted_count) const;

    bool enabled_ = false;
    std::size_t total_ = 0;
    std::size_t probe_draws_ = 5;
    std::size_t update_every_ = 1;
    std::ostream* os_ = nullptr;
    std::chrono::steady_clock::time_point start_;
};

#endif
