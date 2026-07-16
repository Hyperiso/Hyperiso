#ifndef STATISTIC_PROGRESS_H
#define STATISTIC_PROGRESS_H

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>

/** Snapshot shared with non-terminal frontends (Python/Dash, notebooks, etc.). */
struct StatisticProgressEvent {
    std::string phase = "idle";
    std::string message;
    std::size_t completed = 0;
    std::size_t total = 0;
    std::size_t attempts = 0;
    std::size_t failures = 0;
    double fraction = 0.0;
    double elapsed_seconds = 0.0;
    double eta_seconds = -1.0; // negative means that no reliable ETA is available
    bool finished = false;
    std::uint64_t sequence = 0;
};

/** Thread-safe latest-value monitor used by graphical frontends. */
class StatisticProgressMonitor {
public:
    void reset(const std::string& phase = "preparing", const std::string& message = "Preparing statistic workflow");
    void update(StatisticProgressEvent event);
    StatisticProgressEvent snapshot() const;

private:
    mutable std::mutex mutex_;
    StatisticProgressEvent latest_ {};
    std::uint64_t next_sequence_ = 1;
};

/** Terminal MC progress reporter with measured ETA and optional GUI monitor. */
class StatisticProgressReporter {
public:
    StatisticProgressReporter(bool enabled,
                              std::size_t total,
                              std::size_t probe_draws = 5,
                              std::size_t update_every = 1,
                              std::ostream& os = std::cout,
                              std::string label = "Monte-Carlo",
                              std::string item_label = "accepted",
                              std::shared_ptr<StatisticProgressMonitor> monitor = nullptr,
                              std::string phase = "monte_carlo");

    void accepted(std::size_t accepted_count, std::size_t attempts = 0, std::size_t failures = 0);
    void finish(std::size_t accepted_count, std::size_t attempts, std::size_t failures);

private:
    double elapsed_seconds() const;
    double eta_seconds(std::size_t accepted_count) const;
    std::string eta_string(std::size_t accepted_count) const;
    std::string bar(std::size_t accepted_count) const;
    void publish(std::size_t accepted_count,
                 std::size_t attempts,
                 std::size_t failures,
                 bool finished,
                 const std::string& message);

    bool enabled_ = false;
    std::size_t total_ = 0;
    std::size_t probe_draws_ = 5;
    std::size_t update_every_ = 1;
    std::ostream* os_ = nullptr;
    std::chrono::steady_clock::time_point start_;
    std::string label_ = "Monte-Carlo";
    std::string item_label_ = "accepted";
    std::shared_ptr<StatisticProgressMonitor> monitor_;
    std::string phase_ = "monte_carlo";
};

/** Stage-based progress reporter for covariance / likelihood / fit steps. */
class StatisticStageProgressReporter {
public:
    StatisticStageProgressReporter(bool enabled,
                                   std::string label,
                                   std::size_t total_steps,
                                   std::ostream& os = std::cout,
                                   std::shared_ptr<StatisticProgressMonitor> monitor = nullptr,
                                   std::string phase = "chi2_pipeline");

    void start(const std::string& message);
    void step(std::size_t completed_steps, const std::string& message);
    void finish(const std::string& message);

private:
    std::string bar(std::size_t completed_steps) const;
    double elapsed_seconds() const;
    std::string elapsed_string() const;
    void publish(std::size_t completed_steps, const std::string& message, bool finished);

    bool enabled_ = false;
    std::string label_;
    std::size_t total_steps_ = 0;
    std::ostream* os_ = nullptr;
    std::chrono::steady_clock::time_point start_;
    std::shared_ptr<StatisticProgressMonitor> monitor_;
    std::string phase_ = "chi2_pipeline";
};

#endif
