#include "StatisticProgress.h"

#include <algorithm>

namespace {

std::string format_seconds(double seconds) {
    std::ostringstream oss;
    if (seconds < 60.0) {
        oss << std::fixed << std::setprecision(1) << seconds << "s";
    } else if (seconds < 3600.0) {
        const int minutes = static_cast<int>(seconds / 60.0);
        const int sec = static_cast<int>(seconds) % 60;
        oss << minutes << "m" << std::setw(2) << std::setfill('0')
            << sec << std::setfill(' ') << "s";
    } else {
        const int hours = static_cast<int>(seconds / 3600.0);
        const int minutes = (static_cast<int>(seconds) / 60) % 60;
        const int sec = static_cast<int>(seconds) % 60;
        oss << hours << "h"
            << std::setw(2) << std::setfill('0') << minutes << "m"
            << std::setw(2) << std::setfill('0') << sec
            << std::setfill(' ') << "s";
    }
    return oss.str();
}

} // namespace

void StatisticProgressMonitor::reset(const std::string& phase, const std::string& message) {
    std::lock_guard<std::mutex> lock(mutex_);
    latest_ = {};
    latest_.phase = phase;
    latest_.message = message;
    latest_.sequence = next_sequence_++;
}

void StatisticProgressMonitor::update(StatisticProgressEvent event) {
    std::lock_guard<std::mutex> lock(mutex_);
    event.fraction = std::clamp(event.fraction, 0.0, 1.0);
    event.sequence = next_sequence_++;
    latest_ = std::move(event);
}

StatisticProgressEvent StatisticProgressMonitor::snapshot() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return latest_;
}

StatisticProgressReporter::StatisticProgressReporter(bool enabled,
                                                     std::size_t total,
                                                     std::size_t probe_draws,
                                                     std::size_t update_every,
                                                     std::ostream& os,
                                                     std::string label,
                                                     std::string item_label,
                                                     std::shared_ptr<StatisticProgressMonitor> monitor,
                                                     std::string phase)
    : enabled_(enabled), total_(total), probe_draws_(probe_draws == 0 ? 1 : probe_draws),
      update_every_(update_every == 0 ? 1 : update_every), os_(&os),
      start_(std::chrono::steady_clock::now()), label_(std::move(label)),
      item_label_(std::move(item_label)), monitor_(std::move(monitor)), phase_(std::move(phase))
{
    publish(0, 0, 0, false, "Starting Monte-Carlo sampling");
}

double StatisticProgressReporter::elapsed_seconds() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
}

double StatisticProgressReporter::eta_seconds(std::size_t accepted_count) const {
    if (accepted_count == 0 || accepted_count < probe_draws_ || accepted_count >= total_) return -1.0;
    return (elapsed_seconds() / static_cast<double>(accepted_count)) * static_cast<double>(total_ - accepted_count);
}

std::string StatisticProgressReporter::bar(std::size_t accepted_count) const {
    constexpr std::size_t width = 28;
    const double raw_frac = total_ == 0 ? 1.0 : static_cast<double>(accepted_count) / static_cast<double>(total_);
    const double frac = std::clamp(raw_frac, 0.0, 1.0);
    const std::size_t filled = static_cast<std::size_t>(frac * width);
    std::string out("[");
    for (std::size_t i = 0; i < width; ++i) out.push_back(i < filled ? '#' : '.');
    out.push_back(']');
    return out;
}

std::string StatisticProgressReporter::eta_string(std::size_t accepted_count) const {
    const double eta = eta_seconds(accepted_count);
    return eta < 0.0 ? "estimating" : format_seconds(eta);
}

void StatisticProgressReporter::publish(std::size_t accepted_count,
                                        std::size_t attempts,
                                        std::size_t failures,
                                        bool finished,
                                        const std::string& message) {
    if (!monitor_) return;
    StatisticProgressEvent ev;
    ev.phase = phase_;
    ev.message = message;
    ev.completed = accepted_count;
    ev.total = total_;
    ev.attempts = attempts;
    ev.failures = failures;
    ev.fraction = total_ == 0 ? 1.0 : static_cast<double>(accepted_count) / static_cast<double>(total_);
    ev.elapsed_seconds = elapsed_seconds();
    ev.eta_seconds = finished ? 0.0 : eta_seconds(accepted_count);
    ev.finished = finished;
    monitor_->update(std::move(ev));
}

void StatisticProgressReporter::accepted(std::size_t accepted_count, std::size_t attempts, std::size_t failures) {
    if (!enabled_ && !monitor_) return;
    if (total_ == 0) return;
    if (accepted_count != total_ && accepted_count % update_every_ != 0 && accepted_count != probe_draws_) return;

    if (attempts == 0) attempts = accepted_count;
    publish(accepted_count, attempts, failures, false,
            label_ + ": " + std::to_string(accepted_count) + "/" + std::to_string(total_) + " " + item_label_);
    if (!enabled_) return;
    const double pct = 100.0 * static_cast<double>(accepted_count) / static_cast<double>(total_);
    (*os_) << "\r\033[K[" << label_ << "] " << bar(accepted_count)
           << " " << accepted_count << "/" << total_ << " " << item_label_
           << " (" << std::fixed << std::setprecision(1) << pct << "%)"
           << " | ETA " << eta_string(accepted_count) << std::flush;
}

void StatisticProgressReporter::finish(std::size_t accepted_count, std::size_t attempts, std::size_t failures) {
    publish(accepted_count, attempts, failures, true,
            label_ + " completed: " + std::to_string(accepted_count) + "/" + std::to_string(total_) + " " + item_label_);
    if (!enabled_) return;
    (*os_) << "\r\033[K[" << label_ << "] " << bar(accepted_count)
           << " " << accepted_count << "/" << total_ << " " << item_label_
           << " completed in " << format_seconds(elapsed_seconds())
           << " (attempts=" << attempts << ", rejected=" << failures << ")\n";
}

StatisticStageProgressReporter::StatisticStageProgressReporter(bool enabled,
                                                               std::string label,
                                                               std::size_t total_steps,
                                                               std::ostream& os,
                                                               std::shared_ptr<StatisticProgressMonitor> monitor,
                                                               std::string phase)
    : enabled_(enabled), label_(std::move(label)), total_steps_(total_steps == 0 ? 1 : total_steps),
      os_(&os), start_(std::chrono::steady_clock::now()), monitor_(std::move(monitor)), phase_(std::move(phase)) {}

std::string StatisticStageProgressReporter::bar(std::size_t completed_steps) const {
    constexpr std::size_t width = 28;
    const double frac = std::clamp(static_cast<double>(completed_steps) / static_cast<double>(total_steps_), 0.0, 1.0);
    const std::size_t filled = static_cast<std::size_t>(frac * width);
    std::string out("[");
    for (std::size_t i = 0; i < width; ++i) out.push_back(i < filled ? '#' : '.');
    out.push_back(']');
    return out;
}

double StatisticStageProgressReporter::elapsed_seconds() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
}

std::string StatisticStageProgressReporter::elapsed_string() const { return format_seconds(elapsed_seconds()); }

void StatisticStageProgressReporter::publish(std::size_t completed_steps, const std::string& message, bool finished) {
    if (!monitor_) return;
    StatisticProgressEvent ev;
    ev.phase = phase_;
    ev.message = message;
    ev.completed = std::min(completed_steps, total_steps_);
    ev.total = total_steps_;
    ev.fraction = static_cast<double>(ev.completed) / static_cast<double>(total_steps_);
    ev.elapsed_seconds = elapsed_seconds();
    // Stage durations (covariance assembly, inversion, minimization, contours)
    // are not homogeneous, so extrapolating from completed stage count would be
    // misleading.  Keep the measured ETA exclusively for Monte-Carlo sampling.
    ev.eta_seconds = finished ? 0.0 : -1.0;
    ev.finished = finished;
    monitor_->update(std::move(ev));
}

void StatisticStageProgressReporter::start(const std::string& message) {
    publish(0, message, false);
    if (enabled_) (*os_) << "[" << label_ << "] " << message << '\n';
}

void StatisticStageProgressReporter::step(std::size_t completed_steps, const std::string& message) {
    publish(completed_steps, message, false);
    if (!enabled_) return;
    const std::size_t shown = std::min(completed_steps, total_steps_);
    (*os_) << "[" << label_ << "] " << bar(shown) << " " << shown << "/" << total_steps_
           << " | " << message << " | elapsed " << elapsed_string() << '\n';
}

void StatisticStageProgressReporter::finish(const std::string& message) {
    publish(total_steps_, message, true);
    if (enabled_) (*os_) << "[" << label_ << "] " << bar(total_steps_) << " " << total_steps_ << "/" << total_steps_
                         << " | " << message << " | total " << elapsed_string() << '\n';
}
    