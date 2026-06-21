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

StatisticProgressReporter::StatisticProgressReporter(bool enabled,
                                                     std::size_t total,
                                                     std::size_t probe_draws,
                                                     std::size_t update_every,
                                                     std::ostream& os,
                                                     std::string label,
                                                     std::string item_label)
    : enabled_(enabled),
      total_(total),
      probe_draws_(probe_draws == 0 ? 1 : probe_draws),
      update_every_(update_every == 0 ? 1 : update_every),
      os_(&os),
      start_(std::chrono::steady_clock::now()),
      label_(std::move(label)),
      item_label_(std::move(item_label))
{}

std::string StatisticProgressReporter::bar(std::size_t accepted_count) const {
    constexpr std::size_t width = 28;
    const double raw_frac = total_ == 0 ? 1.0 : static_cast<double>(accepted_count) / static_cast<double>(total_);
    const double frac = std::clamp(raw_frac, 0.0, 1.0);
    const std::size_t filled = static_cast<std::size_t>(frac * width);

    std::string out;
    out.reserve(width + 2);
    out.push_back('[');
    for (std::size_t i = 0; i < width; ++i) out.push_back(i < filled ? '#' : '.');
    out.push_back(']');
    return out;
}

std::string StatisticProgressReporter::eta_string(std::size_t accepted_count) const {
    if (accepted_count == 0 || accepted_count < probe_draws_) {
        return "estimating";
    }

    const auto now = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(now - start_).count();
    const double sec_per_draw = elapsed / static_cast<double>(accepted_count);
    const double remaining = sec_per_draw * static_cast<double>(total_ - accepted_count);

    return format_seconds(remaining);
}

void StatisticProgressReporter::accepted(std::size_t accepted_count) {
    if (!enabled_ || total_ == 0) return;
    if (accepted_count != total_ && accepted_count % update_every_ != 0 && accepted_count != probe_draws_) return;

    const double pct = 100.0 * static_cast<double>(accepted_count) / static_cast<double>(total_);
    (*os_) << "\r\033[K[" << label_ << "] " << bar(accepted_count)
           << " " << accepted_count << "/" << total_ << " " << item_label_
           << " (" << std::fixed << std::setprecision(1) << pct << "%)"
           << " | ETA " << eta_string(accepted_count)
           << std::flush;
}

void StatisticProgressReporter::finish(std::size_t accepted_count, std::size_t attempts, std::size_t failures) {
    if (!enabled_) return;
    const auto now = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(now - start_).count();
    (*os_) << "\r\033[K[" << label_ << "] " << bar(accepted_count)
           << " " << accepted_count << "/" << total_ << " " << item_label_
           << " completed in " << format_seconds(elapsed)
           << " (attempts=" << attempts << ", rejected=" << failures << ")\n";
}

StatisticStageProgressReporter::StatisticStageProgressReporter(bool enabled,
                                                               std::string label,
                                                               std::size_t total_steps,
                                                               std::ostream& os)
    : enabled_(enabled),
      label_(std::move(label)),
      total_steps_(total_steps == 0 ? 1 : total_steps),
      os_(&os),
      start_(std::chrono::steady_clock::now())
{}

std::string StatisticStageProgressReporter::bar(std::size_t completed_steps) const {
    constexpr std::size_t width = 28;
    const double raw_frac = static_cast<double>(completed_steps) / static_cast<double>(total_steps_);
    const double frac = std::clamp(raw_frac, 0.0, 1.0);
    const std::size_t filled = static_cast<std::size_t>(frac * width);

    std::string out;
    out.reserve(width + 2);
    out.push_back('[');
    for (std::size_t i = 0; i < width; ++i) out.push_back(i < filled ? '#' : '.');
    out.push_back(']');
    return out;
}

std::string StatisticStageProgressReporter::elapsed_string() const {
    const auto now = std::chrono::steady_clock::now();
    return format_seconds(std::chrono::duration<double>(now - start_).count());
}

void StatisticStageProgressReporter::start(const std::string& message) {
    if (!enabled_) return;
    (*os_) << "[" << label_ << "] " << message << '\n';
}

void StatisticStageProgressReporter::step(std::size_t completed_steps, const std::string& message) {
    if (!enabled_) return;
    const std::size_t shown = std::min(completed_steps, total_steps_);
    (*os_) << "[" << label_ << "] " << bar(shown)
           << " " << shown << "/" << total_steps_
           << " | " << message
           << " | elapsed " << elapsed_string() << '\n';
}

void StatisticStageProgressReporter::finish(const std::string& message) {
    if (!enabled_) return;
    (*os_) << "[" << label_ << "] " << bar(total_steps_)
           << " " << total_steps_ << "/" << total_steps_
           << " | " << message
           << " | total " << elapsed_string() << '\n';
}
