#include "StatisticProgress.h"

StatisticProgressReporter::StatisticProgressReporter(bool enabled,
                                                     std::size_t total,
                                                     std::size_t probe_draws,
                                                     std::size_t update_every,
                                                     std::ostream& os)
    : enabled_(enabled),
      total_(total),
      probe_draws_(probe_draws == 0 ? 1 : probe_draws),
      update_every_(update_every == 0 ? 1 : update_every),
      os_(&os),
      start_(std::chrono::steady_clock::now())
{}

std::string StatisticProgressReporter::bar(std::size_t accepted_count) const {
    constexpr std::size_t width = 28;
    const double frac = total_ == 0 ? 1.0 : static_cast<double>(accepted_count) / static_cast<double>(total_);
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

    std::ostringstream oss;
    if (remaining < 60.0) {
        oss << std::fixed << std::setprecision(1) << remaining << "s";
    } else if (remaining < 3600.0) {
        const int minutes = static_cast<int>(remaining / 60.0);
        const int seconds = static_cast<int>(remaining) % 60;
        oss << minutes << "m" << std::setw(2) << std::setfill('0')
            << seconds << std::setfill(' ') << "s";
    } else {
        const int hours = static_cast<int>(remaining / 3600.0);
        const int minutes = (static_cast<int>(remaining) / 60) % 60;
        const int seconds = static_cast<int>(remaining) % 60;
        oss << hours << "h"
            << std::setw(2) << std::setfill('0') << minutes << "m"
            << std::setw(2) << std::setfill('0') << seconds
            << std::setfill(' ') << "s";
    }
    return oss.str();
}

void StatisticProgressReporter::accepted(std::size_t accepted_count) {
    if (!enabled_ || total_ == 0) return;
    if (accepted_count != total_ && accepted_count % update_every_ != 0 && accepted_count != probe_draws_) return;

    const double pct = 100.0 * static_cast<double>(accepted_count) / static_cast<double>(total_);
    (*os_) << "\r\033[K[MC] " << bar(accepted_count)
           << " " << accepted_count << "/" << total_
           << " (" << std::fixed << std::setprecision(1) << pct << "%)"
           << " ETA " << eta_string(accepted_count)
           << std::flush;
}

void StatisticProgressReporter::finish(std::size_t accepted_count, std::size_t attempts, std::size_t failures) {
    if (!enabled_) return;
    const auto now = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(now - start_).count();
    (*os_) << "\r\033[K[MC] " << bar(accepted_count)
           << " " << accepted_count << "/" << total_
           << " done in " << std::fixed << std::setprecision(2) << elapsed << "s"
           << " (attempts=" << attempts << ", rejected=" << failures << ")\n";
}
