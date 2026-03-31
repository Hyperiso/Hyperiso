#ifndef BINNEDOBSERVABLEID_H
#define BINNEDOBSERVABLEID_H

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "observable_ids.hpp"

static inline double norm_zero(double x) noexcept {
    return (x == 0.0) ? 0.0 : x;
}

static inline std::uint64_t bits_norm_zero(double x) noexcept {
    x = norm_zero(x);
    return std::bit_cast<std::uint64_t>(x);
}

struct EncodedBin {
    long int_part;
    long frac_part;
    long frac_digits;
};

static inline long pow10_long(long n) noexcept {
    long v = 1;
    while (n-- > 0) {
        v *= 10;
    }
    return v;
}

static inline std::string trim_decimal_string(std::string s) {
    auto pos = s.find('.');
    if (pos == std::string::npos) {
        return s;
    }

    while (!s.empty() && s.back() == '0') {
        s.pop_back();
    }

    if (!s.empty() && s.back() == '.') {
        s.pop_back();
    }

    if (s.empty()) {
        return "0";
    }

    return s;
}

static inline EncodedBin encode_bin_gev(double x) {
    if (!std::isfinite(x)) {
        throw std::runtime_error("Bin value must be finite");
    }

    x = norm_zero(x);

    constexpr int kMaxFracDigits = 12;
    constexpr long double kScale = 1.0e12L;

    long double xr = static_cast<long double>(x);
    xr = std::round(xr * kScale) / kScale;

    const bool neg = std::signbit(static_cast<double>(xr));
    const long double ax = std::fabs(xr);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(kMaxFracDigits) << ax;

    std::string s = trim_decimal_string(oss.str());

    const auto pos = s.find('.');
    std::string int_str  = (pos == std::string::npos) ? s : s.substr(0, pos);
    std::string frac_str = (pos == std::string::npos) ? "" : s.substr(pos + 1);

    if (int_str.empty()) {
        int_str = "0";
    }

    const long int_part = std::stol(int_str);
    const long frac_part = frac_str.empty() ? 0L : std::stol(frac_str);
    const long frac_digits = static_cast<long>(frac_str.size());

    if (neg) {
        if (int_part != 0) {
            return EncodedBin{-int_part, frac_part, frac_digits};
        }
        if (frac_part != 0) {
            return EncodedBin{0, -frac_part, frac_digits};
        }
    }

    return EncodedBin{int_part, frac_part, frac_digits};
}

static inline double decode_bin_gev(long int_part, long frac_part, long frac_digits) {
    if (frac_digits < 0) {
        throw std::runtime_error("decode_bin_gev: negative frac_digits");
    }

    const bool neg = (int_part < 0) || (frac_part < 0);
    const long abs_int  = std::labs(int_part);
    const long abs_frac = std::labs(frac_part);

    long double out = static_cast<long double>(abs_int);

    if (frac_digits > 0) {
        const long scale = pow10_long(frac_digits);
        out += static_cast<long double>(abs_frac) / static_cast<long double>(scale);
    }

    if (neg) {
        out = -out;
    }

    return norm_zero(static_cast<double>(out));
}

struct BinnedObservableId {
    ObservableId s;
    std::pair<double, double> p;

    BinnedObservableId() = default;
    BinnedObservableId(ObservableId id) : s(id), p({0.0, 0.0}) {}
    BinnedObservableId(ObservableId id, std::pair<double, double> bin) : s(id), p(bin) {}

    bool operator==(BinnedObservableId const& o) const noexcept {
        return s == o.s
            && norm_zero(p.first)  == norm_zero(o.p.first)
            && norm_zero(p.second) == norm_zero(o.p.second);
    }

    bool operator<(BinnedObservableId const& o) const noexcept {
        const std::uint64_t a = bits_norm_zero(p.first);
        const std::uint64_t b = bits_norm_zero(p.second);
        const std::uint64_t c = bits_norm_zero(o.p.first);
        const std::uint64_t d = bits_norm_zero(o.p.second);

        return std::tie(s, a, b) < std::tie(o.s, c, d);
    }

    LhaID flha() const {
        auto flha_opt = ObservableMapper::flha_of(this->s);
        if (!flha_opt.has_value()) {
            throw std::runtime_error("ObservableId to flha mapping unknown");
        }

        LhaID unbinned_id = flha_opt.value();

        const EncodedBin low  = encode_bin_gev(this->p.first);
        const EncodedBin high = encode_bin_gev(this->p.second);

        std::vector<long> bin_parts{
            low.int_part,
            low.frac_part,
            low.frac_digits,
            high.int_part,
            high.frac_part,
            high.frac_digits
        };

        auto parts = unbinned_id.get_parts();
        parts.insert(parts.begin() + 2, bin_parts.begin(), bin_parts.end());

        return LhaID(parts);
    }

    static BinnedObservableId from_flha(LhaID const& id) {
        auto parts = id.get_parts();

        if (parts.size() < 8) {
            throw std::runtime_error("from_flha: LhaID has not enough parts to contain robust binning");
        }

        const long low_int      = parts.at(2);
        const long low_frac     = parts.at(3);
        const long low_ndigits  = parts.at(4);
        const long high_int     = parts.at(5);
        const long high_frac    = parts.at(6);
        const long high_ndigits = parts.at(7);

        parts.erase(parts.begin() + 2, parts.begin() + 8);
        LhaID unbinned_id(parts);

        LOG_INFO(unbinned_id);

        auto obs_opt = ObservableMapper::from_flha(unbinned_id);
        if (!obs_opt.has_value()) {
            throw std::runtime_error("from_flha: flha to ObservableId mapping unknown");
        }

        BinnedObservableId out;
        out.s = obs_opt.value();
        out.p = {
            decode_bin_gev(low_int, low_frac, low_ndigits),
            decode_bin_gev(high_int, high_frac, high_ndigits)
        };

        return out;
    }

    std::string str() const {
        std::stringstream ss;
        ss << s.str() << " [" << p.first << ", " << p.second << "]";
        return ss.str();
    }
};

template<>
struct std::hash<BinnedObservableId> {
    std::size_t operator()(BinnedObservableId const& id) const noexcept {
        std::size_t h = std::hash<ObservableId>{}(id.s);

        auto mix = [](std::size_t& seed, std::size_t v) {
            seed ^= v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        };

        mix(h, std::hash<std::uint64_t>{}(bits_norm_zero(id.p.first)));
        mix(h, std::hash<std::uint64_t>{}(bits_norm_zero(id.p.second)));

        return h;
    }
};

#endif // BINNEDOBSERVABLEID_H