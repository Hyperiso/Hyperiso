#ifndef __BINNEDOBSERVABLEID_H__
#define __BINNEDOBSERVABLEID_H__

#include <bit>
#include <cstdint>
#include <cmath>
#include <string>
#include <tuple>
#include <functional>
#include "observable_ids.hpp"

static inline double norm_zero(double x) noexcept {
    // Rend -0.0 et +0.0 identiques
    return (x == 0.0) ? 0.0 : x;
}

static inline std::uint64_t bits_norm_zero(double x) noexcept {
    x = norm_zero(x);
    return std::bit_cast<std::uint64_t>(x);
}

struct BinnedObservableId {
    ObservableId s;
    std::pair<double,double> p;

    BinnedObservableId() = default;
    BinnedObservableId(ObservableId id) : s(id), p({0, 0}) {}
    BinnedObservableId(ObservableId id, std::pair<double, double> bin) : s(id), p(bin) {}

    bool operator==(BinnedObservableId const& o) const noexcept {
        return s == o.s
            && norm_zero(p.first)  == norm_zero(o.p.first)
            && norm_zero(p.second) == norm_zero(o.p.second);
    }

    bool operator<(BinnedObservableId const& o) const noexcept {
        std::uint64_t a = bits_norm_zero(p.first);
        std::uint64_t b = bits_norm_zero(p.second);
        std::uint64_t c =  bits_norm_zero(o.p.first);
        std::uint64_t d = bits_norm_zero(o.p.second);
        return std::tie(s,a
                        ,b
                        )
             < std::tie(o.s,c
                       ,d)
                        ;
    }

    LhaID flha() const {
        if (!ObservableMapper::flha_of(this->s).has_value())
            throw std::runtime_error("ObservableId to flha mapping unknown");

        LhaID unbinned_id = ObservableMapper::flha_of(this->s).value();
        std::vector<long> bin {
            static_cast<long>(std::floor(1000 * this->p.first)),
            static_cast<long>(std::floor(1000 * this->p.second))
        };
        auto parts = unbinned_id.get_parts();
        parts.insert(parts.begin() + 2, bin.begin(), bin.end());
        return LhaID(parts);
    }

    static BinnedObservableId from_flha(LhaID const& id) {
        auto parts = id.get_parts();

        if (parts.size() < 4)
            throw std::runtime_error("from_flha: LhaID has not enough parts to contain binning");

        const long bin0 = parts.at(2);
        const long bin1 = parts.at(3);

        parts.erase(parts.begin() + 2, parts.begin() + 4);
        LhaID unbinned_id(parts);

        LOG_INFO(unbinned_id);
        auto obs_opt = ObservableMapper::from_flha(unbinned_id);
        if (!obs_opt.has_value())
            throw std::runtime_error("from_flha: flha to ObservableId mapping unknown");

        BinnedObservableId out;
        out.s = obs_opt.value();
        out.p = { static_cast<double>(bin0) / 1000.0,
                  static_cast<double>(bin1) / 1000.0 };
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
            seed ^= v + 0x9e3779b97f4a7c15ULL + (seed<<6) + (seed>>2);
        };

        mix(h, std::hash<std::uint64_t>{}(bits_norm_zero(id.p.first)));
        mix(h, std::hash<std::uint64_t>{}(bits_norm_zero(id.p.second)));

        return h;
    }
};

#endif // __BINNEDOBSERVABLEID_H__
