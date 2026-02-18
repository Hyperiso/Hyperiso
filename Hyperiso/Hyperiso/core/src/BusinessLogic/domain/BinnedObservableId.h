#ifndef __BINNEDOBSERVABLEID_H__
#define __BINNEDOBSERVABLEID_H__

#include <bit>
#include <cstdint>
#include <cmath>
#include <string>
#include <tuple>
#include <functional>
#include "Include.h"

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

    bool operator==(BinnedObservableId const& o) const noexcept {
        return s == o.s
            && norm_zero(p.first)  == norm_zero(o.p.first)
            && norm_zero(p.second) == norm_zero(o.p.second);
    }

    bool operator<(BinnedObservableId const& o) const noexcept {
        return std::tie(s,
                        bits_norm_zero(p.first),
                        bits_norm_zero(p.second))
             < std::tie(o.s,
                        bits_norm_zero(o.p.first),
                        bits_norm_zero(o.p.second));
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
