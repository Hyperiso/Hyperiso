#ifndef INTERPRETED_PARAM_H
#define INTERPRETED_PARAM_H

#include "Include.h"

struct InterpretedParam {
    std::string block;
    LhaID code;
    bool is_bsm;
    bool is_complex;

    bool operator==(const InterpretedParam& other) const {
        return block == other.block &&
                code == other.code &&
                is_bsm == other.is_bsm &&
                is_complex == other.is_complex;
    }
};

inline void hash_combine(std::size_t& seed, std::size_t value) noexcept {
    seed ^= value + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
}

namespace std {
    template<>
    struct hash<InterpretedParam> {
        std::size_t operator()(const InterpretedParam& p) const noexcept {
            std::size_t seed = 0;
            hash_combine(seed, std::hash<std::string>{}(p.block));
            hash_combine(seed, hash<LhaID>{}(p.code));
            hash_combine(seed, std::hash<bool>{}(p.is_bsm));
            hash_combine(seed, std::hash<bool>{}(p.is_complex));
            return seed;
        }
    };
}


#endif