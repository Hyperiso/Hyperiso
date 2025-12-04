#pragma once
#include <utility>
#include <cstddef>

struct PairHash {
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        std::size_t h = static_cast<std::size_t>(p.first);
        h ^= (static_cast<std::size_t>(p.second) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
        return h;
    }
};
