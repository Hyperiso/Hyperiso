#ifndef HASH_UTILS_H
#define HASH_UTILS_H

#include <utility>
#include <cstddef>

/**
 * @file PairHash.h
 * @brief Hash functor for std::pair<int,int>.
 *
 * This header defines a small utility struct PairHash which can be used
 * as a custom hash functor for std::pair<int,int>, making it convenient
 * to use pairs of integers as keys in unordered containers such as
 * std::unordered_map and std::unordered_set.
 */

/**
 * @struct PairHash
 * @brief Hash functor for std::pair<int,int>.
 *
 * This functor implements a 64-bit hash combine pattern inspired by
 * Boost's hash_combine. The two integer components are mixed into a single
 * std::size_t value, providing a reasonably well-distributed hash for
 * use in hash-based containers.
 *
 * Typical usage:
 * @code
 * #include <unordered_map>
 *
 * std::unordered_map<std::pair<int,int>, double, PairHash> cache;
 * cache[{1,2}] = 3.14;
 * @endcode
 */
struct PairHash {
    /**
     * @brief Computes a hash value for a pair of integers.
     *
     * The hash is built from the first component, then combined with
     * the second using a Boost-style mixing step:
     *   h ^= (second + constant + (h << 6) + (h >> 2));
     *
     * @param p Pair of integers to hash.
     * @return Combined hash value.
     */
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        std::size_t h = static_cast<std::size_t>(p.first);
        h ^= (static_cast<std::size_t>(p.second) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
        return h;
    }
};

#endif