#include "StandardNormal.h"

StandardNormal::StandardNormal(unsigned int seed)
    : eng_(seed), dist_(0.0, 1.0) {}

Vector StandardNormal::sample(std::size_t n) {
    Vector z(n);
    for (std::size_t i = 0; i < n; ++i) z[i] = dist_(eng_);
    return z;
}