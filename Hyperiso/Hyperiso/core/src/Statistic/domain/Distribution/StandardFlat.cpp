#include "StandardFlat.h"

StandardFlat::StandardFlat(unsigned int seed)
    : eng_(seed),
      dist_(-std::sqrt(3.0), std::sqrt(3.0)) // E=0, Var=1
{}

Vector StandardFlat::sample(std::size_t n) {
    Vector z(n);
    for (std::size_t i = 0; i < n; ++i) z[i] = dist_(eng_);
    return z;
}
