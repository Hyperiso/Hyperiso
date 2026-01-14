#ifndef __ICOPULA_H__
#define __ICOPULA_H__

#include "RNGHelper.h"

class ICopula {
public:
    virtual ~ICopula() = default;

    virtual std::vector<Vector> sample_u(std::size_t n) = 0;
    virtual Vector sample_u() = 0;
    virtual double log_density(Vector u) = 0;
};

#endif // __ICOPULA_H__
