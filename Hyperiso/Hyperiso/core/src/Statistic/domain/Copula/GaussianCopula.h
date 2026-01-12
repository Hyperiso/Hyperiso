#ifndef __GAUSSIANCOPULA_H__
#define __GAUSSIANCOPULA_H__

#include "GenericCopula.h"
#include "Matrix.h"

struct GaussianCopulaConfig : public AbstractConfig {
    RealMatrix R;
};

class GaussianCopula final : public GenericCopula {
public:
    explicit GaussianCopula(unsigned int seed, RealMatrix R);

    std::vector<Vector> sample_u(std::size_t n) override;
    double log_density(Vector u) override;

private:
    RealMatrix R, R_inv;
    RealMatrix L;
    double logdet;
};

#endif // __GAUSSIANCOPULA_H__
