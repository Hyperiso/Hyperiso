#ifndef __STUDENTTCOPULA_H__
#define __STUDENTTCOPULA_H__

#include "GenericCopula.h"
#include "Matrix.h"
#include "AbstractConfig.h"
#include <gsl/gsl_sf.h>

struct StudentTCopulaConfig : public AbstractConfig {
    RealMatrix R;
    int nu;
};

class StudentTCopula : public GenericCopula {
public:
    explicit StudentTCopula(unsigned int seed, RealMatrix R, int nu);

    std::vector<Vector> sample_u(std::size_t n) override;
    double log_density(Vector u) override;

private:
    int nu;
    RealMatrix R, R_inv;
    RealMatrix L;
    double logdet;
};

#endif // __STUDENTTCOPULA_H__
