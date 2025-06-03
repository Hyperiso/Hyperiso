#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <cmath>
#include "common.h"

namespace c2_thdm {

void setMu(const double mu)
{
    setmudim(mu * mu);
}

void setLambda2(const double lambda2)
{
    setlambda(lambda2);
}

void setUVDiv(const double x)
{
    setuvdiv(x);
}

} // End of namespace c2_thdm

