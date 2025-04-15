#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include "common.h"

namespace c5_sm {

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

} // End of namespace c5_sm

