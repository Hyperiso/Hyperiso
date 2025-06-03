#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <cmath>
#include "marty/core/looptools_interface.h"
#include "common.h"

namespace c8_thdm {

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

} // End of namespace c8_thdm

