#if !defined(HYPERISO_WILSON_SUSY_H)
#define HYPERISO_WILSON_SUSY_H

#include "Parameters.h"
#include "Logger.h"
#include "Wilson.h"
#include "epsilon_calculator.h"
#include "Math_SUSY.h"
#include "Math_THDM.h"
#include "Wilson_THDM.h"
#include "Math.h"
#include <algorithm>

class SUSY_LO_Strategy : public SM_LO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
    void init_prime(double scale_W,double scale,int gen, WilsonSet& C_match);
    void init_scalar(double scale_W,double scale,int gen, WilsonSet& C_match);

    
};

class SUSY_NLO_Strategy : public SUSY_LO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
    void init_scalar(double Q_match,double Q,int gen, WilsonSet& C_match);
};

class SUSY_NNLO_Strategy : public SUSY_NLO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
};

#endif // HYPERISO_WILSON_SUSY_H