#include "../../Core/Parameters.h"
#include "../../Core/Logger.h"
#include "../Wilson.h"
#include "epsilon_calculator.h"
#include "../Math/Math_BSM/Math_SUSY.h"
#include "Math.h"
#include "Math_SUSY.h"
#include <algorithm>

class SUSY_LO_Strategy : public InitializationStrategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match) override;
    void init_prime(double scale_W,double scale,int gen, WilsonSet& C_match) override;
    void init_scalar(double scale_W,double scale,int gen, WilsonSet& C_match) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override {}

};

class SUSY_NLO_Strategy : public SUSY_LO_Strategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override {}

};

class SUSY_NNLO_Strategy : public SUSY_NLO_Strategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override {}

};