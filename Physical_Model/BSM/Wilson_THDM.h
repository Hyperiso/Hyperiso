#include "../../Core/Parameters.h"
#include "../../Core/Logger.h"
#include "../Wilson.h"
#include "../../Math/Math.h"

class THDM_LO_Strategy : public InitializationStrategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override {}

};


class THDM_NLO_Strategy : public InitializationStrategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override {}

};


class THDM_NNLO_Strategy : public InitializationStrategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override {}

};