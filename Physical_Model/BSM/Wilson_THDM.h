#include "Parameters.h"
#include "Logger.h"
#include "Wilson.h"
#include "Math.h"
#include "Math_THDM.h"

class THDM_LO_Strategy : public SM_LO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
    void init_prime(double scale_W,double scale,int gen, WilsonSet& C_match) override {}
    void init_scalar(double scale_W,double scale,int gen, WilsonSet& C_match) override;

    void set_lu(double lu) {this->lu = lu;
    Logger::getInstance()->info("lu in THDM " + std::to_string(lu));}
    void set_ld(double ld) {this->ld = ld;
    Logger::getInstance()->info("ld in THDM " + std::to_string(ld));}

protected:
    double lu{-1.};
    double ld{-1.};
};


class THDM_NLO_Strategy : public THDM_LO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
};


class THDM_NNLO_Strategy : public THDM_NLO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
};