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
    LOG_INFO("lu in THDM " + std::to_string(lu)); is_thdm=false;}
    void set_ld(double ld) {this->ld = ld;
    LOG_INFO("ld in THDM " + std::to_string(ld)); is_thdm=false;}
    void set_le(double le) {this->le = le;
    LOG_INFO("le in THDM " + std::to_string(le)); is_thdm=false;}

protected:
    double lu{-1.};
    double ld{-1.};
    double le{-1};
    bool is_thdm{true};
};


class THDM_NLO_Strategy : public THDM_LO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
};


class THDM_NNLO_Strategy : public THDM_NLO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
};