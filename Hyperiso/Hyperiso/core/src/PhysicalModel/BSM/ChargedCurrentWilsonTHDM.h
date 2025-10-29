#ifndef CHARGEDCURRENTWILSON_THDM_H
#define CHARGEDCURRENTWILSON_THDM_H

#include "Wilson.h"

class C_V1_bc_THDM : public WilsonCoefficient {
public:
    C_V1_bc_THDM() : WilsonCoefficient("C_V1_bc_THDM", GroupMapper::str(WGroup::BCC_bc) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_bc_THDM>(*this);
    }

};

class C_V2_bc_THDM : public WilsonCoefficient {
public:
    C_V2_bc_THDM() : WilsonCoefficient("C_V2_bc_THDM", GroupMapper::str(WGroup::BCC_bc) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_bc_THDM>(*this);
    }

};

class C_S1_bc_THDM : public WilsonCoefficient {
public:
    C_S1_bc_THDM();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_bc_THDM>(*this);
    }

};

class C_S2_bc_THDM : public WilsonCoefficient {
public:
    C_S2_bc_THDM();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_bc_THDM>(*this);
    }

};

class C_T_bc_THDM : public WilsonCoefficient {
public:
    C_T_bc_THDM() : WilsonCoefficient("C_T_bc_THDM", GroupMapper::str(WGroup::BCC_bc) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_bc_THDM>(*this);
    }

};



class C_V1_bu_THDM : public WilsonCoefficient {
public:
    C_V1_bu_THDM() : WilsonCoefficient("C_V1_bu_THDM", GroupMapper::str(WGroup::BCC_bu) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_bu_THDM>(*this);
    }

};

class C_V2_bu_THDM : public WilsonCoefficient {
public:
    C_V2_bu_THDM() : WilsonCoefficient("C_V2_bu_THDM", GroupMapper::str(WGroup::BCC_bu) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_bu_THDM>(*this);
    }

};

class C_S1_bu_THDM : public WilsonCoefficient {
public:
    C_S1_bu_THDM();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_bu_THDM>(*this);
    }

};

class C_S2_bu_THDM : public WilsonCoefficient {
public:
    C_S2_bu_THDM();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_bu_THDM>(*this);
    }

};

class C_T_bu_THDM : public WilsonCoefficient {
public:
    C_T_bu_THDM() : WilsonCoefficient("C_T_bu_THDM", GroupMapper::str(WGroup::BCC_bu) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_bu_THDM>(*this);
    }

};

#endif