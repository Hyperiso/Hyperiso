#ifndef CHARGEDCURRENTWILSON_SUSY_H
#define CHARGEDCURRENTWILSON_SUSY_H

#include "Wilson.h"

class C_V1_bc_SUSY : public WilsonCoefficient {
public:
    C_V1_bc_SUSY() : WilsonCoefficient("C_V1_bc_SUSY", GroupMapper::str(WGroup::CC_bc) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_bc_SUSY>(*this);
    }

};

class C_V2_bc_SUSY : public WilsonCoefficient {
public:
    C_V2_bc_SUSY() : WilsonCoefficient("C_V2_bc_SUSY", GroupMapper::str(WGroup::CC_bc) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_bc_SUSY>(*this);
    }

};

class C_S1_bc_SUSY : public WilsonCoefficient {
public:
    C_S1_bc_SUSY();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_bc_SUSY>(*this);
    }

};

class C_S2_bc_SUSY : public WilsonCoefficient {
public:
    C_S2_bc_SUSY();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_bc_SUSY>(*this);
    }

};

class C_T_bc_SUSY : public WilsonCoefficient {
public:
    C_T_bc_SUSY() : WilsonCoefficient("C_T_bc_SUSY", GroupMapper::str(WGroup::CC_bc) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_bc_SUSY>(*this);
    }

};


class C_V1_bu_SUSY : public WilsonCoefficient {
public:
    C_V1_bu_SUSY() : WilsonCoefficient("C_V1_bu_SUSY", GroupMapper::str(WGroup::CC_bu) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_bu_SUSY>(*this);
    }

};

class C_V2_bu_SUSY : public WilsonCoefficient {
public:
    C_V2_bu_SUSY() : WilsonCoefficient("C_V2_bu_SUSY", GroupMapper::str(WGroup::CC_bu) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_bu_SUSY>(*this);
    }

};

class C_S1_bu_SUSY : public WilsonCoefficient {
public:
    C_S1_bu_SUSY();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_bu_SUSY>(*this);
    }

};

class C_S2_bu_SUSY : public WilsonCoefficient {
public:
    C_S2_bu_SUSY();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_bu_SUSY>(*this);
    }

};

class C_T_bu_SUSY : public WilsonCoefficient {
public:
    C_T_bu_SUSY() : WilsonCoefficient("C_T_bu_SUSY", GroupMapper::str(WGroup::CC_bu) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_bu_SUSY>(*this);
    }

};


#endif