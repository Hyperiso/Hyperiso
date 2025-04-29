#ifndef __CHARGEDCURRENTWILSON_H__
#define __CHARGEDCURRENTWILSON_H__

#include "Wilson.h"

class C_Blnu_A : public WilsonCoefficient {
public:
    C_Blnu_A() : WilsonCoefficient("C_Blnu_A", GroupMapper::str(WGroup::Blnu) + "_MATCH") {}

    void LO_calculation() override;
    void NLO_calculation() override {} 
    void NNLO_calculation() override {} 
};

class C_Blnu_P : public WilsonCoefficient {
public:
    C_Blnu_P() : WilsonCoefficient("C_Blnu_P", GroupMapper::str(WGroup::Blnu) + "_MATCH") {}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_V1 : public WilsonCoefficient {
public:
    C_V1() : WilsonCoefficient("C_V1", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_V2 : public WilsonCoefficient {
public:
    C_V2() : WilsonCoefficient("C_V2", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_S1 : public WilsonCoefficient {
public:
    C_S1() : WilsonCoefficient("C_S1", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_S2 : public WilsonCoefficient {
public:
    C_S2() : WilsonCoefficient("C_S2", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_T : public WilsonCoefficient {
public:
    C_T() : WilsonCoefficient("C_T", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 
};

#endif // __CHARGEDCURRENTWILSON_H__
