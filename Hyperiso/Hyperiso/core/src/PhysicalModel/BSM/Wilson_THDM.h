#if !defined(HYPERISO_WILSON_THDM_H)
#define HYPERISO_WILSON_THDM_H
#include "Wilson.h"
#include "WilsonGroup.h"
#include "thdm_parameters.h"
#include "Math_THDM.h"
#include "Utils.h"

class WilsonCoefficient_THDM : public WilsonCoefficient {
protected:
    WilsonCoefficient_THDM() = default;

public:
    void init(QCDOrder order);
};

class C1_THDM : public WilsonCoefficient_THDM {
public:
    C1_THDM() : WilsonCoefficient_THDM() { this->set_name("C1_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 
};

class C2_THDM : public WilsonCoefficient_THDM {
public:
    C2_THDM() : WilsonCoefficient_THDM() { this->set_name("C2_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 
};

class C3_THDM : public WilsonCoefficient_THDM {
public:
    C3_THDM() : WilsonCoefficient_THDM() { this->set_name("C3_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation()  {} 
    void NNLO_calculation();

};

class C4_THDM : public WilsonCoefficient_THDM {
public:
    C4_THDM() : WilsonCoefficient_THDM() { this->set_name("C4_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation();
    void NNLO_calculation();

};

class C5_THDM : public WilsonCoefficient_THDM {
public:
    C5_THDM() :  WilsonCoefficient_THDM() { this->set_name("C5_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

};

class C6_THDM : public WilsonCoefficient_THDM {
public:
    C6_THDM() : WilsonCoefficient_THDM() { this->set_name("C6_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

};

class C7_THDM : public WilsonCoefficient_THDM {
public:
    C7_THDM() : WilsonCoefficient_THDM() { this->set_name("C7_THDM"); }

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();
};

class C8_THDM : public WilsonCoefficient_THDM {
public:
    C8_THDM() : WilsonCoefficient_THDM() { this->set_name("C8_THDM"); }

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();
};

class C9_THDM : public WilsonCoefficient_THDM {
public:
    C9_THDM() : WilsonCoefficient_THDM() { this->set_name("C9_THDM"); }

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation() {} 
};

class C10_THDM : public WilsonCoefficient_THDM {
public:
    C10_THDM() : WilsonCoefficient_THDM() { this->set_name("C10_THDM"); }

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation() {} 
};

class CQ1_THDM : public WilsonCoefficient_THDM {
public:
    CQ1_THDM() : WilsonCoefficient_THDM() { this->set_name("CQ1_THDM"); }

    void LO_calculation();
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CQ2_THDM : public WilsonCoefficient_THDM {
public:
    CQ2_THDM() : WilsonCoefficient_THDM() { this->set_name("CQ2_THDM"); }

    void LO_calculation();
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP1_THDM : public WilsonCoefficient_THDM {
public:
    CP1_THDM() : WilsonCoefficient_THDM() { this->set_name("CP1_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP2_THDM : public WilsonCoefficient_THDM {
public:
    CP2_THDM() : WilsonCoefficient_THDM() { this->set_name("CP2_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP3_THDM : public WilsonCoefficient_THDM {
public:
    CP3_THDM() : WilsonCoefficient_THDM() { this->set_name("CP3_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP4_THDM : public WilsonCoefficient_THDM {
public:
    CP4_THDM() : WilsonCoefficient_THDM() { this->set_name("CP4_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP5_THDM : public WilsonCoefficient_THDM {
public:
    CP5_THDM() : WilsonCoefficient_THDM() { this->set_name("CP5_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP6_THDM : public WilsonCoefficient_THDM {
public:
    CP6_THDM() : WilsonCoefficient_THDM() { this->set_name("CP6_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP7_THDM : public WilsonCoefficient_THDM {
public:
    CP7_THDM() : WilsonCoefficient_THDM() { this->set_name("CP7_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP8_THDM : public WilsonCoefficient_THDM {
public:
    CP8_THDM() : WilsonCoefficient_THDM() { this->set_name("CP8_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP9_THDM : public WilsonCoefficient_THDM {
public:
    CP9_THDM() : WilsonCoefficient_THDM() { this->set_name("CP9_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CP10_THDM : public WilsonCoefficient_THDM {
public:
    CP10_THDM() : WilsonCoefficient_THDM() { this->set_name("CP10_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CPQ1_THDM : public WilsonCoefficient_THDM {
public:
    CPQ1_THDM() : WilsonCoefficient_THDM() { this->set_name("CPQ1_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class CPQ2_THDM : public WilsonCoefficient_THDM {
public:
    CPQ2_THDM() : WilsonCoefficient_THDM() { this->set_name("CPQ2_THDM"); }

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {} 
};

class C_Blnu_A_THDM : public WilsonCoefficient_THDM {
public:
    C_Blnu_A_THDM() : WilsonCoefficient_THDM() { this->set_name("C_Blnu_A_THDM"); }

    void LO_calculation() { }
    void NLO_calculation() { } 
    void NNLO_calculation() { } 

};

class C_Blnu_P_THDM : public WilsonCoefficient_THDM {
public:
    C_Blnu_P_THDM() : WilsonCoefficient_THDM() { this->set_name("C_Blnu_P_THDM"); }

    void LO_calculation();
    void NLO_calculation() { } 
    void NNLO_calculation() { } 

};

class C_V1_THDM : public WilsonCoefficient_THDM {
public:
    C_V1_THDM() : WilsonCoefficient_THDM() { this->set_name("C_V1_THDM"); }

    void LO_calculation() { } 
    void NLO_calculation() { } 
    void NNLO_calculation() { } 

};

class C_V2_THDM : public WilsonCoefficient_THDM {
public:
    C_V2_THDM() : WilsonCoefficient_THDM() { this->set_name("C_V2_THDM"); }

    void LO_calculation() { } 
    void NLO_calculation() { } 
    void NNLO_calculation() { } 

};

class C_S1_THDM : public WilsonCoefficient_THDM {
public:
    C_S1_THDM() : WilsonCoefficient_THDM() { this->set_name("C_S1_THDM"); }

    void LO_calculation();
    void NLO_calculation() { } 
    void NNLO_calculation() { } 

};

class C_S2_THDM : public WilsonCoefficient_THDM {
public:
    C_S2_THDM() : WilsonCoefficient_THDM() { this->set_name("C_S2_THDM"); }

    void LO_calculation();
    void NLO_calculation() { } 
    void NNLO_calculation() { } 

};

class C_T_THDM : public WilsonCoefficient_THDM {
public:
    C_T_THDM() : WilsonCoefficient_THDM() { this->set_name("C_T_THDM"); }

    void LO_calculation() { } 
    void NLO_calculation() { } 
    void NNLO_calculation() { } 

};

class BCoefficientGroup_THDM : public BCoefficientGroup {

public:
    BCoefficientGroup_THDM() {this->clear();
        this->insert(std::make_pair("C1", std::make_shared<C1_THDM>())); this->insert(std::make_pair("C2", std::make_shared<C2_THDM>())); this->insert(std::make_pair("C3", std::make_shared<C3_THDM>()));
        this->insert(std::make_pair("C4", std::make_shared<C4_THDM>()));  this->insert(std::make_pair("C5", std::make_shared<C5_THDM>())); this->insert(std::make_pair("C6", std::make_shared<C6_THDM>())); 
        this->insert(std::make_pair("C7", std::make_shared<C7_THDM>()));  this->insert(std::make_pair("C8", std::make_shared<C8_THDM>()));  this->insert(std::make_pair("C9", std::make_shared<C9_THDM>())); 
        this->insert(std::make_pair("C10", std::make_shared<C10_THDM>())); 

        this->storage_block = "B_HADRONIC_THDM";
    }
};

class BPrimeCoefficientGroup_THDM : public BPrimeCoefficientGroup {
public:
    BPrimeCoefficientGroup_THDM() { this->clear();
        this->insert(std::make_pair("CP1", std::make_shared<CP1_THDM>())); this->insert(std::make_pair("CP2", std::make_shared<CP2_THDM>())); this->insert(std::make_pair("CP3", std::make_shared<CP3_THDM>()));
        this->insert(std::make_pair("CP4", std::make_shared<CP4_THDM>()));  this->insert(std::make_pair("CP5", std::make_shared<CP5_THDM>())); this->insert(std::make_pair("CP6", std::make_shared<CP6_THDM>())); 
        this->insert(std::make_pair("CP7", std::make_shared<CP7_THDM>()));  this->insert(std::make_pair("CP8", std::make_shared<CP8_THDM>()));  this->insert(std::make_pair("CP9", std::make_shared<CP9_THDM>())); 
        this->insert(std::make_pair("CP10", std::make_shared<CP10_THDM>())); this->insert(std::make_pair("CPQ1", std::make_shared<CPQ1_THDM>())); this->insert(std::make_pair("CPQ2", std::make_shared<CPQ2_THDM>())); 

        this->storage_block = "B_PRIME_HADRONIC_THDM";
    }

    void set_base_1();
    void set_base_2();
};

class BScalarCoefficientGroup_THDM : public BScalarCoefficientGroup {
public:
    BScalarCoefficientGroup_THDM() { this->clear();
        this->insert(std::make_pair("CQ1", std::make_shared<CQ1_THDM>())); this->insert(std::make_pair("CQ2", std::make_shared<CQ2_THDM>()));

        this->storage_block = "B_SCALAR_HADRONIC_THDM";
    }

    void set_base_1();
    void set_base_2();
};

class BlnuCoefficientGroup_THDM : public BlnuCoefficientGroup {
public:
    BlnuCoefficientGroup_THDM() { this->clear();
        this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A_THDM>())); 
        this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P_THDM>()));

        this->storage_block = "B_LNU_HADRONIC_THDM";
    }

    void set_base_1() {}
    void set_base_2() {}
};

class BclnuCoefficientGroup_THDM : public BclnuCoefficientGroup {
public:
    BclnuCoefficientGroup_THDM() { this->clear();
        this->insert(std::make_pair("C_V1", std::make_shared<C_V1_THDM>()));
        this->insert(std::make_pair("C_V2", std::make_shared<C_V2_THDM>()));
        this->insert(std::make_pair("C_S1", std::make_shared<C_S1_THDM>()));
        this->insert(std::make_pair("C_S2", std::make_shared<C_S2_THDM>()));
        this->insert(std::make_pair("C_T", std::make_shared<C_T_THDM>()));

        this->storage_block = "B_CLNU_HADRONIC_THDM";
    }

    void set_base_1() {}
    void set_base_2() {}
};

#endif