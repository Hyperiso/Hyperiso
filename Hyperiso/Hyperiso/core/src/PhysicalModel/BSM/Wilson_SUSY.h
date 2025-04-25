#if !defined(HYPERISO_WILSON_SUSY_H)
#define HYPERISO_WILSON_SUSY_H
#include "Wilson.h"
#include "WilsonGroup.h"
#include "susy_parameters.h"
#include "Math_THDM.h"
#include "Utils.h"

class WilsonCoefficient_susy : public WilsonCoefficient {
protected:
    WilsonCoefficient_susy() = default;
    int gen{3};

public:
    void init(QCDOrder order);
};

// TODO : adapt with storage_block constructor as in THDM
class C1_susy : public WilsonCoefficient_susy {
public:
    C1_susy() : WilsonCoefficient_susy() { this->set_name("C1_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation()  {} 
    void NNLO_calculation();
};

class C2_susy : public WilsonCoefficient_susy {
public:
    C2_susy() : WilsonCoefficient_susy() { this->set_name("C2_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation()  {} 
    void NNLO_calculation() {} 

};

class C3_susy : public WilsonCoefficient_susy {
public:
    C3_susy() : WilsonCoefficient_susy() { this->set_name("C3_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

};

class C4_susy : public WilsonCoefficient_susy {
public:
    C4_susy() : WilsonCoefficient_susy() { this->set_name("C4_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation();
    void NNLO_calculation(); 
};

class C5_susy : public WilsonCoefficient_susy {
public:
    C5_susy() : WilsonCoefficient_susy() { this->set_name("C5_SUSY"); }
    
    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();
};

class C6_susy : public WilsonCoefficient_susy {
public:
    C6_susy() : WilsonCoefficient_susy() { this->set_name("C6_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

};

class C7_susy : public WilsonCoefficient_susy {
public:
    C7_susy() : WilsonCoefficient_susy() { this->set_name("C7_SUSY"); }

    void LO_calculation() override;
    void NLO_calculation() override;
    void NNLO_calculation() override;

};

class C8_susy : public WilsonCoefficient_susy {
public:
    C8_susy() : WilsonCoefficient_susy() { this->set_name("C8_SUSY"); }

    void LO_calculation()override;
    void NLO_calculation() override;
    void NNLO_calculation() override; 

};

class C9_susy : public WilsonCoefficient_susy {
public:
    C9_susy() : WilsonCoefficient_susy() { this->set_name("C9_SUSY"); }

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation() {} 

};

class C10_susy : public WilsonCoefficient_susy {
public:
    C10_susy() : WilsonCoefficient_susy() { this->set_name("C10_SUSY"); }

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation() {} 

};

class CQ1_susy : public WilsonCoefficient_susy {
public:
    CQ1_susy() : WilsonCoefficient_susy() { this->set_name("CQ1_SUSY"); }

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation() {} 
};

class CQ2_susy : public WilsonCoefficient_susy {
public:
    CQ2_susy() : WilsonCoefficient_susy() { this->set_name("CQ2_SUSY"); }

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation() {} 
};

class CP1_susy : public WilsonCoefficient_susy {
public:
    CP1_susy() : WilsonCoefficient_susy() { this->set_name("CP1_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CP2_susy : public WilsonCoefficient_susy {
public:
    CP2_susy() : WilsonCoefficient_susy() { this->set_name("CP2_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 
};

class CP3_susy : public WilsonCoefficient_susy {
public:
    CP3_susy() : WilsonCoefficient_susy() { this->set_name("CP3_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CP4_susy : public WilsonCoefficient_susy {
public:
    CP4_susy() : WilsonCoefficient_susy() { this->set_name("CP4_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CP5_susy : public WilsonCoefficient_susy {
public:
    CP5_susy() : WilsonCoefficient_susy() { this->set_name("CP5_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CP6_susy : public WilsonCoefficient_susy {
public:
    CP6_susy() : WilsonCoefficient_susy() { this->set_name("CP6_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CP7_susy : public WilsonCoefficient_susy {
public:
    CP7_susy() : WilsonCoefficient_susy() { this->set_name("CP7_SUSY"); }

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CP8_susy : public WilsonCoefficient_susy {
public:
    CP8_susy() : WilsonCoefficient_susy() { this->set_name("CP8_SUSY"); }

    void LO_calculation(); 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CP9_susy : public WilsonCoefficient_susy {
public:
    CP9_susy() : WilsonCoefficient_susy() { this->set_name("CP9_SUSY"); }

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CP10_susy : public WilsonCoefficient_susy {
public:
    CP10_susy() : WilsonCoefficient_susy() { this->set_name("CP10_SUSY"); }

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CPQ1_susy : public WilsonCoefficient_susy {
public:
    CPQ1_susy() : WilsonCoefficient_susy() { this->set_name("CPQ1_SUSY"); }

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class CPQ2_susy : public WilsonCoefficient_susy {
public:
    CPQ2_susy() : WilsonCoefficient_susy() { this->set_name("CPQ2_SUSY"); }

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 
};

class C_Blnu_A_SUSY : public WilsonCoefficient_susy {
public:
    C_Blnu_A_SUSY() : WilsonCoefficient_susy() { this->set_name("C_Blnu_A_SUSY"); }

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_Blnu_P_SUSY : public WilsonCoefficient_susy {
public:
    C_Blnu_P_SUSY() : WilsonCoefficient_susy() { this->set_name("C_Blnu_P_SUSY"); }

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_V1_SUSY : public WilsonCoefficient_susy {
public:
    C_V1_SUSY() : WilsonCoefficient_susy() { this->set_name("C_V1_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_V2_SUSY : public WilsonCoefficient_susy {
public:
    C_V2_SUSY() : WilsonCoefficient_susy() { this->set_name("C_V2_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_S1_SUSY : public WilsonCoefficient_susy {
public:
    C_S1_SUSY() : WilsonCoefficient_susy() { this->set_name("C_S1_SUSY"); }

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_S2_SUSY : public WilsonCoefficient_susy {
public:
    C_S2_SUSY() : WilsonCoefficient_susy() { this->set_name("C_S2_SUSY"); }

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_T_SUSY : public WilsonCoefficient_susy {
public:
    C_T_SUSY() : WilsonCoefficient_susy() { this->set_name("C_T_SUSY"); }

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class BCoefficientGroup_susy : public BCoefficientGroup {

public:
    BCoefficientGroup_susy() { this->clear();
        this->insert(std::make_pair("C1", std::make_shared<C1_susy>())); this->insert(std::make_pair("C2", std::make_shared<C2_susy>())); this->insert(std::make_pair("C3", std::make_shared<C3_susy>()));
        this->insert(std::make_pair("C4", std::make_shared<C4_susy>()));  this->insert(std::make_pair("C5", std::make_shared<C5_susy>())); this->insert(std::make_pair("C6", std::make_shared<C6_susy>())); 
        this->insert(std::make_pair("C7", std::make_shared<C7_susy>()));  this->insert(std::make_pair("C8", std::make_shared<C8_susy>()));  this->insert(std::make_pair("C9", std::make_shared<C9_susy>())); 
        this->insert(std::make_pair("C10", std::make_shared<C10_susy>())); 

        this->storage_block = "B_HADRONIC_SUSY";
        this->wilson_type = ContributionType::BSM;
    }
};

class BPrimeCoefficientGroup_susy : public BPrimeCoefficientGroup {
public:
    BPrimeCoefficientGroup_susy() { this->clear();
        this->insert(std::make_pair("CP1", std::make_shared<CP1_susy>())); this->insert(std::make_pair("CP2", std::make_shared<CP2_susy>())); this->insert(std::make_pair("CP3", std::make_shared<CP3_susy>()));
        this->insert(std::make_pair("CP4", std::make_shared<CP4_susy>()));  this->insert(std::make_pair("CP5", std::make_shared<CP5_susy>())); this->insert(std::make_pair("CP6", std::make_shared<CP6_susy>())); 
        this->insert(std::make_pair("CP7", std::make_shared<CP7_susy>()));  this->insert(std::make_pair("CP8", std::make_shared<CP8_susy>()));  this->insert(std::make_pair("CP9", std::make_shared<CP9_susy>())); 
        this->insert(std::make_pair("CP10", std::make_shared<CP10_susy>())); this->insert(std::make_pair("CPQ1", std::make_shared<CPQ1_susy>())); this->insert(std::make_pair("CPQ2", std::make_shared<CPQ2_susy>())); 

        this->storage_block = "B_PRIME_HADRONIC_SUSY";
        this->wilson_type = ContributionType::BSM;
    }
};

class BScalarCoefficientGroup_susy : public BScalarCoefficientGroup {
public:
    BScalarCoefficientGroup_susy() : BScalarCoefficientGroup() { this->clear();
        this->insert(std::make_pair("CQ1", std::make_shared<CQ1_susy>())); this->insert(std::make_pair("CQ2", std::make_shared<CQ2_susy>()));

        this->storage_block = "B_SCALAR_HADRONIC_SUSY";
        this->wilson_type = ContributionType::BSM;
    }

    void set_base_1_LO();

private:
    // std::shared_ptr<Parameters> susy = Parameters::GetInstance(ParameterType::BSM);
    susy_parameters* sus_param;
};

class BlnuCoefficientGroup_SUSY : public BlnuCoefficientGroup {
public:
    BlnuCoefficientGroup_SUSY() { this->clear();
        this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A_SUSY>())); 
        this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P_SUSY>()));

        this->storage_block = "B_LNU_HADRONIC_SUSY";
        this->wilson_type = ContributionType::BSM;
    }
};

class BclnuCoefficientGroup_SUSY : public BclnuCoefficientGroup {
public:
    BclnuCoefficientGroup_SUSY() { this->clear();
        this->insert(std::make_pair("C_V1", std::make_shared<C_V1_SUSY>()));
        this->insert(std::make_pair("C_V2", std::make_shared<C_V2_SUSY>()));
        this->insert(std::make_pair("C_S1", std::make_shared<C_S1_SUSY>()));
        this->insert(std::make_pair("C_S2", std::make_shared<C_S2_SUSY>()));
        this->insert(std::make_pair("C_T", std::make_shared<C_T_SUSY>()));

        this->storage_block = "B_CLNU_HADRONIC_SUSY";
        this->wilson_type = ContributionType::BSM;
    }

    void set_base_1() {}
    void set_base_2() {}
};

#endif