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
    WilsonCoefficient_THDM(const std::string& name, const std::string& storage_block) : WilsonCoefficient(name, storage_block) {};

public:
    void init(QCDOrder order);

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {}
    
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<WilsonCoefficient_THDM>(*this);
    }
};

class C1_THDM : public WilsonCoefficient_THDM {
public:
    C1_THDM() : WilsonCoefficient_THDM("C1_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C1_THDM>(*this);
    }
};

class C2_THDM : public WilsonCoefficient_THDM {
public:
    C2_THDM() : WilsonCoefficient_THDM("C2_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C2_THDM>(*this);
    }
};

class C3_THDM : public WilsonCoefficient_THDM {
public:
    C3_THDM() : WilsonCoefficient_THDM("C3_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation()  {} 
    void NNLO_calculation();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C3_THDM>(*this);
    }
};

class C4_THDM : public WilsonCoefficient_THDM {
public:
    C4_THDM() : WilsonCoefficient_THDM("C4_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation();
    void NNLO_calculation();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C4_THDM>(*this);
    }
};

class C5_THDM : public WilsonCoefficient_THDM {
public:
    C5_THDM() : WilsonCoefficient_THDM("C5_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C5_THDM>(*this);
    }
};

class C6_THDM : public WilsonCoefficient_THDM {
public:
    C6_THDM() : WilsonCoefficient_THDM("C6_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C6_THDM>(*this);
    }
};

class C7_THDM : public WilsonCoefficient_THDM {
public:
    C7_THDM() : WilsonCoefficient_THDM("C7_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C7_THDM>(*this);
    }
};

class C8_THDM : public WilsonCoefficient_THDM {
public:
    C8_THDM() : WilsonCoefficient_THDM("C8_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C8_THDM>(*this);
    }
};

class C9_THDM : public WilsonCoefficient_THDM {
public:
    C9_THDM() : WilsonCoefficient_THDM("C9_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C9_THDM>(*this);
    }
};

class C10_THDM : public WilsonCoefficient_THDM {
public:
    C10_THDM() : WilsonCoefficient_THDM("C10_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C10_THDM>(*this);
    }
};

class CQ1_THDM : public WilsonCoefficient_THDM {
public:
    CQ1_THDM() : WilsonCoefficient_THDM("CQ1_THDM", GroupMapper::str(WGroup::BScalar) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ1_THDM>(*this);
    }
};

class CQ2_THDM : public WilsonCoefficient_THDM {
public:
    CQ2_THDM() : WilsonCoefficient_THDM("CQ2_THDM", GroupMapper::str(WGroup::BScalar) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ2_THDM>(*this);
    }
};

class CP1_THDM : public WilsonCoefficient_THDM {
public:
    CP1_THDM() : WilsonCoefficient_THDM("CP1_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP1_THDM>(*this);
    }
};

class CP2_THDM : public WilsonCoefficient_THDM {
public:
    CP2_THDM() : WilsonCoefficient_THDM("CP2_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP2_THDM>(*this);
    }
};

class CP3_THDM : public WilsonCoefficient_THDM {
public:
    CP3_THDM() : WilsonCoefficient_THDM("CP3_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP3_THDM>(*this);
    }
};

class CP4_THDM : public WilsonCoefficient_THDM {
public:
    CP4_THDM() : WilsonCoefficient_THDM("CP4_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP4_THDM>(*this);
    }
};

class CP5_THDM : public WilsonCoefficient_THDM {
public:
    CP5_THDM() : WilsonCoefficient_THDM("CP5_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP5_THDM>(*this);
    }
};

class CP6_THDM : public WilsonCoefficient_THDM {
public:
    CP6_THDM() : WilsonCoefficient_THDM("CP6_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP6_THDM>(*this);
    }
};

class CP7_THDM : public WilsonCoefficient_THDM {
public:
    CP7_THDM() : WilsonCoefficient_THDM("CP7_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP7_THDM>(*this);
    }
};

class CP8_THDM : public WilsonCoefficient_THDM {
public:
    CP8_THDM() : WilsonCoefficient_THDM("CP8_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP8_THDM>(*this);
    }
};

class CP9_THDM : public WilsonCoefficient_THDM {
public:
    CP9_THDM() : WilsonCoefficient_THDM("CP9_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP9_THDM>(*this);
    }
};

class CP10_THDM : public WilsonCoefficient_THDM {
public:
    CP10_THDM() : WilsonCoefficient_THDM("CP10_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP10_THDM>(*this);
    }
};

class CPQ1_THDM : public WilsonCoefficient_THDM {
public:
    CPQ1_THDM() : WilsonCoefficient_THDM("CPQ1_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ1_THDM>(*this);
    }
};

class CPQ2_THDM : public WilsonCoefficient_THDM {
public:
    CPQ2_THDM() : WilsonCoefficient_THDM("CPQ2_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {}
    void NNLO_calculation() {}

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ2_THDM>(*this);
    }
};

class C_Blnu_A_THDM : public WilsonCoefficient_THDM {
public:
    C_Blnu_A_THDM() : WilsonCoefficient_THDM("C_Blnu_A_THDM", GroupMapper::str(WGroup::Blnu) + "_MATCH") {}

    void LO_calculation() { }
    void NLO_calculation() { } 
    void NNLO_calculation() { }
    
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_Blnu_A_THDM>(*this);
    }

};

class C_Blnu_P_THDM : public WilsonCoefficient_THDM {
public:
    C_Blnu_P_THDM() : WilsonCoefficient_THDM("C_Blnu_P_THDM", GroupMapper::str(WGroup::Blnu) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation() { } 
    void NNLO_calculation() { }

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_Blnu_P_THDM>(*this);
    }

};

class C_V1_THDM : public WilsonCoefficient_THDM {
public:
    C_V1_THDM() : WilsonCoefficient_THDM("C_V1_THDM", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation() { } 
    void NLO_calculation() { } 
    void NNLO_calculation() { }

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_THDM>(*this);
    }

};

class C_V2_THDM : public WilsonCoefficient_THDM {
public:
    C_V2_THDM() : WilsonCoefficient_THDM("C_V2_THDM", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation() { } 
    void NLO_calculation() { } 
    void NNLO_calculation() { }

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_THDM>(*this);
    }

};

class C_S1_THDM : public WilsonCoefficient_THDM {
public:
    C_S1_THDM() : WilsonCoefficient_THDM("C_S1_THDM", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation() { } 
    void NNLO_calculation() { }

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_THDM>(*this);
    }

};

class C_S2_THDM : public WilsonCoefficient_THDM {
public:
    C_S2_THDM() : WilsonCoefficient_THDM("C_S2_THDM", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation();
    void NLO_calculation() { } 
    void NNLO_calculation() { }

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_THDM>(*this);
    }

};

class C_T_THDM : public WilsonCoefficient_THDM {
public:
    C_T_THDM() : WilsonCoefficient_THDM("C_T_THDM", GroupMapper::str(WGroup::BCLNU) + "_MATCH") {}

    void LO_calculation() { } 
    void NLO_calculation() { } 
    void NNLO_calculation() { }

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_THDM>(*this);
    }

};

class BCoefficientGroup_THDM : public BCoefficientGroup {

public:
    BCoefficientGroup_THDM() {this->clear();

        if (!thdm_parameters::is_initialized()) {
            thdm_parameters::init();
        }
        this->insert(std::make_pair("C1", std::make_shared<C1_THDM>())); this->insert(std::make_pair("C2", std::make_shared<C2_THDM>())); this->insert(std::make_pair("C3", std::make_shared<C3_THDM>()));
        this->insert(std::make_pair("C4", std::make_shared<C4_THDM>()));  this->insert(std::make_pair("C5", std::make_shared<C5_THDM>())); this->insert(std::make_pair("C6", std::make_shared<C6_THDM>())); 
        this->insert(std::make_pair("C7", std::make_shared<C7_THDM>()));  this->insert(std::make_pair("C8", std::make_shared<C8_THDM>()));  this->insert(std::make_pair("C9", std::make_shared<C9_THDM>())); 
        this->insert(std::make_pair("C10", std::make_shared<C10_THDM>())); 

        this->id = WGroup::B;
        this->wilson_type = ContributionType::BSM;
    }
    
};

class BPrimeCoefficientGroup_THDM : public BPrimeCoefficientGroup {
public:
    BPrimeCoefficientGroup_THDM() { this->clear();

        if (!thdm_parameters::is_initialized()) {
            thdm_parameters::init();
        }
        this->insert(std::make_pair("CP1", std::make_shared<CP1_THDM>())); this->insert(std::make_pair("CP2", std::make_shared<CP2_THDM>())); this->insert(std::make_pair("CP3", std::make_shared<CP3_THDM>()));
        this->insert(std::make_pair("CP4", std::make_shared<CP4_THDM>()));  this->insert(std::make_pair("CP5", std::make_shared<CP5_THDM>())); this->insert(std::make_pair("CP6", std::make_shared<CP6_THDM>())); 
        this->insert(std::make_pair("CP7", std::make_shared<CP7_THDM>()));  this->insert(std::make_pair("CP8", std::make_shared<CP8_THDM>()));  this->insert(std::make_pair("CP9", std::make_shared<CP9_THDM>())); 
        this->insert(std::make_pair("CP10", std::make_shared<CP10_THDM>())); this->insert(std::make_pair("CPQ1", std::make_shared<CPQ1_THDM>())); this->insert(std::make_pair("CPQ2", std::make_shared<CPQ2_THDM>())); 

        this->id = WGroup::BPrime;
        this->wilson_type = ContributionType::BSM;
    }

    void set_base_1();
    void set_base_2();
};

class BScalarCoefficientGroup_THDM : public BScalarCoefficientGroup {
public:
    BScalarCoefficientGroup_THDM() { this->clear();

        if (!thdm_parameters::is_initialized()) {
            thdm_parameters::init();
        }
        this->insert(std::make_pair("CQ1", std::make_shared<CQ1_THDM>())); this->insert(std::make_pair("CQ2", std::make_shared<CQ2_THDM>()));

        this->id = WGroup::BScalar;
        this->wilson_type = ContributionType::BSM;
    }

    void set_base_1();
    void set_base_2();
};

class BlnuCoefficientGroup_THDM : public BlnuCoefficientGroup {
public:
    BlnuCoefficientGroup_THDM() { this->clear();

        if (!thdm_parameters::is_initialized()) {
            thdm_parameters::init();
        }
        this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A_THDM>())); 
        this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P_THDM>()));

        this->id = WGroup::Blnu;
        this->wilson_type = ContributionType::BSM;
    }

    void set_base_1() {}
    void set_base_2() {}
};

class BclnuCoefficientGroup_THDM : public BclnuCoefficientGroup {
public:
    BclnuCoefficientGroup_THDM() { this->clear();

        if (!thdm_parameters::is_initialized()) {
            thdm_parameters::init();
        }
        this->insert(std::make_pair("C_V1", std::make_shared<C_V1_THDM>()));
        this->insert(std::make_pair("C_V2", std::make_shared<C_V2_THDM>()));
        this->insert(std::make_pair("C_S1", std::make_shared<C_S1_THDM>()));
        this->insert(std::make_pair("C_S2", std::make_shared<C_S2_THDM>()));
        this->insert(std::make_pair("C_T", std::make_shared<C_T_THDM>()));

        this->id = WGroup::BCLNU;
        this->wilson_type = ContributionType::BSM;
    }

    void set_base_1() {}
    void set_base_2() {}
};

#endif