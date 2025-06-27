#ifndef WILSON_THDM_SUPER_H
#define WILSON_THDM_SUPER_H
#include "Wilson.h"
#include "WilsonGroup.h"
#include "BWilsonGroup.h"
#include "thdm_parameters.h"
#include "Math_THDM.h"
#include "Utils.h"
#include "ChargedCurrentsWilsonGroup.h"


class C1_THDM : public WilsonCoefficient {
public:
    C1_THDM() : WilsonCoefficient("C1_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C1_THDM>(*this);
    }
};

class C2_THDM : public WilsonCoefficient {
public:
    C2_THDM() : WilsonCoefficient("C2_THDM", GroupMapper::str(WGroup::B) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C2_THDM>(*this);
    }
};

class C3_THDM : public WilsonCoefficient {
public:
    C3_THDM();

    
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C3_THDM>(*this);
    }

    static double compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class C4_THDM : public WilsonCoefficient {
public:
    C4_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C4_THDM>(*this);
    }

    static double compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class C5_THDM : public WilsonCoefficient {
public:
    C5_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C5_THDM>(*this);
    }

    static double compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class C6_THDM : public WilsonCoefficient {
public:
    C6_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C6_THDM>(*this);
    }

    static double compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class C7_THDM : public WilsonCoefficient {
public:
    C7_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C7_THDM>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class C8_THDM : public WilsonCoefficient {
public:
    C8_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C8_THDM>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class C9_THDM : public WilsonCoefficient {
public:
    C9_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C9_THDM>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class C10_THDM : public WilsonCoefficient {
public:
    C10_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C10_THDM>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class CQ1_THDM : public WilsonCoefficient {
public:
    CQ1_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ1_THDM>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class CQ2_THDM : public WilsonCoefficient {
public:
    CQ2_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ2_THDM>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class CP1_THDM : public WilsonCoefficient {
public:
    CP1_THDM() : WilsonCoefficient("CP1_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP1_THDM>(*this);
    }
};

class CP2_THDM : public WilsonCoefficient {
public:
    CP2_THDM() : WilsonCoefficient("CP2_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP2_THDM>(*this);
    }
};

class CP3_THDM : public WilsonCoefficient {
public:
    CP3_THDM() : WilsonCoefficient("CP3_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP3_THDM>(*this);
    }
};

class CP4_THDM : public WilsonCoefficient {
public:
    CP4_THDM() : WilsonCoefficient("CP4_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP4_THDM>(*this);
    }
};

class CP5_THDM : public WilsonCoefficient {
public:
    CP5_THDM() : WilsonCoefficient("CP5_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP5_THDM>(*this);
    }
};

class CP6_THDM : public WilsonCoefficient {
public:
    CP6_THDM() : WilsonCoefficient("CP6_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP6_THDM>(*this);
    }
};

// TODO : Implement CP7, CP8, CP9, CP10, CPQ1 and CPQ2 contributions

class CP7_THDM : public WilsonCoefficient {
public:
    CP7_THDM() : WilsonCoefficient("CP7_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP7_THDM>(*this);
    }
};

class CP8_THDM : public WilsonCoefficient {
public:
    CP8_THDM() : WilsonCoefficient("CP8_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP8_THDM>(*this);
    }
};

class CP9_THDM : public WilsonCoefficient {
public:
    CP9_THDM() : WilsonCoefficient("CP9_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP9_THDM>(*this);
    }
};

class CP10_THDM : public WilsonCoefficient {
public:
    CP10_THDM() : WilsonCoefficient("CP10_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP10_THDM>(*this);
    }
};

class CPQ1_THDM : public WilsonCoefficient {
public:
    CPQ1_THDM() : WilsonCoefficient("CPQ1_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ1_THDM>(*this);
    }
};

class CPQ2_THDM : public WilsonCoefficient {
public:
    CPQ2_THDM() : WilsonCoefficient("CPQ2_THDM", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ2_THDM>(*this);
    }
};

class C_V1_THDM : public WilsonCoefficient {
public:
    C_V1_THDM() : WilsonCoefficient("C_V1_THDM", GroupMapper::str(WGroup::BCC) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_THDM>(*this);
    }

};

class C_V2_THDM : public WilsonCoefficient {
public:
    C_V2_THDM() : WilsonCoefficient("C_V2_THDM", GroupMapper::str(WGroup::BCC) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_THDM>(*this);
    }

};

class C_S1_THDM : public WilsonCoefficient {
public:
    C_S1_THDM();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_THDM>(*this);
    }

};

class C_S2_THDM : public WilsonCoefficient {
public:
    C_S2_THDM();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_THDM>(*this);
    }

};

class C_T_THDM : public WilsonCoefficient {
public:
    C_T_THDM() : WilsonCoefficient("C_T_THDM", GroupMapper::str(WGroup::BCC) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_THDM>(*this);
    }

};

class BCoefficientGroup_THDM : public BCoefficientGroup {

public:
    BCoefficientGroup_THDM() {
        this->id = WGroup::B;
        this->wilson_type = ContributionType::BSM;

        if (UseMarty().get()) {
            this->wilson_type = ContributionType::TOTAL;
            return;
        }
        this->clear();

        if (!thdm_parameters::is_initialized()) {
            thdm_parameters::init();
        }
        this->insert(std::make_pair("C1", std::make_shared<C1_THDM>())); this->insert(std::make_pair("C2", std::make_shared<C2_THDM>())); this->insert(std::make_pair("C3", std::make_shared<C3_THDM>()));
        this->insert(std::make_pair("C4", std::make_shared<C4_THDM>()));  this->insert(std::make_pair("C5", std::make_shared<C5_THDM>())); this->insert(std::make_pair("C6", std::make_shared<C6_THDM>())); 
        this->insert(std::make_pair("C7", std::make_shared<C7_THDM>()));  this->insert(std::make_pair("C8", std::make_shared<C8_THDM>()));  this->insert(std::make_pair("C9", std::make_shared<C9_THDM>())); 
        this->insert(std::make_pair("C10", std::make_shared<C10_THDM>())); 

        
    }
    
};

class BPrimeCoefficientGroup_THDM : public BPrimeCoefficientGroup {
public:
    BPrimeCoefficientGroup_THDM() {
        this->id = WGroup::BPrime;
        this->wilson_type = ContributionType::BSM;
        if (UseMarty().get()) {
            this->wilson_type = ContributionType::TOTAL;
            return;
        }
        this->clear();

        if (!thdm_parameters::is_initialized()) {
            thdm_parameters::init();
        }
        this->insert(std::make_pair("CP1", std::make_shared<CP1_THDM>())); this->insert(std::make_pair("CP2", std::make_shared<CP2_THDM>())); this->insert(std::make_pair("CP3", std::make_shared<CP3_THDM>()));
        this->insert(std::make_pair("CP4", std::make_shared<CP4_THDM>()));  this->insert(std::make_pair("CP5", std::make_shared<CP5_THDM>())); this->insert(std::make_pair("CP6", std::make_shared<CP6_THDM>())); 
        this->insert(std::make_pair("CP7", std::make_shared<CP7_THDM>()));  this->insert(std::make_pair("CP8", std::make_shared<CP8_THDM>()));  this->insert(std::make_pair("CP9", std::make_shared<CP9_THDM>())); 
        this->insert(std::make_pair("CP10", std::make_shared<CP10_THDM>())); this->insert(std::make_pair("CPQ1", std::make_shared<CPQ1_THDM>())); this->insert(std::make_pair("CPQ2", std::make_shared<CPQ2_THDM>())); 

        
    }

    void set_base_1();
    void set_base_2();
};

class BScalarCoefficientGroup_THDM : public BScalarCoefficientGroup {
public:
    BScalarCoefficientGroup_THDM() {
        this->id = WGroup::BScalar;
        this->wilson_type = ContributionType::BSM;
        if (UseMarty().get()) {
            this->wilson_type = ContributionType::TOTAL;
            return;
        }
        this->clear();

        if (!thdm_parameters::is_initialized()) {
            thdm_parameters::init();
        }
        this->insert(std::make_pair("CQ1", std::make_shared<CQ1_THDM>())); this->insert(std::make_pair("CQ2", std::make_shared<CQ2_THDM>()));

        
    }

    void set_base_1();
    void set_base_2();
};

// class BlnuCoefficientGroup_THDM : public BlnuCoefficientGroup {
// public:
//     BlnuCoefficientGroup_THDM() { 
//         this->id = WGroup::Blnu;
//         this->wilson_type = ContributionType::BSM;
//         if (UseMarty().get()) {
//             this->wilson_type = ContributionType::TOTAL;
//             return;
//         }
//         this->clear();

//         if (!thdm_parameters::is_initialized()) {
//             thdm_parameters::init();
//         }
//         this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A_THDM>())); 
//         this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P_THDM>()));
//     }

//     void set_base_1() {}
//     void set_base_2() {}
// };

class BclnuCoefficientGroup_THDM : public BclnuCoefficientGroup {
public:
    BclnuCoefficientGroup_THDM() { 
        this->id = WGroup::BCC;
        this->wilson_type = ContributionType::BSM;
        if (UseMarty().get()) {
            this->wilson_type = ContributionType::TOTAL;
            return;
        }
        this->clear();

        if (!thdm_parameters::is_initialized()) {
            thdm_parameters::init();
        }
        this->insert(std::make_pair("C_V1", std::make_shared<C_V1_THDM>()));
        this->insert(std::make_pair("C_V2", std::make_shared<C_V2_THDM>()));
        this->insert(std::make_pair("C_S1", std::make_shared<C_S1_THDM>()));
        this->insert(std::make_pair("C_S2", std::make_shared<C_S2_THDM>()));
        this->insert(std::make_pair("C_T", std::make_shared<C_T_THDM>()));
    }

    void set_base_1() {}
    void set_base_2() {}
};

#endif