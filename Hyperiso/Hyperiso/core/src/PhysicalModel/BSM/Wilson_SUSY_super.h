#ifndef WILSON_SUSY_SUPER_H
#define WILSON_SUSY_SUPER_H
#include "WilsonSuper.h"
#include "BWilsonGroupSuper.h"
#include "susy_parameters.h"
#include "Math_THDM.h"
#include "Utils.h"


// TODO : adapt with storage_block constructor as in THDM
class C1_susy : public WilsonCoefficient {
public:
    C1_susy();

    
     
    
    static scalar_t compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C1_susy>(*this);
    }
};

class C2_susy : public WilsonCoefficient {
public:
    C2_susy() : WilsonCoefficient("C2_SUSY", GroupMapper::str(WGroup::B) + "_MATCH") { }

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C2_susy>(*this);
    }
};

class C3_susy : public WilsonCoefficient {
public:
    C3_susy();


    static scalar_t compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C3_susy>(*this);
    }
};

class C4_susy : public WilsonCoefficient {
public:
    C4_susy() : WilsonCoefficient() { this->set_name("C4_SUSY"); }

    static scalar_t compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static scalar_t compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C4_susy>(*this);
    }
};

class C5_susy : public WilsonCoefficient {
public:
    C5_susy() : WilsonCoefficient() { this->set_name("C5_SUSY"); }

    static scalar_t compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C5_susy>(*this);
    }
};

class C6_susy : public WilsonCoefficient {
public:
    C6_susy() : WilsonCoefficient() { this->set_name("C6_SUSY"); }

    static scalar_t compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C6_susy>(*this);
    }
};

class C7_susy : public WilsonCoefficient {
public:
    C7_susy() : WilsonCoefficient() { this->set_name("C7_SUSY"); }

    static scalar_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static scalar_t compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static scalar_t compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C7_susy>(*this);
    }
};

class C8_susy : public WilsonCoefficient {
public:
    C8_susy() : WilsonCoefficient() { this->set_name("C8_SUSY"); }

    static scalar_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static scalar_t compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static scalar_t compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C8_susy>(*this);
    }
};

class C9_susy : public WilsonCoefficient {
public:
    C9_susy() : WilsonCoefficient() { this->set_name("C9_SUSY"); }

    static scalar_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static scalar_t compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C9_susy>(*this);
    }
};

class C10_susy : public WilsonCoefficient {
public:
    C10_susy() : WilsonCoefficient() { this->set_name("C10_SUSY"); }

    static scalar_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static scalar_t compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C10_susy>(*this);
    }
};

class CQ1_susy : public WilsonCoefficient {
public:
    CQ1_susy() : WilsonCoefficient() { this->set_name("CQ1_SUSY"); }



    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ1_susy>(*this);
    }
};

class CQ2_susy : public WilsonCoefficient {
public:
    CQ2_susy() : WilsonCoefficient() { this->set_name("CQ2_SUSY"); }



    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ2_susy>(*this);
    }
};

class CP1_susy : public WilsonCoefficient {
public:
    CP1_susy() : WilsonCoefficient() { this->set_name("CP1_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP1_susy>(*this);
    }
};

class CP2_susy : public WilsonCoefficient {
public:
    CP2_susy() : WilsonCoefficient() { this->set_name("CP2_SUSY"); }



    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP2_susy>(*this);
    }
};

class CP3_susy : public WilsonCoefficient {
public:
    CP3_susy() : WilsonCoefficient() { this->set_name("CP3_SUSY"); }



    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP3_susy>(*this);
    }

};

class CP4_susy : public WilsonCoefficient {
public:
    CP4_susy() : WilsonCoefficient() { this->set_name("CP4_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP4_susy>(*this);
    }
};

class CP5_susy : public WilsonCoefficient {
public:
    CP5_susy() : WilsonCoefficient() { this->set_name("CP5_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP5_susy>(*this);
    }
};

class CP6_susy : public WilsonCoefficient {
public:
    CP6_susy() : WilsonCoefficient() { this->set_name("CP6_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP6_susy>(*this);
    }
};

class CP7_susy : public WilsonCoefficient {
public:
    CP7_susy() : WilsonCoefficient() { this->set_name("CP7_SUSY"); }

    static scalar_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP7_susy>(*this);
    }
};

class CP8_susy : public WilsonCoefficient {
public:
    CP8_susy() : WilsonCoefficient() { this->set_name("CP8_SUSY"); }

    static scalar_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP8_susy>(*this);
    }
};

class CP9_susy : public WilsonCoefficient {
public:
    CP9_susy() : WilsonCoefficient() { this->set_name("CP9_SUSY"); }

    static scalar_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP9_susy>(*this);
    }
};

class CP10_susy : public WilsonCoefficient {
public:
    CP10_susy() : WilsonCoefficient() { this->set_name("CP10_SUSY"); }

    static scalar_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP10_susy>(*this);
    }
};

class CPQ1_susy : public WilsonCoefficient {
public:
    CPQ1_susy() : WilsonCoefficient() { this->set_name("CPQ1_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ1_susy>(*this);
    }
};

class CPQ2_susy : public WilsonCoefficient {
public:
    CPQ2_susy() : WilsonCoefficient() { this->set_name("CPQ2_SUSY"); }



    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ2_susy>(*this);
    }
};

class C_Blnu_A_SUSY : public WilsonCoefficient {
public:
    C_Blnu_A_SUSY() : WilsonCoefficient() { this->set_name("C_Blnu_A_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_Blnu_A_SUSY>(*this);
    }
};

class C_Blnu_P_SUSY : public WilsonCoefficient {
public:
    C_Blnu_P_SUSY() : WilsonCoefficient() { this->set_name("C_Blnu_P_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_Blnu_P_SUSY>(*this);
    }
};

class C_V1_SUSY : public WilsonCoefficient {
public:
    C_V1_SUSY() : WilsonCoefficient() { this->set_name("C_V1_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_SUSY>(*this);
    }
};

class C_V2_SUSY : public WilsonCoefficient {
public:
    C_V2_SUSY() : WilsonCoefficient() { this->set_name("C_V2_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_SUSY>(*this);
    }
};

class C_S1_SUSY : public WilsonCoefficient {
public:
    C_S1_SUSY() : WilsonCoefficient() { this->set_name("C_S1_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_SUSY>(*this);
    }
};

class C_S2_SUSY : public WilsonCoefficient {
public:
    C_S2_SUSY() : WilsonCoefficient() { this->set_name("C_S2_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_SUSY>(*this);
    }
};

class C_T_SUSY : public WilsonCoefficient {
public:
    C_T_SUSY() : WilsonCoefficient() { this->set_name("C_T_SUSY"); }

 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_SUSY>(*this);
    }
};

class BCoefficientGroup_susy : public BCoefficientGroup {

public:
    BCoefficientGroup_susy() { this->clear();

        if (!susy_parameters::is_initialized()) {
            susy_parameters::init();
        }

        this->insert(std::make_pair("C1", std::make_shared<C1_susy>())); this->insert(std::make_pair("C2", std::make_shared<C2_susy>())); this->insert(std::make_pair("C3", std::make_shared<C3_susy>()));
        this->insert(std::make_pair("C4", std::make_shared<C4_susy>()));  this->insert(std::make_pair("C5", std::make_shared<C5_susy>())); this->insert(std::make_pair("C6", std::make_shared<C6_susy>())); 
        this->insert(std::make_pair("C7", std::make_shared<C7_susy>()));  this->insert(std::make_pair("C8", std::make_shared<C8_susy>()));  this->insert(std::make_pair("C9", std::make_shared<C9_susy>())); 
        this->insert(std::make_pair("C10", std::make_shared<C10_susy>())); 

        this->id = WGroup::B;
        this->wilson_type = ContributionType::BSM;
    }
};

class BPrimeCoefficientGroup_susy : public BPrimeCoefficientGroup {
public:
    BPrimeCoefficientGroup_susy() { this->clear();

        if (!susy_parameters::is_initialized()) {
            susy_parameters::init();
        }
        this->insert(std::make_pair("CP1", std::make_shared<CP1_susy>())); this->insert(std::make_pair("CP2", std::make_shared<CP2_susy>())); this->insert(std::make_pair("CP3", std::make_shared<CP3_susy>()));
        this->insert(std::make_pair("CP4", std::make_shared<CP4_susy>()));  this->insert(std::make_pair("CP5", std::make_shared<CP5_susy>())); this->insert(std::make_pair("CP6", std::make_shared<CP6_susy>())); 
        this->insert(std::make_pair("CP7", std::make_shared<CP7_susy>()));  this->insert(std::make_pair("CP8", std::make_shared<CP8_susy>()));  this->insert(std::make_pair("CP9", std::make_shared<CP9_susy>())); 
        this->insert(std::make_pair("CP10", std::make_shared<CP10_susy>())); this->insert(std::make_pair("CPQ1", std::make_shared<CPQ1_susy>())); this->insert(std::make_pair("CPQ2", std::make_shared<CPQ2_susy>())); 

        this->id = WGroup::BPrime;
        this->wilson_type = ContributionType::BSM;
    }
};

class BScalarCoefficientGroup_susy : public BScalarCoefficientGroup {
public:
    BScalarCoefficientGroup_susy() : BScalarCoefficientGroup() { this->clear();

        if (!susy_parameters::is_initialized()) {
            susy_parameters::init();
        }
        this->insert(std::make_pair("CQ1", std::make_shared<CQ1_susy>())); this->insert(std::make_pair("CQ2", std::make_shared<CQ2_susy>()));

        this->id = WGroup::BScalar;
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

        if (!susy_parameters::is_initialized()) {
            susy_parameters::init();
        }
        this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A_SUSY>())); 
        this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P_SUSY>()));

        this->id = WGroup::Blnu;
        this->wilson_type = ContributionType::BSM;
    }
};

class BclnuCoefficientGroup_SUSY : public BclnuCoefficientGroup {
public:
    BclnuCoefficientGroup_SUSY() { this->clear();

        if (!susy_parameters::is_initialized()) {
            susy_parameters::init();
        }
        this->insert(std::make_pair("C_V1", std::make_shared<C_V1_SUSY>()));
        this->insert(std::make_pair("C_V2", std::make_shared<C_V2_SUSY>()));
        this->insert(std::make_pair("C_S1", std::make_shared<C_S1_SUSY>()));
        this->insert(std::make_pair("C_S2", std::make_shared<C_S2_SUSY>()));
        this->insert(std::make_pair("C_T", std::make_shared<C_T_SUSY>()));

        this->id = WGroup::BCLNU;
        this->wilson_type = ContributionType::BSM;
    }

    void set_base_1() {}
    void set_base_2() {}
};

#endif