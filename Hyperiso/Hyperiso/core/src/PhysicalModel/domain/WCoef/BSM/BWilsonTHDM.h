#ifndef BWILSON_THDM_H
#define BWILSON_THDM_H
#include "Wilson.h"
#include "WilsonGroup.h"
#include "BWilsonGroup.h"
#include "THDMParametersHelper.h"
#include "Math.h"
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

    static double compute_NNLO(const ParamSrc& src);
};

class C4_THDM : public WilsonCoefficient {
public:
    C4_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C4_THDM>(*this);
    }

    static double compute_NLO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);
};

class C5_THDM : public WilsonCoefficient {
public:
    C5_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C5_THDM>(*this);
    }

    static double compute_NNLO(const ParamSrc& src);
};

class C6_THDM : public WilsonCoefficient {
public:
    C6_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C6_THDM>(*this);
    }

    static double compute_NNLO(const ParamSrc& src);
};

class C7_THDM : public WilsonCoefficient {
public:
    C7_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C7_THDM>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);
};

class C8_THDM : public WilsonCoefficient {
public:
    C8_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C8_THDM>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);
};

class C9_THDM : public WilsonCoefficient {
public:
    C9_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C9_THDM>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
};

class C10_THDM : public WilsonCoefficient {
public:
    C10_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C10_THDM>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
};

class CQ1_THDM : public WilsonCoefficient {
public:
    CQ1_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ1_THDM>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class CQ2_THDM : public WilsonCoefficient {
public:
    CQ2_THDM();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ2_THDM>(*this);
    }

    static double compute_LO(const ParamSrc& src);
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

#endif