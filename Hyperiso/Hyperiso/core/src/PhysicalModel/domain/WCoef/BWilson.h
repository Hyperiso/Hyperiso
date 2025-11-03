#ifndef BWILSON_H
#define BWILSON_H

#include "Wilson.h"

class C1 : public WilsonCoefficient {
public:
    C1();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C1>(*this);
    }

};

class C2 : public WilsonCoefficient {
public:
    C2();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C2>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);

};

class C3 : public WilsonCoefficient {
public:
    C3();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C3>(*this);
    }

    static double compute_NNLO(const ParamSrc& src);
};

class C4 : public WilsonCoefficient {
public:
    C4();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C4>(*this);
    }

    static double compute_NLO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);
};

class C5 : public WilsonCoefficient {
public:
    C5();

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C5>(*this);
    }

    static double compute_NNLO(const ParamSrc& src);
};

class C6 : public WilsonCoefficient {
public:
    C6();

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C6>(*this);
    }


    static double compute_NNLO(const ParamSrc& src);
};

class C7 : public WilsonCoefficient {
public:
    C7();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C7>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);
};

class C8 : public WilsonCoefficient {
public:
    C8();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C8>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);
};

class C9 : public WilsonCoefficient {
public:
    C9();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C9>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
};

class C10 : public WilsonCoefficient {
public:
    C10();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C10>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);

};


class CQ1 : public WilsonCoefficient {
public:
    CQ1();
    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ1>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class CQ2 : public WilsonCoefficient {
public:
    CQ2();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ2>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class CP1 : public WilsonCoefficient {
public:
    CP1() : WilsonCoefficient("CP1", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP1>(*this);
    }

};

class CP2 : public WilsonCoefficient {
public:
    CP2() : WilsonCoefficient("CP2", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP2>(*this);
    }

};

class CP3 : public WilsonCoefficient {
public:
    CP3() : WilsonCoefficient("CP3", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP3>(*this);
    }

};

class CP4 : public WilsonCoefficient {
public:
    CP4() : WilsonCoefficient("CP4", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP4>(*this);
    }

};

class CP5 : public WilsonCoefficient {
public:
    CP5() : WilsonCoefficient("CP5", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP5>(*this);
    }

};

class CP6 : public WilsonCoefficient {
public:
    CP6() : WilsonCoefficient("CP6", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP6>(*this);
    }

};

class CP7 : public WilsonCoefficient {
public:
    CP7();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP7>(*this);
    }
    static double compute_LO(const ParamSrc& src);
};

class CP8 : public WilsonCoefficient {
public:
    CP8();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP8>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class CP9 : public WilsonCoefficient {
public:
    CP9() : WilsonCoefficient("CP9", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP9>(*this);
    }

};

class CP10 : public WilsonCoefficient {
public:
    CP10() : WilsonCoefficient("CP10", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP10>(*this);
    }

};

class CPQ1 : public WilsonCoefficient {
public:
    CPQ1() : WilsonCoefficient("CPQ1", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ1>(*this);
    }

    int gen{2};
};

class CPQ2 : public WilsonCoefficient {
public:
    CPQ2() : WilsonCoefficient("CPQ2", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ2>(*this);
    }

    int gen{2};
};

#endif