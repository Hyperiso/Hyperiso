#include "Wilson.h"

class CK9 : public WilsonCoefficient {
public:
    CK9();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CK9>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
};

class CK10 : public WilsonCoefficient {
public:
    CK10();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CK10>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);

};


class CKQ1 : public WilsonCoefficient {
public:
    CKQ1();
    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CKQ1>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class CKQ2 : public WilsonCoefficient {
public:
    CKQ2();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CKQ2>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};


class CPK9 : public WilsonCoefficient {
public:
    CPK9();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPK9>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
};

class CPK10 : public WilsonCoefficient {
public:
    CPK10();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPK10>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
    static double compute_NNLO(const ParamSrc& src);

};


class CPKQ1 : public WilsonCoefficient {
public:
    CPKQ1();
    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPKQ1>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class CPKQ2 : public WilsonCoefficient {
public:
    CPKQ2();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPKQ2>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class CK_L : public WilsonCoefficient {
public:
    CK_L();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CK_L>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_NLO(const ParamSrc& src);
};