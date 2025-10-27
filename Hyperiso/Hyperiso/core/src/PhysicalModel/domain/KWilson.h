#include "Wilson.h"

class CK9 : public WilsonCoefficient {
public:
    CK9();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CK9>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class CK10 : public WilsonCoefficient {
public:
    CK10();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CK10>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

};


class CKQ1 : public WilsonCoefficient {
public:
    CKQ1();
    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CKQ1>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class CKQ2 : public WilsonCoefficient {
public:
    CKQ2();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CKQ2>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};


class CPK9 : public WilsonCoefficient {
public:
    CPK9();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPK9>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class CPK10 : public WilsonCoefficient {
public:
    CPK10();

    
    
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPK10>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
    static double compute_NNLO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

};


class CPKQ1 : public WilsonCoefficient {
public:
    CPKQ1();
    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPKQ1>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};

class CPKQ2 : public WilsonCoefficient {
public:
    CPKQ2();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPKQ2>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);
};
