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


#endif