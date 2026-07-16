#ifndef MESON_MIXING_WILSON_THDM_H
#define MESON_MIXING_WILSON_THDM_H

#include "Wilson.h"
#include "Math.h"

class C_mix_bd_1_THDM : public WilsonCoefficient {
public:
    C_mix_bd_1_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_1_THDM>(*this);
    }
};

class C_mix_bd_1_tilde_THDM : public WilsonCoefficient {
public:
    C_mix_bd_1_tilde_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_1_tilde_THDM>(*this);
    }
};

class C_mix_bd_2_THDM : public WilsonCoefficient {
public:
    C_mix_bd_2_THDM();

    static double compute_LO(const ParamSrc& src);
    
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_2_THDM>(*this);
    }
};

class C_mix_bd_2_tilde_THDM : public WilsonCoefficient {
public:
    C_mix_bd_2_tilde_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_2_tilde_THDM>(*this);
    }
};

class C_mix_bd_3_THDM : public WilsonCoefficient {
public:
    C_mix_bd_3_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_3_THDM>(*this);
    }
};

class C_mix_bd_3_tilde_THDM : public WilsonCoefficient {
public:
    C_mix_bd_3_tilde_THDM();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_3_tilde_THDM>(*this);
    }
};

class C_mix_bd_4_THDM : public WilsonCoefficient {
public:
    C_mix_bd_4_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_4_THDM>(*this);
    }
};

class C_mix_bd_5_THDM : public WilsonCoefficient {
public:
    C_mix_bd_5_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_5_THDM>(*this);
    }
};


/* Bs */

class C_mix_bs_1_THDM : public WilsonCoefficient {
public:
    C_mix_bs_1_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_1_THDM>(*this);
    }
};

class C_mix_bs_1_tilde_THDM : public WilsonCoefficient {
public:
    C_mix_bs_1_tilde_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_1_tilde_THDM>(*this);
    }
};

class C_mix_bs_2_THDM : public WilsonCoefficient {
public:
    C_mix_bs_2_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_2_THDM>(*this);
    }
};

class C_mix_bs_2_tilde_THDM : public WilsonCoefficient {
public:
    C_mix_bs_2_tilde_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_2_tilde_THDM>(*this);
    }
};

class C_mix_bs_3_THDM : public WilsonCoefficient {
public:
    C_mix_bs_3_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_3_THDM>(*this);
    }
};

class C_mix_bs_3_tilde_THDM : public WilsonCoefficient {
public:
    C_mix_bs_3_tilde_THDM();

    static double compute_LO(const ParamSrc& src);
    
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_3_tilde_THDM>(*this);
    }
};

class C_mix_bs_4_THDM : public WilsonCoefficient {
public:
    C_mix_bs_4_THDM();

    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_4_THDM>(*this);
    }
};

class C_mix_bs_5_THDM : public WilsonCoefficient {
public:
    C_mix_bs_5_THDM();

    static double compute_LO(const ParamSrc& src);
    
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_5_THDM>(*this);
    }
};

#endif