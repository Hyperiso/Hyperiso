#ifndef MESON_MIXING_WILSON_SUSY_H
#define MESON_MIXING_WILSON_SUSY_H

#include "Wilson.h"
#include "Math.h"

class C_mix_bd_1_SUSY : public WilsonCoefficient {
public:
    C_mix_bd_1_SUSY();

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_1_SUSY>(*this);
    }
};

class C_mix_bd_1_tilde_SUSY : public WilsonCoefficient {
public:
    C_mix_bd_1_tilde_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_1_tilde_SUSY>(*this);
    }
};

class C_mix_bd_2_SUSY : public WilsonCoefficient {
public:
    C_mix_bd_2_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_2_SUSY>(*this);
    }
};

class C_mix_bd_2_tilde_SUSY : public WilsonCoefficient {
public:
    C_mix_bd_2_tilde_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_2_tilde_SUSY>(*this);
    }
};

class C_mix_bd_3_SUSY : public WilsonCoefficient {
public:
    C_mix_bd_3_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_3_SUSY>(*this);
    }
};

class C_mix_bd_3_tilde_SUSY : public WilsonCoefficient {
public:
    C_mix_bd_3_tilde_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_3_tilde_SUSY>(*this);
    }
};

class C_mix_bd_4_SUSY : public WilsonCoefficient {
public:
    C_mix_bd_4_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_4_SUSY>(*this);
    }
};

class C_mix_bd_5_SUSY : public WilsonCoefficient {
public:
    C_mix_bd_5_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_5_SUSY>(*this);
    }
};


/* Bs */

class C_mix_bs_1_SUSY : public WilsonCoefficient {
public:
    C_mix_bs_1_SUSY();

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_1_SUSY>(*this);
    }
};

class C_mix_bs_1_tilde_SUSY : public WilsonCoefficient {
public:
    C_mix_bs_1_tilde_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_1_tilde_SUSY>(*this);
    }
};

class C_mix_bs_2_SUSY : public WilsonCoefficient {
public:
    C_mix_bs_2_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_2_SUSY>(*this);
    }
};

class C_mix_bs_2_tilde_SUSY : public WilsonCoefficient {
public:
    C_mix_bs_2_tilde_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_2_tilde_SUSY>(*this);
    }
};

class C_mix_bs_3_SUSY : public WilsonCoefficient {
public:
    C_mix_bs_3_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_3_SUSY>(*this);
    }
};

class C_mix_bs_3_tilde_SUSY : public WilsonCoefficient {
public:
    C_mix_bs_3_tilde_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_3_tilde_SUSY>(*this);
    }
};

class C_mix_bs_4_SUSY : public WilsonCoefficient {
public:
    C_mix_bs_4_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_4_SUSY>(*this);
    }
};

class C_mix_bs_5_SUSY : public WilsonCoefficient {
public:
    C_mix_bs_5_SUSY();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_5_SUSY>(*this);
    }
};

#endif