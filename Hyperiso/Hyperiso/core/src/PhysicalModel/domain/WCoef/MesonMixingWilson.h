#ifndef __MESONMIXINGWILSON_H__
#define __MESONMIXINGWILSON_H__

#include "Wilson.h"

// All coefficients are relative to the SUSY basis for Delta F = 2 operators

/* Bd */

class C_mix_bd_1 : public WilsonCoefficient {
public:
    C_mix_bd_1();

    static complex_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_1>(*this);
    }
};

class C_mix_bd_1_tilde : public WilsonCoefficient {
public:
    C_mix_bd_1_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_1_tilde>(*this);
    }
};

class C_mix_bd_2 : public WilsonCoefficient {
public:
    C_mix_bd_2();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_2>(*this);
    }
};

class C_mix_bd_2_tilde : public WilsonCoefficient {
public:
    C_mix_bd_2_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_2_tilde>(*this);
    }
};

class C_mix_bd_3 : public WilsonCoefficient {
public:
    C_mix_bd_3();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_3>(*this);
    }
};

class C_mix_bd_3_tilde : public WilsonCoefficient {
public:
    C_mix_bd_3_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_3_tilde>(*this);
    }
};

class C_mix_bd_4 : public WilsonCoefficient {
public:
    C_mix_bd_4();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_4>(*this);
    }
};

class C_mix_bd_5 : public WilsonCoefficient {
public:
    C_mix_bd_5();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bd_5>(*this);
    }
};


/* Bs */

class C_mix_bs_1 : public WilsonCoefficient {
public:
    C_mix_bs_1();

    static complex_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_1>(*this);
    }
};

class C_mix_bs_1_tilde : public WilsonCoefficient {
public:
    C_mix_bs_1_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_1_tilde>(*this);
    }
};

class C_mix_bs_2 : public WilsonCoefficient {
public:
    C_mix_bs_2();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_2>(*this);
    }
};

class C_mix_bs_2_tilde : public WilsonCoefficient {
public:
    C_mix_bs_2_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_2_tilde>(*this);
    }
};

class C_mix_bs_3 : public WilsonCoefficient {
public:
    C_mix_bs_3();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_3>(*this);
    }
};

class C_mix_bs_3_tilde : public WilsonCoefficient {
public:
    C_mix_bs_3_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_3_tilde>(*this);
    }
};

class C_mix_bs_4 : public WilsonCoefficient {
public:
    C_mix_bs_4();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_4>(*this);
    }
};

class C_mix_bs_5 : public WilsonCoefficient {
public:
    C_mix_bs_5();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_bs_5>(*this);
    }
};

/* K0 */

class C_mix_sd_1 : public WilsonCoefficient {
public:
    C_mix_sd_1();

    static complex_t compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_sd_1>(*this);
    }
};

class C_mix_sd_1_tilde : public WilsonCoefficient {
public:
    C_mix_sd_1_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_sd_1_tilde>(*this);
    }
};

class C_mix_sd_2 : public WilsonCoefficient {
public:
    C_mix_sd_2();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_sd_2>(*this);
    }
};

class C_mix_sd_2_tilde : public WilsonCoefficient {
public:
    C_mix_sd_2_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_sd_2_tilde>(*this);
    }
};

class C_mix_sd_3 : public WilsonCoefficient {
public:
    C_mix_sd_3();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_sd_3>(*this);
    }
};

class C_mix_sd_3_tilde : public WilsonCoefficient {
public:
    C_mix_sd_3_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_sd_3_tilde>(*this);
    }
};

class C_mix_sd_4 : public WilsonCoefficient {
public:
    C_mix_sd_4();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_sd_4>(*this);
    }
};

class C_mix_sd_5 : public WilsonCoefficient {
public:
    C_mix_sd_5();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_sd_5>(*this);
    }
};


/* D0 */

class C_mix_cu_1 : public WilsonCoefficient {
public:
    C_mix_cu_1();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_cu_1>(*this);
    }
};

class C_mix_cu_1_tilde : public WilsonCoefficient {
public:
    C_mix_cu_1_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_cu_1_tilde>(*this);
    }
};

class C_mix_cu_2 : public WilsonCoefficient {
public:
    C_mix_cu_2();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_cu_2>(*this);
    }
};

class C_mix_cu_2_tilde : public WilsonCoefficient {
public:
    C_mix_cu_2_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_cu_2_tilde>(*this);
    }
};

class C_mix_cu_3 : public WilsonCoefficient {
public:
    C_mix_cu_3();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_cu_3>(*this);
    }
};

class C_mix_cu_3_tilde : public WilsonCoefficient {
public:
    C_mix_cu_3_tilde();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_cu_3_tilde>(*this);
    }
};

class C_mix_cu_4 : public WilsonCoefficient {
public:
    C_mix_cu_4();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_cu_4>(*this);
    }
};

class C_mix_cu_5 : public WilsonCoefficient {
public:
    C_mix_cu_5();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_mix_cu_5>(*this);
    }
};


#endif // __MESONMIXINGWILSON_H__
