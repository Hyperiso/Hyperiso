#include "WilsonCoefficientRegistry.h"
#include "BWilson.h"
#include "ChargedCurrentsWilsonGroup.h"
#include "MesonMixingWilsonGroup.h" //TODO
#include "MartyWilson.h"
#include "Wilson_SUSY_super.h"
#include "Wilson_THDM_super.h"

#include "Include.h"

#define REG(c,m,b,body) \
  reg.register_creator((c),(m),(b), [](const BuildContext& ctx, WCoef coef) -> CoefPtr { return (body); })


static LhaID flhaid(WCoef c, QCDOrder ord, ContributionType ct) {
    return WCoefMapper::flha_full(c, ord, ct);
}

static CoefPtr make_marty(const BuildContext& ctx, WCoef c) {
    std::string name = (ctx.contrib==ContributionType::SM) ? "SM" : ctx.adapters.marty_model_name->get();
    fs::path path    = (ctx.contrib==ContributionType::SM) ? ctx.adapters.sm_path : ctx.adapters.marty_model_path->get();
    std::string block= GroupMapper::str(ctx.group_id, ScaleType::MATCHING);
    LhaID id        = flhaid(c, QCDOrder::LO, ctx.contrib);
    MartyWilsonConfig cfg { name, id, block, path, ctx.adapters.marty_proxy };
    return std::make_shared<MartyWilson>(cfg);
}

void register_B(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::C1,  SM, Builtin, std::make_shared<C1>());
    REG(WCoef::C2,  SM, Builtin, std::make_shared<C2>());
    REG(WCoef::C3,  SM, Builtin, std::make_shared<C3>());
    REG(WCoef::C4,  SM, Builtin, std::make_shared<C4>());
    REG(WCoef::C5,  SM, Builtin, std::make_shared<C5>());
    REG(WCoef::C6,  SM, Builtin, std::make_shared<C6>());
    REG(WCoef::C7,  SM, Builtin, std::make_shared<C7>());
    REG(WCoef::C8,  SM, Builtin, std::make_shared<C8>());
    REG(WCoef::C9,  SM, Builtin, std::make_shared<C9>());
    REG(WCoef::C10, SM, Builtin, std::make_shared<C10>());

    REG(WCoef::C1,  SUSY, Builtin, std::make_shared<C1_susy>());
    REG(WCoef::C2,  SUSY, Builtin, std::make_shared<C2_susy>());
    REG(WCoef::C3,  SUSY, Builtin, std::make_shared<C3_susy>());
    REG(WCoef::C4,  SUSY, Builtin, std::make_shared<C4_susy>());
    REG(WCoef::C5,  SUSY, Builtin, std::make_shared<C5_susy>());
    REG(WCoef::C6,  SUSY, Builtin, std::make_shared<C6_susy>());
    REG(WCoef::C7,  SUSY, Builtin, std::make_shared<C7_susy>());
    REG(WCoef::C8,  SUSY, Builtin, std::make_shared<C8_susy>());
    REG(WCoef::C9,  SUSY, Builtin, std::make_shared<C9_susy>());
    REG(WCoef::C10, SUSY, Builtin, std::make_shared<C10_susy>());

    REG(WCoef::C1,  THDM, Builtin, std::make_shared<C1_THDM>());
    REG(WCoef::C2,  THDM, Builtin, std::make_shared<C2_THDM>());
    REG(WCoef::C3,  THDM, Builtin, std::make_shared<C3_THDM>());
    REG(WCoef::C4,  THDM, Builtin, std::make_shared<C4_THDM>());
    REG(WCoef::C5,  THDM, Builtin, std::make_shared<C5_THDM>());
    REG(WCoef::C6,  THDM, Builtin, std::make_shared<C6_THDM>());
    REG(WCoef::C7,  THDM, Builtin, std::make_shared<C7_THDM>());
    REG(WCoef::C8,  THDM, Builtin, std::make_shared<C8_THDM>());
    REG(WCoef::C9,  THDM, Builtin, std::make_shared<C9_THDM>());
    REG(WCoef::C10, THDM, Builtin, std::make_shared<C10_THDM>());

    REG(WCoef::C1,  SM,   Marty, make_marty(ctx, WCoef::C1));
    REG(WCoef::C2,  SM,   Marty, make_marty(ctx, WCoef::C2));
    REG(WCoef::C3,  SM,   Marty, make_marty(ctx, WCoef::C3));
    REG(WCoef::C4,  SM,   Marty, make_marty(ctx, WCoef::C4));
    REG(WCoef::C5,  SM,   Marty, make_marty(ctx, WCoef::C5));
    REG(WCoef::C6,  SM,   Marty, make_marty(ctx, WCoef::C6));
    REG(WCoef::C7,  SM,   Marty, make_marty(ctx, WCoef::C7));
    REG(WCoef::C8,  SM,   Marty, make_marty(ctx, WCoef::C8));
    REG(WCoef::C9,  SM,   Marty, make_marty(ctx, WCoef::C9));
    REG(WCoef::C10, SM,   Marty, make_marty(ctx, WCoef::C10));

    // TODO Marty custom, THDM, SUSY
}

void register_BPrime(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::CP1,  SM, Builtin, std::make_shared<CP1>());
    REG(WCoef::CP2,  SM, Builtin, std::make_shared<CP2>());
    REG(WCoef::CP3,  SM, Builtin, std::make_shared<CP3>());
    REG(WCoef::CP4,  SM, Builtin, std::make_shared<CP4>());
    REG(WCoef::CP5,  SM, Builtin, std::make_shared<CP5>());
    REG(WCoef::CP6,  SM, Builtin, std::make_shared<CP6>());
    REG(WCoef::CP7,  SM, Builtin, std::make_shared<CP7>());
    REG(WCoef::CP8,  SM, Builtin, std::make_shared<CP8>());
    REG(WCoef::CP9,  SM, Builtin, std::make_shared<CP9>());
    REG(WCoef::CP10, SM, Builtin, std::make_shared<CP10>());
    REG(WCoef::CPQ1, SM, Builtin, std::make_shared<CPQ1>());
    REG(WCoef::CPQ2, SM, Builtin, std::make_shared<CPQ2>());

    REG(WCoef::CP1,  SUSY, Builtin, std::make_shared<CP1_susy>());
    REG(WCoef::CP2,  SUSY, Builtin, std::make_shared<CP2_susy>());
    REG(WCoef::CP3,  SUSY, Builtin, std::make_shared<CP3_susy>());
    REG(WCoef::CP4,  SUSY, Builtin, std::make_shared<CP4_susy>());
    REG(WCoef::CP5,  SUSY, Builtin, std::make_shared<CP5_susy>());
    REG(WCoef::CP6,  SUSY, Builtin, std::make_shared<CP6_susy>());
    REG(WCoef::CP7,  SUSY, Builtin, std::make_shared<CP7_susy>());
    REG(WCoef::CP8,  SUSY, Builtin, std::make_shared<CP8_susy>());
    REG(WCoef::CP9,  SUSY, Builtin, std::make_shared<CP9_susy>());
    REG(WCoef::CP10, SUSY, Builtin, std::make_shared<CP10_susy>());
    REG(WCoef::CPQ1, SUSY, Builtin, std::make_shared<CPQ1_susy>());
    REG(WCoef::CPQ2, SUSY, Builtin, std::make_shared<CPQ2_susy>());

    // THDM builtin
    REG(WCoef::CP1,  THDM, Builtin, std::make_shared<CP1_THDM>());
    REG(WCoef::CP2,  THDM, Builtin, std::make_shared<CP2_THDM>());
    REG(WCoef::CP3,  THDM, Builtin, std::make_shared<CP3_THDM>());
    REG(WCoef::CP4,  THDM, Builtin, std::make_shared<CP4_THDM>());
    REG(WCoef::CP5,  THDM, Builtin, std::make_shared<CP5_THDM>());
    REG(WCoef::CP6,  THDM, Builtin, std::make_shared<CP6_THDM>());
    REG(WCoef::CP7,  THDM, Builtin, std::make_shared<CP7_THDM>());
    REG(WCoef::CP8,  THDM, Builtin, std::make_shared<CP8_THDM>());
    REG(WCoef::CP9,  THDM, Builtin, std::make_shared<CP9_THDM>());
    REG(WCoef::CP10, THDM, Builtin, std::make_shared<CP10_THDM>());
    REG(WCoef::CPQ1, THDM, Builtin, std::make_shared<CPQ1_THDM>());
    REG(WCoef::CPQ2, THDM, Builtin, std::make_shared<CPQ2_THDM>());

    for (WCoef c : {WCoef::CP1,WCoef::CP2,WCoef::CP3,WCoef::CP4,WCoef::CP5,
                    WCoef::CP6,WCoef::CP7,WCoef::CP8,WCoef::CP9,WCoef::CP10,
                    WCoef::CPQ1,WCoef::CPQ2}) {
        REG(c, Model::SM, Backend::Marty, make_marty(ctx, coef));
    }

    //TODO Marty other cases
}

void register_BScalar(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::CQ1, SM, Builtin, std::make_shared<CQ1>());
    REG(WCoef::CQ2, SM, Builtin, std::make_shared<CQ2>());

    REG(WCoef::CQ1, SUSY, Builtin, std::make_shared<CQ1_susy>());
    REG(WCoef::CQ2, SUSY, Builtin, std::make_shared<CQ2_susy>());

    REG(WCoef::CQ1, THDM, Builtin, std::make_shared<CQ1_THDM>());
    REG(WCoef::CQ2, THDM, Builtin, std::make_shared<CQ2_THDM>());

    REG(WCoef::CQ1, SM, Marty, make_marty(ctx, WCoef::CQ1));
    REG(WCoef::CQ2, SM, Marty, make_marty(ctx, WCoef::CQ2));

    //TODO : marty other cases
}

void register_BCC(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::C_V1, SM, Builtin, std::make_shared<C_V1>());
    REG(WCoef::C_V2, SM, Builtin, std::make_shared<C_V2>());
    REG(WCoef::C_S1, SM, Builtin, std::make_shared<C_S1>());
    REG(WCoef::C_S2, SM, Builtin, std::make_shared<C_S2>());
    REG(WCoef::C_T,  SM, Builtin, std::make_shared<C_T>());

    REG(WCoef::C_V1, SUSY, Builtin, std::make_shared<C_V1_SUSY>());
    REG(WCoef::C_V2, SUSY, Builtin, std::make_shared<C_V2_SUSY>());
    REG(WCoef::C_S1, SUSY, Builtin, std::make_shared<C_S1_SUSY>());
    REG(WCoef::C_S2, SUSY, Builtin, std::make_shared<C_S2_SUSY>());
    REG(WCoef::C_T,  SUSY, Builtin, std::make_shared<C_T_SUSY>());

    REG(WCoef::C_V1, THDM, Builtin, std::make_shared<C_V1_THDM>());
    REG(WCoef::C_V2, THDM, Builtin, std::make_shared<C_V2_THDM>());
    REG(WCoef::C_S1, THDM, Builtin, std::make_shared<C_S1_THDM>());
    REG(WCoef::C_S2, THDM, Builtin, std::make_shared<C_S2_THDM>());
    REG(WCoef::C_T,  THDM, Builtin, std::make_shared<C_T_THDM>());

    for (WCoef c : {WCoef::C_V1, WCoef::C_V2, WCoef::C_S1, WCoef::C_S2, WCoef::C_T})
        REG(c, Model::SM, Backend::Marty, make_marty(ctx, coef));

        //TODO : marty other cases
}
