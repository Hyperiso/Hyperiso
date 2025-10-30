#include "WilsonCoefficientRegistry.h"
#include "BWilson.h"
#include "ChargedCurrentsWilsonGroup.h"
#include "DChargedCurrentsWilsonGroup.h"
#include "KChargedCurrentsWilsonGroup.h"
#include "PIChargedCurrentsWilsonGroup.h"
#include "MesonMixingWilsonGroup.h"
#include "MartyWilson.h"
#include "BWilsonSUSY.h"
#include "BWilsonTHDM.h"
#include "ChargedCurrentWilsonTHDM.h"
#include "DChargedCurrentWilsonTHDM.h"
#include "KChargedCurrentWilsonTHDM.h"
#include "PIChargedCurrentWilsonTHDM.h"
#include "ChargedCurrentWilsonSUSY.h"
#include "DChargedCurrentWilsonSUSY.h"
#include "KChargedCurrentWilsonSUSY.h"
#include "PIChargedCurrentWilsonSUSY.h"
#include "KWilson.h"

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

void register_CC_bc(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::C_V1_bc, SM, Builtin, std::make_shared<C_V1_bc>());
    REG(WCoef::C_V2_bc, SM, Builtin, std::make_shared<C_V2_bc>());
    REG(WCoef::C_S1_bc, SM, Builtin, std::make_shared<C_S1_bc>());
    REG(WCoef::C_S2_bc, SM, Builtin, std::make_shared<C_S2_bc>());
    REG(WCoef::C_T_bc,  SM, Builtin, std::make_shared<C_T_bc>());

    REG(WCoef::C_V1_bc, SUSY, Builtin, std::make_shared<C_V1_bc_SUSY>());
    REG(WCoef::C_V2_bc, SUSY, Builtin, std::make_shared<C_V2_bc_SUSY>());
    REG(WCoef::C_S1_bc, SUSY, Builtin, std::make_shared<C_S1_bc_SUSY>());
    REG(WCoef::C_S2_bc, SUSY, Builtin, std::make_shared<C_S2_bc_SUSY>());
    REG(WCoef::C_T_bc,  SUSY, Builtin, std::make_shared<C_T_bc_SUSY>());

    REG(WCoef::C_V1_bc, THDM, Builtin, std::make_shared<C_V1_bc_THDM>());
    REG(WCoef::C_V2_bc, THDM, Builtin, std::make_shared<C_V2_bc_THDM>());
    REG(WCoef::C_S1_bc, THDM, Builtin, std::make_shared<C_S1_bc_THDM>());
    REG(WCoef::C_S2_bc, THDM, Builtin, std::make_shared<C_S2_bc_THDM>());
    REG(WCoef::C_T_bc,  THDM, Builtin, std::make_shared<C_T_bc_THDM>());

    for (WCoef c : {WCoef::C_V1_bc, WCoef::C_V2_bc, WCoef::C_S1_bc, WCoef::C_S2_bc, WCoef::C_T_bc})
        REG(c, Model::SM, Backend::Marty, make_marty(ctx, coef));

        //TODO : marty other cases
}

void register_CC_bu(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::C_V1_bu, SM, Builtin, std::make_shared<C_V1_bu>());
    REG(WCoef::C_V2_bu, SM, Builtin, std::make_shared<C_V2_bu>());
    REG(WCoef::C_S1_bu, SM, Builtin, std::make_shared<C_S1_bu>());
    REG(WCoef::C_S2_bu, SM, Builtin, std::make_shared<C_S2_bu>());
    REG(WCoef::C_T_bu,  SM, Builtin, std::make_shared<C_T_bu>());

    REG(WCoef::C_V1_bu, SUSY, Builtin, std::make_shared<C_V1_bu_SUSY>());
    REG(WCoef::C_V2_bu, SUSY, Builtin, std::make_shared<C_V2_bu_SUSY>());
    REG(WCoef::C_S1_bu, SUSY, Builtin, std::make_shared<C_S1_bu_SUSY>());
    REG(WCoef::C_S2_bu, SUSY, Builtin, std::make_shared<C_S2_bu_SUSY>());
    REG(WCoef::C_T_bu,  SUSY, Builtin, std::make_shared<C_T_bu_SUSY>());

    REG(WCoef::C_V1_bu, THDM, Builtin, std::make_shared<C_V1_bu_THDM>());
    REG(WCoef::C_V2_bu, THDM, Builtin, std::make_shared<C_V2_bu_THDM>());
    REG(WCoef::C_S1_bu, THDM, Builtin, std::make_shared<C_S1_bu_THDM>());
    REG(WCoef::C_S2_bu, THDM, Builtin, std::make_shared<C_S2_bu_THDM>());
    REG(WCoef::C_T_bu,  THDM, Builtin, std::make_shared<C_T_bu_THDM>());

    for (WCoef c : {WCoef::C_V1_bu, WCoef::C_V2_bu, WCoef::C_S1_bu, WCoef::C_S2_bu, WCoef::C_T_bu})
        REG(c, Model::SM, Backend::Marty, make_marty(ctx, coef));

        //TODO : marty other cases
}

void register_CC_cs(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::C_V1_cs, SM, Builtin, std::make_shared<C_V1_cs>());
    REG(WCoef::C_V2_cs, SM, Builtin, std::make_shared<C_V2_cs>());
    REG(WCoef::C_S1_cs, SM, Builtin, std::make_shared<C_S1_cs>());
    REG(WCoef::C_S2_cs, SM, Builtin, std::make_shared<C_S2_cs>());
    REG(WCoef::C_T_cs,  SM, Builtin, std::make_shared<C_T_cs>());

    REG(WCoef::C_V1_cs, SUSY, Builtin, std::make_shared<C_V1_cs_SUSY>());
    REG(WCoef::C_V2_cs, SUSY, Builtin, std::make_shared<C_V2_cs_SUSY>());
    REG(WCoef::C_S1_cs, SUSY, Builtin, std::make_shared<C_S1_cs_SUSY>());
    REG(WCoef::C_S2_cs, SUSY, Builtin, std::make_shared<C_S2_cs_SUSY>());
    REG(WCoef::C_T_cs,  SUSY, Builtin, std::make_shared<C_T_cs_SUSY>());

    REG(WCoef::C_V1_cs, THDM, Builtin, std::make_shared<C_V1_cs_THDM>());
    REG(WCoef::C_V2_cs, THDM, Builtin, std::make_shared<C_V2_cs_THDM>());
    REG(WCoef::C_S1_cs, THDM, Builtin, std::make_shared<C_S1_cs_THDM>());
    REG(WCoef::C_S2_cs, THDM, Builtin, std::make_shared<C_S2_cs_THDM>());
    REG(WCoef::C_T_cs,  THDM, Builtin, std::make_shared<C_T_cs_THDM>());

    for (WCoef c : {WCoef::C_V1_cs, WCoef::C_V2_cs, WCoef::C_S1_cs, WCoef::C_S2_cs, WCoef::C_T_cs})
        REG(c, Model::SM, Backend::Marty, make_marty(ctx, coef));

        //TODO : marty other cases
}

void register_CC_cd(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::C_V1_cd, SM, Builtin, std::make_shared<C_V1_cd>());
    REG(WCoef::C_V2_cd, SM, Builtin, std::make_shared<C_V2_cd>());
    REG(WCoef::C_S1_cd, SM, Builtin, std::make_shared<C_S1_cd>());
    REG(WCoef::C_S2_cd, SM, Builtin, std::make_shared<C_S2_cd>());
    REG(WCoef::C_T_cd,  SM, Builtin, std::make_shared<C_T_cd>());

    REG(WCoef::C_V1_cd, SUSY, Builtin, std::make_shared<C_V1_cd_SUSY>());
    REG(WCoef::C_V2_cd, SUSY, Builtin, std::make_shared<C_V2_cd_SUSY>());
    REG(WCoef::C_S1_cd, SUSY, Builtin, std::make_shared<C_S1_cd_SUSY>());
    REG(WCoef::C_S2_cd, SUSY, Builtin, std::make_shared<C_S2_cd_SUSY>());
    REG(WCoef::C_T_cd,  SUSY, Builtin, std::make_shared<C_T_cd_SUSY>());

    REG(WCoef::C_V1_cd, THDM, Builtin, std::make_shared<C_V1_cd_THDM>());
    REG(WCoef::C_V2_cd, THDM, Builtin, std::make_shared<C_V2_cd_THDM>());
    REG(WCoef::C_S1_cd, THDM, Builtin, std::make_shared<C_S1_cd_THDM>());
    REG(WCoef::C_S2_cd, THDM, Builtin, std::make_shared<C_S2_cd_THDM>());
    REG(WCoef::C_T_cd,  THDM, Builtin, std::make_shared<C_T_cd_THDM>());

    for (WCoef c : {WCoef::C_V1_cd, WCoef::C_V2_cd, WCoef::C_S1_cd, WCoef::C_S2_cd, WCoef::C_T_cd})
        REG(c, Model::SM, Backend::Marty, make_marty(ctx, coef));

        //TODO : marty other cases
}

void register_CC_su(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::C_V1_su, SM, Builtin, std::make_shared<C_V1_su>());
    REG(WCoef::C_V2_su, SM, Builtin, std::make_shared<C_V2_su>());
    REG(WCoef::C_S1_su, SM, Builtin, std::make_shared<C_S1_su>());
    REG(WCoef::C_S2_su, SM, Builtin, std::make_shared<C_S2_su>());
    REG(WCoef::C_T_su,  SM, Builtin, std::make_shared<C_T_su>());

    REG(WCoef::C_V1_su, SUSY, Builtin, std::make_shared<C_V1_su_SUSY>());
    REG(WCoef::C_V2_su, SUSY, Builtin, std::make_shared<C_V2_su_SUSY>());
    REG(WCoef::C_S1_su, SUSY, Builtin, std::make_shared<C_S1_su_SUSY>());
    REG(WCoef::C_S2_su, SUSY, Builtin, std::make_shared<C_S2_su_SUSY>());
    REG(WCoef::C_T_su,  SUSY, Builtin, std::make_shared<C_T_su_SUSY>());

    REG(WCoef::C_V1_su, THDM, Builtin, std::make_shared<C_V1_su_THDM>());
    REG(WCoef::C_V2_su, THDM, Builtin, std::make_shared<C_V2_su_THDM>());
    REG(WCoef::C_S1_su, THDM, Builtin, std::make_shared<C_S1_su_THDM>());
    REG(WCoef::C_S2_su, THDM, Builtin, std::make_shared<C_S2_su_THDM>());
    REG(WCoef::C_T_su,  THDM, Builtin, std::make_shared<C_T_su_THDM>());

    for (WCoef c : {WCoef::C_V1_su, WCoef::C_V2_su, WCoef::C_S1_su, WCoef::C_S2_su, WCoef::C_T_su})
        REG(c, Model::SM, Backend::Marty, make_marty(ctx, coef));

        //TODO : marty other cases
}

void register_CC_du(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::C_V1_du, SM, Builtin, std::make_shared<C_V1_du>());
    REG(WCoef::C_V2_du, SM, Builtin, std::make_shared<C_V2_du>());
    REG(WCoef::C_S1_du, SM, Builtin, std::make_shared<C_S1_du>());
    REG(WCoef::C_S2_du, SM, Builtin, std::make_shared<C_S2_du>());
    REG(WCoef::C_T_du,  SM, Builtin, std::make_shared<C_T_du>());

    REG(WCoef::C_V1_du, SUSY, Builtin, std::make_shared<C_V1_du_SUSY>());
    REG(WCoef::C_V2_du, SUSY, Builtin, std::make_shared<C_V2_du_SUSY>());
    REG(WCoef::C_S1_du, SUSY, Builtin, std::make_shared<C_S1_du_SUSY>());
    REG(WCoef::C_S2_du, SUSY, Builtin, std::make_shared<C_S2_du_SUSY>());
    REG(WCoef::C_T_du,  SUSY, Builtin, std::make_shared<C_T_du_SUSY>());

    REG(WCoef::C_V1_du, THDM, Builtin, std::make_shared<C_V1_du_THDM>());
    REG(WCoef::C_V2_du, THDM, Builtin, std::make_shared<C_V2_du_THDM>());
    REG(WCoef::C_S1_du, THDM, Builtin, std::make_shared<C_S1_du_THDM>());
    REG(WCoef::C_S2_du, THDM, Builtin, std::make_shared<C_S2_du_THDM>());
    REG(WCoef::C_T_du,  THDM, Builtin, std::make_shared<C_T_du_THDM>());

    for (WCoef c : {WCoef::C_V1_du, WCoef::C_V2_du, WCoef::C_S1_du, WCoef::C_S2_du, WCoef::C_T_du})
        REG(c, Model::SM, Backend::Marty, make_marty(ctx, coef));

        //TODO : marty other cases
}


void register_MesonMixing(CoefficientRegistry& reg) {
    using enum Backend;using enum Model;
    REG(WCoef::C_BD_1, SM, Builtin, std::make_shared<C_mix_bd_1>());
    REG(WCoef::CT_BD_1, SM, Builtin, std::make_shared<C_mix_bd_1_tilde>());
    REG(WCoef::C_BD_2, SM, Builtin, std::make_shared<C_mix_bd_2>());
    REG(WCoef::CT_BD_2, SM, Builtin, std::make_shared<C_mix_bd_2_tilde>());
    REG(WCoef::C_BD_3,  SM, Builtin, std::make_shared<C_mix_bd_3>());
    REG(WCoef::CT_BD_3,  SM, Builtin, std::make_shared<C_mix_bd_3_tilde>());
    REG(WCoef::C_BD_4,  SM, Builtin, std::make_shared<C_mix_bd_4>());
    REG(WCoef::C_BD_5,  SM, Builtin, std::make_shared<C_mix_bd_5>());

    REG(WCoef::C_BS_1, SM, Builtin, std::make_shared<C_mix_bs_1>());
    REG(WCoef::CT_BS_1, SM, Builtin, std::make_shared<C_mix_bs_1_tilde>());
    REG(WCoef::C_BS_2, SM, Builtin, std::make_shared<C_mix_bs_2>());
    REG(WCoef::CT_BS_2, SM, Builtin, std::make_shared<C_mix_bs_2_tilde>());
    REG(WCoef::C_BS_3,  SM, Builtin, std::make_shared<C_mix_bs_3>());
    REG(WCoef::CT_BS_3,  SM, Builtin, std::make_shared<C_mix_bs_3_tilde>());
    REG(WCoef::C_BS_4,  SM, Builtin, std::make_shared<C_mix_bs_4>());
    REG(WCoef::C_BS_5,  SM, Builtin, std::make_shared<C_mix_bs_5>());

    REG(WCoef::C_SD_1, SM, Builtin, std::make_shared<C_mix_sd_1>());
    REG(WCoef::CT_SD_1, SM, Builtin, std::make_shared<C_mix_sd_1_tilde>());
    REG(WCoef::C_SD_2, SM, Builtin, std::make_shared<C_mix_sd_2>());
    REG(WCoef::CT_SD_2,  SM, Builtin, std::make_shared<C_mix_sd_2_tilde>());
    REG(WCoef::C_SD_3,  SM, Builtin, std::make_shared<C_mix_sd_3>());
    REG(WCoef::CT_SD_3,  SM, Builtin, std::make_shared<C_mix_sd_3_tilde>());
    REG(WCoef::C_SD_4,  SM, Builtin, std::make_shared<C_mix_sd_4>());
    REG(WCoef::C_SD_5, SM, Builtin, std::make_shared<C_mix_sd_5>());

    REG(WCoef::C_CU_1, SM, Builtin, std::make_shared<C_mix_cu_1>());
    REG(WCoef::CT_CU_1, SM, Builtin, std::make_shared<C_mix_cu_1_tilde>());
    REG(WCoef::C_CU_2, SM, Builtin, std::make_shared<C_mix_cu_2>());
    REG(WCoef::CT_CU_2,  SM, Builtin, std::make_shared<C_mix_cu_2_tilde>());
    REG(WCoef::C_CU_3,  SM, Builtin, std::make_shared<C_mix_cu_3>());
    REG(WCoef::CT_CU_3,  SM, Builtin, std::make_shared<C_mix_cu_3_tilde>());
    REG(WCoef::C_CU_4,  SM, Builtin, std::make_shared<C_mix_cu_4>());
    REG(WCoef::C_CU_5, SM, Builtin, std::make_shared<C_mix_cu_5>());
}

void register_K(CoefficientRegistry& reg) {
    using enum Model; using enum Backend;

    REG(WCoef::CK9, SM, Builtin, std::make_shared<CK9>());
    REG(WCoef::CK10, SM, Builtin, std::make_shared<CK10>());
    REG(WCoef::CKQ1, SM, Builtin, std::make_shared<CKQ1>());
    REG(WCoef::CKQ2, SM, Builtin, std::make_shared<CKQ2>());
    REG(WCoef::CPK9,  SM, Builtin, std::make_shared<CPK9>());
    REG(WCoef::CPK10,  SM, Builtin, std::make_shared<CPK10>());
    REG(WCoef::CPKQ1,  SM, Builtin, std::make_shared<CPKQ1>());
    REG(WCoef::CPKQ2,  SM, Builtin, std::make_shared<CPKQ2>());
    REG(WCoef::CK_L,  SM, Builtin, std::make_shared<CK_L>());

    for (WCoef c : {WCoef::CK9, WCoef::CK10, WCoef::CKQ1, WCoef::CKQ2, WCoef::CPK9, WCoef::CPK10, WCoef::CPKQ1, WCoef::CPKQ2, WCoef::CK_L})
        REG(c, Model::SM, Backend::Marty, make_marty(ctx, coef));

        //TODO : marty other cases
}
