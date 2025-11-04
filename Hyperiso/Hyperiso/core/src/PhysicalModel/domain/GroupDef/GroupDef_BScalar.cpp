#include "GroupDefinition.h"
#include "BScalarWilsonGroup.h"
#include "Math_SUSY.h"

using CGS = CoefficientGroupSources;

static std::unordered_map<WCoef, scalar_t>
BScalar_SUSY_Base1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const BlockSrc& src)
{
    const auto ids = WCoefMapper::get_group(WGroup::BScalar);

    const auto itLO = coef_matching.find(QCDOrder::LO);
    const auto* matchLO = (itLO != coef_matching.end()) ? &itLO->second : nullptr;

    auto getM = [&](WCoef c) -> scalar_t {
        if (!matchLO) return scalar_t(0);
        auto it = matchLO->find(c);
        return (it != matchLO->end()) ? it->second : scalar_t(0);
    };

    const double eta    = src.get_val("WPARAM_RUN_SM",2);
    const double beta_0 = src.get_val("WPARAM_SI_SM",5); // TODO: passer aux vrais params QCD si dispo

    const double fact = std::pow(eta, -4.0 / beta_0);

    std::unordered_map<WCoef, scalar_t> out;
    out.reserve(ids.size());
    for (auto c : ids) {
        out[c] = fact * getM(c);
    }

    complex_t coeff_temp2 = 0;
    if(src.get_val("MASS",46)!=0.||src.get_val("MASS",45)!=0.) {
            double mass_b_2 = src.get_val("QCD",LhaID(5, 2));
            double mA = src.get_val("MASS",36);
			if(mA < mass_b_2) {	
            double lambdaNMSSM = 1;
            double lambdaSNMSSM = 1;
            double AlambdaNSSM = 1;
            double kappaNMSSM = 1;
            double m_Bs = 1;
            double mass_nutl = 1;
            
            double sw2 = src.get_val("WPARAM_SI_SM",4);
            double mH = src.get_val("MASS",37);
		    
            double mass_top_muW = src.get_val("WPARAM_MATCH_SM",6);
            double g2 = src.get_val("GAUGE",2);
            double tanb = src.get_val("HMIX",2);
            double mW = src.get_val("MASS",24);

            double mH0[4],mA0[3],mstop[3];
        
            mstop[0]=src.get_val("MASS",2000002); //mass upr, is that right ?
            mstop[1]=src.get_val("MASS",1000006);
            mstop[2]=src.get_val("MASS",2000006);

            complex_t CAH={0,-lambdaNMSSM*AlambdaNSSM/g2/mW*tanb*f30(mH*mH/mass_top_muW/mass_top_muW,mW*mW/mass_top_muW/mass_top_muW)};
            complex_t CAc{};
            double s=lambdaSNMSSM/lambdaNMSSM;
            double v=sqrt(1./sqrt(2.)/src.get_val("SMINPUTS",2));
            double v_deltam_s=v/s*(sqrt(2.)*AlambdaNSSM-2.*kappaNMSSM*s)/(sqrt(2.)*AlambdaNSSM+kappaNMSSM*s);

            double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
            double T2[4][4][4];
            std::array<std::array<double,4>,4> TU;
            double vu=sqrt(pow(sin(atan(tanb)),2.)/sqrt(2.)/src.get_val("SMINPUTS",2));
            double vd=vu/tanb;

            TU[1][1]=1.;
            for(int ie=0;ie<2;ie++){
                for(int je=0;je<2;je++) {
                    TU[ie+1][je+1]=src.get_val("STOPMIX",{ie+1, je+1});
                }
            }

            for(int je=0;je<2;je++) {
                for(int le=0;le<2;le++) {
                    for(int ae=0;ae<3;ae++) {
                        if (ae <3 ){
                            Ralj[ae][le][je]=-g2/sqrt(2.)*(src.get_val("NMAMIX",{ae+1, 1+1})*src.get_val("UMIX",{2+1, le+1})*src.get_val("VMIX",{2+1, je+1})+src.get_val("NMAMIX",{ae+1, 2+1})*src.get_val("UMIX",{1+1, le+1})*src.get_val("VMIX",{2+1, je+1}))-lambdaNMSSM/sqrt(2.)*src.get_val("NMAMIX",{ae+1, 3+1})*src.get_val("UMIX",{2+1, le+1})*src.get_val("VMIX",{2+1, je+1});
                        }
                        Qalj[ae][le][je]=g2/sqrt(2.)*(src.get_val("NMHMIX",{ae+1,1+1})*src.get_val("UMIX",{2+1, le+1})*src.get_val("VMIX",{2+1, je+1})+src.get_val("NMHMIX",{ae+1, 2+1})*src.get_val("UMIX",{1+1, le+1})*src.get_val("VMIX",{2+1, je+1}))-lambdaNMSSM/sqrt(2.)*src.get_val("NMHMIX",{ae+1, 3+1})*src.get_val("UMIX",{2+1, le+1})*src.get_val("VMIX",{2+1, je+1});
                        for(int ke=1;ke<=3;ke++) {
                            G1[ae][ke][je][le]=(TU[ae][2]*TU[ke][2]-kron(ae,1)*kron(ke,1))*src.get_val("VMIX",{1+1, le+1})*src.get_val("UMIX",{2+1, je+1})-mass_top_muW/sqrt(2.)/sin(atan(tanb))/mW*TU[ae][3]*TU[ke][2]*src.get_val("VMIX",{2+1, le+1})*src.get_val("UMIX",{2+1, je+1});
                        }
                    }
                }
            }
            for(int ae=0;ae<3;ae++) {
                for(int je=0;je<2;je++) {
                    for(int le=0;le<2;le++) {
                        CAc = complex_t(CAc.real(), CAc.imag()+(tanb)/sqrt(2.)*G1[ae][ae][je][le]*(v_deltam_s*kron(le,je)*fabs(src.get_val("WPARAM_SI_BSM",je)/mW)*f80(pow(mstop[ae-1]/src.get_val("WPARAM_SI_BSM",je),2.))-(Ralj[1][je][le]*fabs(src.get_val("WPARAM_SI_BSM",je)/src.get_val("WPARAM_SI_BSM",le))*f30(pow(mstop[ae-1]/src.get_val("WPARAM_SI_BSM",le),2.),pow(src.get_val("WPARAM_SI_BSM",je)/src.get_val("WPARAM_SI_BSM",le),2.))-Ralj[1][le][je]*f40(pow(mstop[ae-1]/src.get_val("WPARAM_SI_BSM",le),2.),pow(src.get_val("WPARAM_SI_BSM",je)/src.get_val("WPARAM_SI_BSM",le),2.)))));
                    }
                }
            }
            complex_t CA=CAH+CAc;
            double width_A0=1.e-6;
            coeff_temp2+=complex_t{v_deltam_s/2.*mass_b_2/sw2*src.get_val("WPARAM_SI_SM",3)*CA/(m_Bs*m_Bs-mA*mA,mA*width_A0)};
        }
    }
    out[WCoef::CQ2] += coeff_temp2;
    return out;
}

static void Setup_BScalar_SUSY_Base1_LO(const BuildContext& ctx, CoefficientGroup& grp) {
    std::map<QCDOrder, CGS> m;

    CGS lo;
    lo.sources = {
        { ParameterType::WILSON, { GroupMapper::str(ctx.group_id, ScaleType::MATCHING),
                                   "WPARAM_RUN_SM", "WPARAM_SI_SM", "WPARAM_MATCH_SM" } },
        { ParameterType::SM,     { "MASS", "SMINPUTS" } },
        { ParameterType::BSM,    { "GAUGE", "HMIX", "STOPMIX", "UMIX", "VMIX", "NMAMIX", "NMHMIX", "MASS" } }
    };
    lo.func = &BScalar_SUSY_Base1_LO_calculation;
    m[QCDOrder::LO] = lo;
    
    CGS nlo;
    nlo.sources = {
        { ParameterType::WILSON, { GroupMapper::str(ctx.group_id, ScaleType::MATCHING),
                                   "WPARAM_RUN_SM", "WPARAM_SI_SM" } }
    };
    nlo.func = &BScalarCoefficientGroup::base_1_NLO_calculation;
    m[QCDOrder::NLO] = nlo;

    grp.add_sources(WilsonBasis::B_STANDARD, std::move(m));


}

namespace GroupDefinitions {
    const GroupDefinition& BScalar() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::BScalar);
            d.members = WCoefMapper::get_group(WGroup::BScalar);

            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER, "WPARAM_RUN_SM", "WPARAM_SI_SM" } }
            };
            lo.func = &BScalarCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            CGS nlo = lo; nlo.func = &BScalarCoefficientGroup::base_1_NLO_calculation;
            m[QCDOrder::NLO] = nlo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));

            d.setup[Model::SUSY].push_back(&Setup_BScalar_SUSY_Base1_LO);
            return d;
        }();
        return def;
    }
}
