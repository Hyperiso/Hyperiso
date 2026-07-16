#include "marty.h"
#include "../../../Third_party/MARTY/MARTY_INSTALL/include/marty/models/sm.h"
#include <vector>

using namespace mty::sm_input;
using namespace csl;
using namespace std;

namespace mty {

class ZPrime_Model : public mty::SM_Model {

  public:
    ZPrime_Model(bool initialize = true);

    void init();

    void initGauge();
    void initContent();
    void initHiggsPotential();
    void initFermions();
    void initYukawas();
    void getToLowEnergyLagrangian();
    void gaugeSymmetryBreaking();
    void HiggsVEVExpansion();
    void flavorSymmetryBreaking();
    void rotateNeutrinos(); 
    void adjust();

    friend std::ostream &operator<<(std::ostream     &out,
                                    ZPrime_Model const &model);

  protected:
    csl::Expr m_X;
    csl::Expr vphi;
    //mixing constant cam
    csl::Expr lambda_mix; 
};

ZPrime_Model::ZPrime_Model(bool initialize) : SM_Model(false)
{
    if (initialize)
        init();
        refresh();
}

void ZPrime_Model::initHiggsPotential()
{
    Particle H = scalarboson_s("H", *this);
    H->setGroupRep("L", 1);
    H->setGroupRep("Y", {1, 2});
    addParticle(H);

    csl::Index a = generateIndex("L", H);
    //csl::Index dirac_i = DiracIndex();
    //csl::Index f1 = FlavorIndex(*this, "f");
    //csl::Index mu = MinkowskiIndex();

    csl::Expr mh  = sm_input::m_h;
    csl::Expr H2  = csl::GetComplexConjugate(H(a)) * H(a);
    csl::Expr m2  = mh * mh / 2;
    csl::Expr lam = mh * mh / (2 * v * v);

    // Mexican hat potential
    // addLagrangianTerm(m2 / (2 * v * v) * pow_s(H2 - v * v / 2, 2));  Same as below
    addLagrangianTerm(m2 * H2);
    addLagrangianTerm(-lam * csl::pow_s(H2, 2));
    // later on: m   = m_h / sqrt(2)
    //           lam = mh^2 / (2*v^2)
    //           (With H0 -> (v + h0) / sqrt(2))

    Particle Phi = scalarboson_s("Phi", *this);
    std::cout << "higgstest1" << std::endl;
    // Phi->setGroupRep("L", 1);
    Phi->setGroupRep("X", 2);
    std::cout << "higgstest2" << std::endl;
    addParticle(Phi);
    std::cout << "higgstest3" << std::endl;

    m_X = constant_s("m_X");
    vphi = constant_s("v_phi");
    csl::Expr mphi  = m_X;
    std::cout << "higgstest4" << std::endl;
    csl::Expr Phi2  = csl::GetComplexConjugate(Phi()) * Phi();
    std::cout << "higgstest5" << std::endl;
    csl::Expr lamphi = mphi * mphi / (2 * vphi * vphi);
    std::cout << "higgstest6" << std::endl;
    // Mexican hat potential
    //why isn't there a lamphi in the first added term? cam
    addLagrangianTerm(mphi * mphi * Phi2);
    addLagrangianTerm(-lamphi * csl::pow_s(Phi2, 2));

    // Mixing terms cam
    lambda_mix = constant_s("lambda_mix");
    addLagrangianTerm(lambda_mix * H2 *Phi2);
    addLagrangianTerm(lambda_mix * H2 * mphi * mphi);
    addLagrangianTerm(lambda_mix * mphi * mphi *Phi2);
}

void ZPrime_Model::init()
{
    initContent();
    std::cout << "1" << std::endl;
    getToLowEnergyLagrangian();
    std::cout << "2" << std::endl;
}

void ZPrime_Model::getToLowEnergyLagrangian()
{
    ZPrime_Model::gaugeSymmetryBreaking();
    ZPrime_Model::HiggsVEVExpansion();
    SM_Model::diagonalizeSMMassMatrices();
    SM_Model::replaceLeptonYukawa();
    SM_Model::replaceUpYukawa();
    SM_Model::replaceDownYukawa();
    SM_Model::flavorSymmetryBreaking();
    ZPrime_Model::adjust();
}

void ZPrime_Model::initContent()
{
    std::cout << "test1" << std::endl;
    ZPrime_Model::initGauge();
    std::cout << "test2" << std::endl;
    ZPrime_Model::initFermions();
    std::cout << "test3" << std::endl;
    ZPrime_Model::initHiggsPotential();
    std::cout << "test4" << std::endl;
    ZPrime_Model::initYukawas();
}


void ZPrime_Model::initGauge()
{
    addGaugedGroup(group::Type::SU, "C", 3, g_s);
    addGaugedGroup(group::Type::SU, "L", 2, csl::constant_s("g_L"));
    addGaugedGroup(group::Type::U1, "Y", csl::constant_s("g_Y"));
    addGaugedGroup(group::Type::U1, "X", csl::constant_s("g_X"));
    addFlavorGroup("SM_flavor", 3);
    Model::init();
    renameParticle("A_Y", "B");
    renameParticle("A_L", "W");
    renameParticle("A_C", "G");
    renameParticle("A_X", "B_X");
    getParticle("G")->setDrawType(drawer::ParticleType::Gluon);
}

void ZPrime_Model::initFermions()
{
    Particle Q = weylfermion_s("Q", *this, Chirality::Left);
    Q->setGroupRep("C", {1, 0});
    Q->setGroupRep("L", 1);
    Q->setGroupRep("Y", {1, 6});
    Q->setGroupRep("X", {1,3});
    Q->setFundamentalFlavorRep("SM_flavor");

    Particle U = weylfermion_s("U_R", *this, Chirality::Right);
    U->setGroupRep("C", {1, 0});
    U->setGroupRep("Y", {2, 3});
    U->setGroupRep("X", {1,3});
    U->setFundamentalFlavorRep("SM_flavor");

    Particle D = weylfermion_s("D_R", *this, Chirality::Right);
    D->setGroupRep("C", {1, 0});
    D->setGroupRep("Y", {-1, 3});
    D->setGroupRep("X", {1,3});
    D->setFundamentalFlavorRep("SM_flavor");

    Particle L = weylfermion_s("L", *this, Chirality::Left);
    L->setGroupRep("L", 1);
    L->setGroupRep("Y", {-1, 2});
    L->setGroupRep("X", -1);
    L->setFundamentalFlavorRep("SM_flavor");

    Particle E = weylfermion_s("E_R", *this, Chirality::Right);
    E->setGroupRep("Y", -1);
    E->setGroupRep("X", -1);
    E->setFundamentalFlavorRep("SM_flavor");

    Particle N = weylfermion_s("nu_R", *this, Chirality::Right);
    N->setGroupRep("X", -1);
    // Q->setGroupRep("X", -1);
    N->setFundamentalFlavorRep("SM_flavor");

    addParticles({Q, U, D, L, E, N});
}

void ZPrime_Model::initYukawas()
{
    auto  *flavorSpace = getVectorSpace("SM_flavor");
    Tensor Yu("Yu", {flavorSpace, flavorSpace});
    Yu->setComplexProperty(ComplexProperty::Complex);
    Tensor Yd("Yd", {flavorSpace, flavorSpace});
    Yd->setComplexProperty(ComplexProperty::Complex);
    Tensor Ye("Ye", {flavorSpace, flavorSpace});
    Ye->setComplexProperty(ComplexProperty::Complex);
    Tensor Ynu("Ynu", {flavorSpace, flavorSpace});
    Ynu->setComplexProperty(ComplexProperty::Complex);
    Tensor Yx("Yx", {flavorSpace, flavorSpace});
    Yx->setComplexProperty(ComplexProperty::Complex);
    Tensor eps = getVectorSpace("L", "Q")->getEpsilon();
    Index  I   = flavorSpace->generateIndex();
    Index  J   = flavorSpace->generateIndex();
    Index  a   = generateIndex("C", "Q");
    Index  i   = generateIndex("L", "Q");
    Index  j   = generateIndex("L", "Q");
    Index  al  = DiracIndex();

    Particle Q = getParticle("Q");
    Particle U = getParticle("U_R");
    Particle D = getParticle("D_R");
    Particle L = getParticle("L");
    Particle E = getParticle("E_R");
    Particle N = getParticle("nu_R");
    Particle H = getParticle("H");
    Particle Phi = getParticle("Phi");

    addLagrangianTerm(Yu({I, J}) * GetComplexConjugate(H(i)) * eps({i, j})
                          * GetComplexConjugate(Q({I, a, j, al}))
                          * U({J, a, al}),
                      true);
    addLagrangianTerm(-Yd({I, J}) * H(i)
                          * GetComplexConjugate(Q({I, a, i, al}))
                          * D({J, a, al}),
                      true);
    addLagrangianTerm(-Ye({I, J}) * H(i) * GetComplexConjugate(L({I, i, al}))
                          * E({J, al}),
                      true);

    addLagrangianTerm(-Ynu({I, J}) * GetComplexConjugate(H(i)) * eps({i, j}) * GetComplexConjugate(L({I,j, al}))
                          * N({J, al}),
                      true);

    // addLagrangianTerm(-Yx({I, J}) * Phi * GetComplexConjugate(N({I, al}))
    //                       * N({J, al}),
    //                   true);

    addTensorCoupling(Ye);
    addTensorCoupling(Yu);
    addTensorCoupling(Yd);
}

void ZPrime_Model::HiggsVEVExpansion()
{
    ///////////////////////////////////////////////////
    // Actual gauge (spontaneous) symmetry breaking
    ///////////////////////////////////////////////////

    csl::Expr v = sm_input::v;

    Particle H1 = getParticle("H_1");
    Particle H2 = getParticle("H_2");

    Particle h0 = scalarboson_s("h", *this); // SM Higgs boson
    Particle Gp = scalarboson_s("Gp ; G^+", *this);
    Particle G0 = scalarboson_s("G0 ; G^0", *this);
    h0->setSelfConjugate(true);
    G0->setSelfConjugate(true);

    replace(H1, Gp());
    replace(H2, (h0() + CSL_I * G0() + v) / csl::sqrt_s(2));

    Particle phi0 = scalarboson_s("phi", *this);
    phi0->setSelfConjugate(true);

    Particle Phi = getParticle("Phi");
    replace(Phi, (vphi+phi0()) / csl::sqrt_s(2));
}

void ZPrime_Model::gaugeSymmetryBreaking()
{
    ///////////////////////////////////////////////////
    // Breaking gauge SU(2)_L symmetry, renaming
    ///////////////////////////////////////////////////

    BreakGaugeSymmetry(*this, "Y");
    BreakGaugeSymmetry(*this, "L");
    BreakGaugeSymmetry(*this, "X");
    renameParticle("Q_1", "U_L");
    renameParticle("Q_2", "D_L");
    renameParticle("L_1", "Nu_L");
    renameParticle("L_2", "E_L");
    renameParticle("B_X", "Z_X");
    ///////////////////////////////////////////////////
    // Replacements to get SM particles W +-
    ///////////////////////////////////////////////////

    Particle W1   = GetParticle(*this, "W_1");
    Particle W2   = GetParticle(*this, "W_2");
    Particle W_SM = W1->generateSimilar("W");
    W_SM->setSelfConjugate(false);

    Particle cW1 = getParticle("c_W_1");
    Particle cW2 = getParticle("c_W_2");
    Particle cWp = W_SM->getGhostBoson();
    cWp->setName("c_Wp ; c_{+}");
    Particle cWm = ghostboson_s("c_Wm; c_{-}", W_SM, true);
    W_SM->setConjugatedGhostBoson(cWm);

    csl::Index mu    = MinkowskiIndex();
    csl::Index nu    = MinkowskiIndex();
    csl::Expr  W_p   = W_SM(+mu);
    csl::Expr  W_m   = csl::GetComplexConjugate(W_SM(+mu));
    csl::Expr  F_W_p = W_SM({+mu, +nu});
    csl::Expr  F_W_m = csl::GetComplexConjugate(W_SM({+mu, +nu}));

    auto W1_expr = [](csl::Expr const &Wp, csl::Expr const &Wm) {
        return (Wp + Wm) / csl::sqrt_s(2);
    };
    auto W2_expr = [](csl::Expr const &Wp, csl::Expr const &Wm) {
        return CSL_I * (Wp - Wm) / csl::sqrt_s(2);
    };
    replace(W1, W1_expr(W_p, W_m));
    replace(W2, W2_expr(W_p, W_m));
    replace(W1->getFieldStrength(), W1_expr(F_W_p, F_W_m));
    replace(W2->getFieldStrength(), W2_expr(F_W_p, F_W_m));
    replace(cW1, W1_expr(cWp, cWm));
    replace(cW2, W2_expr(cWp, cWm));
}

void ZPrime_Model::adjust()
{
    using namespace sm_input;
    replace(v, (2 * M_W * csl::sin_s(theta_W)) / e_em);
    replace(getParticle("W")->getMass(), M_W);
    getParticle("W")->setMass(M_W);
    replace(getParticle("Z")->getMass(), M_Z);
    getParticle("Z")->setMass(M_Z);
    promoteToGoldstone("Gp", "W");
    promoteToGoldstone("G0", "Z");
    replace(getParticle("Z_X")->getMass(), m_X);
    getParticle("Z_X")->setMass(m_X);
    addGaugeFixingTerms();
}

std::ostream &operator<<(std::ostream &out, ZPrime_Model const &model)
{
    return out << *static_cast<Model const *>(&model);
}

}