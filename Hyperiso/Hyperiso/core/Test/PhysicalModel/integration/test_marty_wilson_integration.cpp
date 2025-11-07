// testMartyWilsonIntegration.cpp
#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <map>

#include "MartyWilson.h"
#include "InterpretedParam.h"
#include "WilsonGroup.h"
#include "IBlockComposer.h"
#include "IParameterProxy.h"
#include "ICoreAPI.h"
#include "Parameter.h"
#include "Include.h"

namespace fs = std::filesystem;

// Spy IParameterProxy (utilisé par get_matching_coefficient / running si besoin)
class SpyProxy : public IParameterProxy<std::string, LhaID> {
public:
    mutable std::string last_block;
    mutable LhaID      last_id{0};
    double ret = 0.0;
    scalar_t operator()(const std::string& b, const LhaID& id) const override {
        last_block = b; last_id = id; return ret;
    }
    bool exist(const std::string&, const LhaID&) const override { return true; }
    double get_scale(const std::string&) const override { return 0.0; }
};

// Spy Block composer
class SpyComposer : public IBlockComposer {
public:
    struct Rec { ParamId target; std::unordered_set<ParamId> sources; };
    std::vector<Rec> calls;

    void compose_block(const std::string&,
                       const std::unordered_map<ParameterType, std::vector<std::string>>&,
                       const DepUpdateFunc&) override {}

    void compose_parameter(const ParamId& id,
                           const std::unordered_set<ParamId>& src,
                           const DepParamUpdateFunc&) override {
        calls.push_back({id, src});
    }

    void remove_block(const std::string&) override {}
    void update(const std::string&) override {}
    void remove_all_composed_blocks() override {}
};

// Dummies CoreAPI
class DummyBoolAPI    : public ICoreAPI<bool>        { public: bool        get() override { return false; } };
class DummyStringAPI  : public ICoreAPI<std::string> { public: std::string get() override { return "SM";  } };
class DummyPathAPI    : public ICoreAPI<fs::path>    { public: fs::path    get() override { return fs::path("/dev/null"); } };

// Fake Marty proxy (écrit le CSV)
class FakeMartyProxy : public IMartyWilsonProxy<InterpretedParam> {
public:
    fs::path csv_path;
    explicit FakeMartyProxy(const fs::path& p): csv_path(p) {}
    void calculate(std::string wilson, std::string /*model*/, double Q_match,
                   std::string /*model_path*/, bool /*new_params*/ = false) override {
        fs::create_directories(csv_path.parent_path());
        std::ofstream out(csv_path);
        out << "Q_match," << wilson << "_real," << wilson << "_img\n";
        // Choix deterministic (C7: -0.2 + 3e-7 i)
        double re = -2.0e-1, im = 3.0e-7;
        out << Q_match << "," << std::scientific << re << "," << std::scientific << im << "\n";
        out.flush();
    }
    std::set<std::string> get_special_blocks() override { return {}; }
    std::unordered_set<InterpretedParam> get_dependencies(std::string) override {
        return {
            InterpretedParam{"SMINPUTS", LhaID(7,1), false, false},
            InterpretedParam{"XBLK",     LhaID(1),   true,  true }
        };
    }
};

// Groupe concret minimal pour fixer l’id et initialiser
class TestGroup : public CoefficientGroup {
public:
    using CoefficientGroup::CoefficientGroup;
    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<TestGroup>(*this);
    }
    void set_id(WGroupId gid) { this->id = gid; }
    void set_type(ContributionType t) { this->wilson_type = t; }
};

int main() {
    std::cout << "== MartyWilson INTEGRATION ==\n";

    // Ports/adapters
    auto proxy   = std::make_shared<SpyProxy>();
    auto comp    = std::make_shared<SpyComposer>();
    auto use_mty = std::make_shared<DummyBoolAPI>();
    auto mname   = std::make_shared<DummyStringAPI>();
    auto mpath   = std::make_shared<DummyPathAPI>();
    WilsonGroupAdapterConfig cfg(proxy, comp, use_mty, mname, mpath);

    // Fake Marty proxy (écrit CSV)
    const fs::path tmpcsv = fs::temp_directory_path() / "mw_integ" / "SM_wilson.csv";
    auto mw_proxy = std::make_shared<FakeMartyProxy>(tmpcsv);

    // Construit un MartyWilson(C7)
    MartyWilsonConfig c7cfg(
        "SM",                          // model name
        LhaID(305,4422,0,0),           // C7 LO id
        GroupMapper::str(WGroup::B, ScaleType::MATCHING),
        "/dev/null",
        mw_proxy
    );
    c7cfg.csv_path = tmpcsv.string();

    auto c7 = std::make_shared<MartyWilson>(c7cfg);

    // On fabrique un groupe B et on y met C7
    std::map<std::string, std::shared_ptr<WilsonCoefficient>> coeffs;
    coeffs["C7"] = c7;

    TestGroup grp(cfg);
    grp.set_id(GroupMapper::to_id(WGroup::B));
    grp.set_type(ContributionType::SM);
    // Injecte puis init en LO
    grp.insert(coeffs.begin(), coeffs.end());
    grp.init(QCDOrder::LO);

    // Au moins 1 appel à compose_parameter (pour C7 LO)
    assert(!comp->calls.empty());

    // On s’attend à une cible = (B_MATCH, LHAID(C7 LO))
    const auto match_blk = GroupMapper::str(WGroup::B, ScaleType::MATCHING);
    bool saw_c7 = false;
    for (const auto& r : comp->calls) {
        if (r.target == ParamId{match_blk, c7->get_lhaid(QCDOrder::LO)}) {
            // Sources incluent EW_SCALE + nos deps
            assert(r.sources.count(ParamId{ParameterType::WILSON, "EW_SCALE", LhaID(1)}) == 1);
            assert(r.sources.count(ParamId{ParameterType::SM,     "SMINPUTS", LhaID(7,1)}) == 1);
            assert(r.sources.count(ParamId{ParameterType::BSM,    "XBLK",     LhaID(1)})   == 1);
            saw_c7 = true;
        }
    }
    assert(saw_c7);

    // Exerce la fonction LO de C7 au Q=81 → le Fake écrit la valeur et on lit
    auto fLO = c7->get_func(QCDOrder::LO);
    ParamId pid{ParameterType::WILSON, "EW_SCALE", 1};
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> src{
        { pid, std::make_shared<Parameter>(pid, 81.0, 0.0, 0.0) }
    };
    scalar_t v = fLO(ParamSrc(src));
    assert(std::abs(v.real() - (-2.0e-1)) < 1e-12);
    assert(std::abs(v.imag() - ( 3.0e-7)) < 1e-12);

    std::cout << "✅ INTEGRATION OK\n";
    return 0;
}
