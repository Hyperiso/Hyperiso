#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "MartyWilson.h"
#include "InterpretedParam.h"
#include "IParameterProxy.h"
#include "IMartyWilsonProxy.h"   
#include "Parameter.h"   
#include "Include.h"  
#include "CSVReader.h"

namespace fs = std::filesystem;

static void ensure_parent(const fs::path& p) {
    fs::create_directories(p.parent_path());
}
static std::string slurp(const fs::path& p){
    std::ifstream f(p); std::stringstream ss; ss << f.rdbuf(); return ss.str();
}

class FakeMartyProxy : public IMartyWilsonProxy<InterpretedParam> {
public:
    fs::path csv_path;
    std::unordered_map<std::string, std::pair<double,double>> values{
        {"C7",  { -1.942999e-01, 3.438710e-07 }},
        {"C10", { -1.936392e+00, -3.200705e-02 }},
    };

    mutable std::string last_wilson, last_model;
    mutable double last_Q = 0.0;
    mutable fs::path last_model_path;

    explicit FakeMartyProxy(const fs::path& csv) : csv_path(csv) {}

    void calculate(std::string wilson, std::string model, double Q_match,
                   std::string model_path, bool /*new_params*/ = false) override
    {
        last_wilson = wilson; last_model = model; last_Q = Q_match; last_model_path = model_path;

        ensure_parent(csv_path);
        std::ofstream out(csv_path);

        out << "Q_match," << wilson << "_real," << wilson << "_img\n";
        auto it = values.find(wilson);
        double re = (it!=values.end() ? it->second.first : 1.23);
        double im = (it!=values.end() ? it->second.second: -4.56);

        out << Q_match << ",  " << std::scientific << re << " , " << std::scientific << im << "\n";
        out.flush();
    }

    std::set<std::string> get_special_blocks() override {
        return {};
    }

    std::unordered_set<InterpretedParam> get_dependencies(std::string wilson) override {

        std::unordered_set<InterpretedParam> s;
        s.insert(InterpretedParam{"SMINPUTS", LhaID(7,1), false, false});
        s.insert(InterpretedParam{"XBLK",     LhaID(1),   true,  true});
        return s;
    }
};

static std::shared_ptr<Parameter> P(double x) {
    return std::make_shared<Parameter>(ParamId("DUMMY", 0), x, 0.0, 0.0);
}

int main() {
    std::cout << "== MartyWilson UNIT ==\n";

    const fs::path tmp = fs::temp_directory_path() / "mw_unit" / "SM_wilson.csv";

    auto mw_proxy = std::make_shared<FakeMartyProxy>(tmp);

    MartyWilsonConfig cfg(
        /*model_name*/ "SM",
        /*coeff_id*/   LhaID(305,4422,0,0),
        /*storage*/    GroupMapper::str(WGroup::B, ScaleType::MATCHING),
        /*model_path*/ "/dev/null",
        /*proxy*/      mw_proxy
    );
    cfg.csv_path = tmp.string();

    MartyWilson mw(cfg);

    auto fLO = mw.get_func(QCDOrder::LO);

    ParamId pid{ParameterType::WILSON, "EW_SCALE", 1};
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> src{
        { pid, P(81.0) } 
    };

    scalar_t val = fLO(ParamSrc(src));
    // (-1.942999e-01, 3.438710e-07)
    assert(std::abs(val.real() - (-1.942999e-01)) < 1e-12);
    assert(std::abs(val.imag() - ( 3.438710e-07)) < 1e-12);

    auto S = mw.get_sources(QCDOrder::LO);
    assert(S.count(ParamId{ParameterType::WILSON, "EW_SCALE", LhaID(1)}) == 1);
    assert(S.count(ParamId{ParameterType::SM,     "SMINPUTS", LhaID(7,1)}) == 1);
    assert(S.count(ParamId{ParameterType::BSM,    "XBLK",     LhaID(1)})   == 1);

    assert(mw_proxy->last_wilson == "C7");
    assert(mw_proxy->last_model  == "SM");
    assert(std::abs(mw_proxy->last_Q - 81.0) < 1e-12);

    std::cout << " UNIT OK\n";
    return 0;
}
