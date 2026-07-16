#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <memory>
#include <set>
#include <iostream>

#include "FileNameManager.h"
#include "GeneralNumModelModifier.h"
#include "TemplateManager.h"
#include "IncludeManager.hpp"
#include "LineProcessor.h"
#include "ModelWriter.h"
#include "MartyFileWriter.h"
#include "SMParamSetter.h"
#include "IMartyParameterProxy.h"
#include "IInterpreterPortsFactory.h"
#include "IParameterResolver.h"
#include "ICoreAPI.h"
#include "Extractor.h"
#include "Include.h"

namespace fs = std::filesystem;


class FakeResolver : public IParameterResolver {
public:
    std::unordered_map<std::string, ResolvedParam> table;

    std::unordered_map<std::string, ResolvedParam>
    resolve(const std::vector<Extractor::Parameter>& params,
            bool /*modelIsSM*/) const override
    {
        std::unordered_map<std::string, ResolvedParam> out;
        for (const auto& p : params) {
            auto it = table.find(p.name);
            if (it != table.end()) out[p.name] = it->second;
            else                   out[p.name] = {"MASS", LhaID(25), false, false};
        }
        return out;
    }

    std::unique_ptr<IParameterResolver> clone() const override {
        auto r = std::make_unique<FakeResolver>();
        r->table = table;
        return r;
    }
};

class FakePortsFactory : public IInterpreterPortsFactory {
public:
    std::unordered_map<std::string, ResolvedParam> preset;

    std::unique_ptr<IParameterResolver>
    makeResolver(const std::string&, const std::string&, const std::string&) const override {
        auto r = std::make_unique<FakeResolver>();
        r->table = preset;
        return r;
    }
};

class DummyCoreAPI : public ICoreAPI<Model> {
    Model m;
public:
    explicit DummyCoreAPI(Model mm): m(mm) {}
    Model get() override { return m; }
};

class DummySMProxy : public IMartyParameterProxy<std::string, LhaID> {
public:
    scalar_t operator()(const std::string& blk, const LhaID& id) const override {
        if (blk=="MASS_EW_SCALE" && id==LhaID(5,1)) return 4.7;    // mb(M_W)
        if (blk=="MASS_EW_SCALE" && id==LhaID(6))   return 173.0;  // mt(M_W)
        if (blk=="SMINPUTS"      && id==LhaID(7,1)) return 0.231;  // sin^2(theta_W)
        return 0.0;
    }
};

class DummyBSMProxy : public IMartyParameterProxy<std::string, LhaID> {
public:
    scalar_t operator()(const std::string& blk, const LhaID& id) const override {
        if (blk=="XBLK"   && id==LhaID(1)) return scalar_t(3.0, 4.0);  // 3 + 4i
        if (blk=="MINPAR" && id==LhaID(3)) return 10.0;
        return 0.0;
    }
};


static void write_file(const fs::path& p, const std::string& s) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p); f << s; f.flush();
}
static std::string slurp(const fs::path& p) {
    std::ifstream f(p); std::stringstream ss; ss << f.rdbuf(); return ss.str();
}


int main() {
    std::cout << "== Numeric Pipeline INTEGRATION ==\n";

    const fs::path root   = fs::temp_directory_path() / "marty_num_integ";
    const std::string templ  = (root / "templ").string()  + "/";
    const std::string base   = (root / "base").string()   + "/";
    const std::string assets = (root / "assets").string() + "/";

    fs::create_directories(root / "templ");
    fs::create_directories(root / "base");
    fs::create_directories(root / "assets");

    FileNameManager::setTestingRoots(templ, base, assets);

    const std::string wilson = "C7";
    const std::string model  = "THDM";
    auto mgr = FileNameManager::getInstance(wilson, model);

    write_file(mgr->getBaseHelperFileName("h"),   "// csv_helper.h (base)\n");
    write_file(mgr->getBaseHelperFileName("cpp"), "// csv_helper.cpp (base)\n");

    write_file(mgr->getNumGeneratedFileName(),
        "#include <iostream>\n"
        "using namespace std;\n"
        "int main(){\n"
        "  return 0;\n"
        "}\n"
    );

    write_file(mgr->getNumParamFileName(),
        "csl::InitSanitizer<real_t> mb { \"mb\" };\n"
        "csl::InitSanitizer<complex_t> Yb { \"Yb\" };\n"
    );

    FakePortsFactory ports_a;
    ports_a.preset["mb"] = {"MASS", LhaID(5), false, false};
    ports_a.preset["Yb"] = {"XBLK", LhaID(1), true,  true};

    std::shared_ptr<IInterpreterPortsFactory> ports = std::make_shared<FakePortsFactory>(ports_a);

    auto api  = std::make_shared<DummyCoreAPI>(Model::THDM);
    auto smp  = std::make_shared<DummySMProxy>();
    auto bsmp = std::make_shared<DummyBSMProxy>();
    std::set<std::string> specials{};

    std::unique_ptr<SMParamSetter> setter = std::make_unique<SMParamSetter>(model, specials, smp, bsmp);

    auto numMod = std::make_unique<GeneralNumModelModifier>(wilson, model, std::move(setter), api, ports, /*force=*/false);

    NumericTemplateManager tm(mgr->getLibDir());
    tm.setModelAndWilson(model, wilson);
    tm.setNumModelModifier(std::move(numMod));

    const std::string outPath = mgr->getNumGeneratedFileName();
    tm.generateTemplate("ignored", outPath);

    assert(std::filesystem::exists(mgr->getHelperFileName("h")));
    assert(std::filesystem::exists(mgr->getHelperFileName("cpp")));

    const auto code = slurp(outPath);
    assert(code.find("#include \"csv_helper.h\"") != std::string::npos);
    assert(code.find("int main(int argc, char** argv)") != std::string::npos);
    assert(code.find("param_t param;") != std::string::npos);

    const fs::path params_csv = mgr->getParamFileName();
    assert(std::filesystem::exists(params_csv));
    const auto csv = slurp(params_csv);

    assert(csv.find("mb,4.7")     != std::string::npos);
    assert(csv.find("Yb_rel,3")   != std::string::npos);
    assert(csv.find("Yb_img,4")   != std::string::npos);

    auto deps = tm.get_dependencies();
    assert(!deps.empty());

    FileNameManager::clearTestingRoots();

    std::cout << "INTEGRATION OK\n";
    return 0;
}
