#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <iostream>
#include "GeneralModelModifier.h"
#include "MartyRuntimeConfig.h"

namespace fs = std::filesystem;

static void write_file(const fs::path& p, const std::string& content) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p);
    f << content;
    f.flush();
}

static std::string slurp(const fs::path& p) {
    std::ifstream f(p);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

int main() {
    std::cout << "== GeneralModelModifier UNIT ==\n";

    const fs::path root = fs::temp_directory_path() / "gmm_unit";
    const fs::path thdm_hdr = root / "models" / "thdm.hpp";
    const fs::path zprime_hdr = root / "models" / "zprime.hpp";
    const fs::path marty_install = root / "marty_install";

    // GeneralModelModifier validates the MARTY runtime path eagerly.  These
    // tests exercise source rewriting only, so provide a minimal fake install
    // tree instead of requiring the full third-party dependency in standard CI.
    write_file(marty_install / "include" / "marty.h", "// test stub\n");
    write_file(marty_install / "lib" / "libmarty.a", "");
    const auto marty = MartyRuntimeConfig::set_external_install_path(marty_install);
    assert(marty.valid);

    write_file(thdm_hdr,
        "template <int N>\n"
        "class THDM_Model : public mty::Model {};\n"
    );

    write_file(zprime_hdr,
        "class ZPrime_Model : public mty::Model {};\n"
    );

    {
        GeneralModelModifier mod(/*wilson*/"none", /*model*/"THDM", thdm_hdr.string(), 2);
        std::string line = "SM_Model sm;";
        mod.modifyLine(line);
        assert(line.find("THDM_Model<2> sm;") != std::string::npos);
    }

    {
        bool threw = false;
        try {
            GeneralModelModifier mod("none", "THDM", thdm_hdr.string());
            (void)mod;
        } catch (const std::runtime_error&) {
            threw = true;
        }
        assert(threw && "Templated MARTY models require an explicit template index");
    }

    {
        GeneralModelModifier mod("none", "ZPrime", zprime_hdr.string());
        std::string line = "SM_Model sm;";
        mod.modifyLine(line);
        assert(line.find("ZPrime_Model sm;") != std::string::npos);
    }

    {
        GeneralModelModifier mod("none", "ZPrime", zprime_hdr.string());
        std::string line = "int flag_SM = 0;";
        mod.modifyLine(line);
        // remplace le premier “SM” rencontré → flag_ZPrime
        assert(line.find("flag_ZPrime = 0;") != std::string::npos);
    }

    {
        GeneralModelModifier mod("none", "ZPrime", zprime_hdr.string());
        std::ostringstream oss;
        fs::path out = root / "out.hpp";
        {
            std::ofstream f(out);
            std::string l = "#include <iostream>";
            mod.addLine(f, l);
            mod.addLine(f, "// nothing");
        }
        std::ifstream ifs(out);
        std::stringstream buf; buf << ifs.rdbuf();
        std::string s = buf.str();

        assert(s.find("#include <iostream>") != std::string::npos);
        assert(s.find(zprime_hdr.string()) != std::string::npos);
        assert(s.find("marty.h") != std::string::npos);
        assert(s.find("// nothing") != std::string::npos);
    }

    {
        // Generic coefficients must keep their original tree/loop policy while
        // filtering the target model down to diagrams containing a non-SM
        // particle.  This is what prevents TOTAL = SM + BSM from doubling SM.
        GeneralModelModifier mod(
            "C2", "ZPrime", "ZPrime", zprime_hdr.string(), std::nullopt,
            false, true
        );

        std::string model_line = "SM_Model sm;";
        mod.modifyLine(model_line);
        assert(model_line.find("ZPrime_Model sm;") != std::string::npos);

        std::string order_line = "auto order = mty::Order::TreeLevel;";
        mod.modifyLine(order_line);
        assert(order_line.find("mty::Order::TreeLevel") != std::string::npos);

        fs::path out = root / "generic_bsm.cpp";
        {
            std::ofstream f(out);
            mod.addLine(f, "#include <iostream>");
            mod.addLine(f, "using namespace sm_input;");
            mod.addLine(f, "    FeynOptions opts;");
        }
        const std::string generated = slurp(out);
        assert(generated.find("HYPERISO_MARTY_BSM_ONLY_FILTER") != std::string::npos);
        assert(generated.find("DiagramParticleType::External") != std::string::npos);
        assert(generated.find("hyperiso_marty_require_non_sm_diagram_particle(opts)") != std::string::npos);
    }

    {
        // Every coefficient template that previously requested OneLoop directly
        // is rewritten into a generic TreeLevel-first builder. Tree-only
        // templates (such as C2) remain untouched, and C9/CP9/CP10 keep their
        // specialised reg_prop implementation.
        GeneralModelModifier mod(
            "C7", "ZPrime", "ZPrime", zprime_hdr.string(), std::nullopt,
            false, true, false, true
        );
        fs::path out = root / "c7_generic_tree_first.cpp";
        std::ofstream f(out);
        const std::vector<std::string> source = {
            "#include <iostream>",
            "using namespace sm_input;",
            "int calculate_C7(Model &model, gauge::Type gauge) {",
            "    FeynOptions opts;",
            "    auto wil = model.computeWilsonCoefficients(mty::Order::OneLoop, process, opts);",
            "    Expr C7 = getWilsonCoefficient(wil, O7);",
            "    mty::Library wilsonLib(\"C7_SM\", \"libs\");",
            "    wilsonLib.cleanExistingSources();",
            "    wilsonLib.addFunction(\"C7\", C7);",
            "    defineLibPath(wilsonLib);",
            "    wilsonLib.print();",
            "    return 0;",
            "}",
            "int main() {",
            "    SM_Model sm;",
            "    return calculate_C7(sm, gauge::Type::Feynman);",
            "}",
        };
        for (auto line : source) {
            mod.modifyLine(line);
            mod.addLine(f, line);
        }
        f.close();

        const std::string generated = slurp(out);
        assert(generated.find("HYPERISO_MARTY_TREE_FIRST") != std::string::npos);
        assert(generated.find("computeWilsonCoefficients(hyperiso_marty_order") != std::string::npos);
        assert(generated.find("model.computeAmplitude(hyperiso_marty_order") != std::string::npos);
        assert(generated.find("hyperiso_marty_tree_probe.empty()") != std::string::npos);
        assert(generated.find("model.getWilsonCoefficients(hyperiso_marty_tree_probe") != std::string::npos);
        assert(generated.find("preserve_explicit_order") == std::string::npos);
        assert(generated.find("orderExternalFermions = true") != std::string::npos);
        assert(generated.find("mty::Order::TreeLevel") != std::string::npos);
        assert(generated.find("if (!hyperiso_marty_use_tree)") != std::string::npos);
        assert(generated.find("mty::Order::OneLoop") != std::string::npos);
        assert(generated.find("selected order=") != std::string::npos);
        const auto generic_undefine = generated.find("undefineNumericalValues();");
        const auto generic_model = generated.find("ZPrime_Model model;");
        assert(generic_undefine != std::string::npos);
        assert(generic_model != std::string::npos);
        assert(generic_undefine < generic_model
               && "SM inputs must be undefined before target-model construction");
        assert(generated.find("hyperiso_marty_require_non_sm_diagram_particle(opts)") != std::string::npos);
        assert(generated.find("wilsonLib.addDefaultParameter(\"s_12\", false)") != std::string::npos);
        assert(generated.find("wilsonLib.addDefaultParameter(\"s_13\", false)") != std::string::npos);
    }

    {
        // Regression test for suffixed scalar coefficients such as CQ1_E.  The
        // closing-parenthesis index must remain valid when `int` is rewritten
        // to the longer `Expr`; otherwise the generated signature becomes
        // `gaug, ... hyperiso_marty_ordere` and does not compile.
        GeneralModelModifier mod(
            "CQ1_E", "THDM", "THDM", thdm_hdr.string(), 2,
            false, true, false, true
        );
        fs::path out = root / "cq1e_generic_tree_first.cpp";
        std::ofstream f(out);
        const std::vector<std::string> source = {
            "#include <iostream>",
            "using namespace sm_input;",
            "int calculate_CQ1e(Model &model, gauge::Type gauge) {",
            "    model.getParticle(\"W\")->setGaugeChoice(gauge);",
            "    FeynOptions opts;",
            "    auto wil = model.computeWilsonCoefficients(mty::Order::OneLoop, process, opts);",
            "    Expr CQ1_e = getWilsonCoefficient(wil, Q1);",
            "    mty::Library wilsonLib(\"CQ1_E_SM\", \"libs\");",
            "    wilsonLib.addFunction(\"CQ1_E\", CQ1_e);",
            "    return 0;",
            "}",
            "int main() {",
            "    SM_Model sm;",
            "    return calculate_CQ1e(sm, gauge::Type::Feynman);",
            "}",
        };
        for (auto line : source) {
            mod.modifyLine(line);
            mod.addLine(f, line);
        }
        f.close();

        const std::string generated = slurp(out);
        assert(generated.find(
            "Expr calculate_CQ1e(Model &model, gauge::Type gauge, mty::Order hyperiso_marty_order, const std::vector<int>& hyperiso_marty_forced_fermion_order = {}) {"
        ) != std::string::npos);
        assert(generated.find("setGaugeChoice(gauge)") != std::string::npos);
        assert(generated.find("computeWilsonCoefficients(hyperiso_marty_order") != std::string::npos);
        assert(generated.find("model.computeAmplitude(hyperiso_marty_order") != std::string::npos);
        assert(generated.find("hyperiso_marty_tree_probe.empty()") != std::string::npos);
        assert(generated.find("gaug,") == std::string::npos);
        assert(generated.find("hyperiso_marty_ordere") == std::string::npos);
    }

    {
        // C10 is a specialised template: it already performs its own tree-first
        // selection and therefore bypasses the generic tree-first source rewrite.
        // The library parameter contract must nevertheless include s_12/s_13,
        // because CSL's generated 1->3 kinematics.cpp references both fields.
        GeneralModelModifier mod(
            "C10", "ZPrime", "ZPrime", zprime_hdr.string(), std::nullopt,
            false, true, false, false, MartyOrderPolicy::TREE_LEVEL_ONLY,
            {1, 0, 2, 3}, {}
        );
        fs::path out = root / "c10_specialised_library_params.cpp";
        std::ofstream f(out);
        const std::vector<std::string> source = {
            "#include <iostream>",
            "using namespace sm_input;",
            "int calculate_C10mu(Model &model, gauge::Type gauge) {",
            "    FeynOptions opts;",
            "    mty::Library wilsonLib(\"C10_SM\", \"libs\");",
            "    wilsonLib.cleanExistingSources();",
            "    wilsonLib.addFunction(\"C10\", C10_mu);",
            "    return 0;",
            "}",
        };
        for (auto line : source) {
            mod.modifyLine(line);
            mod.addLine(f, line);
        }
        f.close();

        const std::string generated = slurp(out);
        assert(generated.find("HYPERISO_MARTY_TREE_FIRST") == std::string::npos);
        const auto clean = generated.find("wilsonLib.cleanExistingSources()");
        const auto s12 = generated.find("wilsonLib.addDefaultParameter(\"s_12\", false)");
        const auto s13 = generated.find("wilsonLib.addDefaultParameter(\"s_13\", false)");
        const auto add = generated.find("wilsonLib.addFunction(\"C10\"");
        assert(clean != std::string::npos);
        assert(s12 != std::string::npos);
        assert(s13 != std::string::npos);
        assert(add != std::string::npos);
        assert(clean < s12 && s12 < s13 && s13 < add);
    }

    {
        // C9 keeps its specialised order/reg_prop source rewrite and generates
        // a tree-first main: a non-zero tree coefficient prevents OneLoop from
        // being evaluated, while loop-only models fall back to OneLoop.
        GeneralModelModifier mod(
            "C9", "ZPrime", "ZPrime", zprime_hdr.string(), std::nullopt,
            false, true
        );
        // A standalone order literal outside computeWilsonCoefficients must be
        // left untouched.  Projection-recipe helpers contain TreeLevel calls of
        // their own and must never inherit the calculate-function local order.
        std::string helper_order_line = "        mty::Order::TreeLevel,";
        mod.modifyLine(helper_order_line);
        assert(helper_order_line.find("mty::Order::TreeLevel,") != std::string::npos);
        assert(helper_order_line.find("hyperiso_marty_order") == std::string::npos);

        // The first argument of an actual multi-line Wilson call is still
        // rewritten exactly as before.
        std::string wilson_call_line = "    auto wil = model.computeWilsonCoefficients(";
        mod.modifyLine(wilson_call_line);
        std::string order_line = "        mty::Order::OneLoop,";
        mod.modifyLine(order_line);
        assert(order_line.find("hyperiso_marty_order,") != std::string::npos);

        // computeAmplitude also takes an mty::Order, but it is not the Wilson
        // call whose perturbative order is controlled by this rewrite.
        std::string amplitude_call_line = "    auto amp = model.computeAmplitude(";
        mod.modifyLine(amplitude_call_line);
        std::string amplitude_order_line = "        mty::Order::TreeLevel,";
        mod.modifyLine(amplitude_order_line);
        assert(amplitude_order_line.find("mty::Order::TreeLevel,") != std::string::npos);
        assert(amplitude_order_line.find("hyperiso_marty_order") == std::string::npos);

        std::string setter_line =
            "    hyperiso_marty_tree_level_matching = (order == mty::Order::TreeLevel);";
        mod.modifyLine(setter_line);
        assert(setter_line.find("order == mty::Order::TreeLevel") != std::string::npos);
        assert(setter_line.find("order == hyperiso_marty_order") == std::string::npos);

        fs::path out = root / "c9_tree_first.cpp";
        {
            std::ofstream f(out);
            const std::vector<std::string> source = {
                "#include <iostream>",
                "using namespace sm_input;",
                "int calculate_C9mu(Model &model, gauge::Type gauge) {",
                "    FeynOptions opts;",
                "    auto wil = model.computeWilsonCoefficients(",
                "        mty::Order::OneLoop,",
                "        process,",
                "        opts",
                "    );",
                "    Expr C9_mu = getWilsonCoefficient(wil, O9);",
                "    mty::Library wilsonLib(\"C9_SM\", \"libs\");",
                "    wilsonLib.addFunction(\"C9\", C9_mu);",
                "    return 0;",
                "}",
                "int main() {",
            };
            for (auto line : source) {
                mod.modifyLine(line);
                mod.addLine(f, line);
            }
        }
        const std::string generated = slurp(out);
        assert(generated.find("mty::Order::TreeLevel") != std::string::npos);
        assert(generated.find("if (!hyperiso_marty_use_tree_level)") != std::string::npos);
        assert(generated.find("mty::Order::OneLoop") != std::string::npos);
        assert(generated.find("selected order=") != std::string::npos);
        const auto split_undefine = generated.find("undefineNumericalValues();");
        const auto split_tree_model = generated.find("ZPrime_Model tree_model;");
        assert(split_undefine != std::string::npos);
        assert(split_tree_model != std::string::npos);
        assert(split_undefine < split_tree_model
               && "SM inputs must be undefined before split target-model construction");
        assert(generated.find("ZPrime_Model loop_model;") != std::string::npos);
        assert(generated.find(
            "undefineNumericalValues();", split_tree_model + 1
        ) != std::string::npos
               && "The isolated C9 fallback model must also start from symbolic SM inputs");
        assert(generated.find("hyperiso_marty_build_C9(tree_model") != std::string::npos);
        assert(generated.find("hyperiso_marty_build_C9(loop_model") != std::string::npos);
        assert(generated.find("model.computeAmplitude(") != std::string::npos);
        assert(generated.find("hyperiso_marty_tree_probe.empty()") != std::string::npos);
        assert(generated.find("return std::make_pair(CSL_0, std::size_t{0})") != std::string::npos);
        assert(generated.find("model.getWilsonCoefficients(hyperiso_marty_tree_probe") != std::string::npos);
        assert(generated.find("preserve_explicit_order") == std::string::npos);
        assert(generated.find("orderExternalFermions = true") != std::string::npos);
        assert(generated.find("hyperiso_marty_tree_fermion_orders()") == std::string::npos);
        assert(generated.find("hyperiso_marty_configured_fermion_order(mty::Order::TreeLevel)") != std::string::npos);
        assert(generated.find("wilsonLib.addDefaultParameter(\"s_12\", false)") != std::string::npos);
        assert(generated.find("wilsonLib.addDefaultParameter(\"s_13\", false)") != std::string::npos);
    }

    {
        // A coefficient may use distinct explicit fermion orders at TreeLevel
        // and OneLoop.  The generator must keep both values and never scan
        // permutations to select the first non-zero projection.
        GeneralModelModifier mod(
            "C9", "ZPrime", "ZPrime", zprime_hdr.string(), std::nullopt,
            false, true, false, false, MartyOrderPolicy::TREE_LEVEL_ONLY,
            {1, 0, 2, 3},
            {1, 0, 3, 2},
            {0, 3, 1, 2},
            {1, 0, 2, 3}
        );
        fs::path out = root / "c9_tree_only_order_scan.cpp";
        std::ofstream f(out);
        const std::vector<std::string> source = {
            "#include <iostream>",
            "using namespace sm_input;",
            "int calculate_C9(Model &model, gauge::Type gauge) {",
            "    FeynOptions opts;",
            "    opts.setFermionOrder({0, 2, 1, 3});",
            "    auto wil = model.computeWilsonCoefficients(",
            "        mty::Order::OneLoop,",
            "        process,",
            "        opts",
            "    );",
            "    Expr C9_mu = getWilsonCoefficient(wil, dimension6Operator(model, wil, DiracCoupling::VL, DiracCoupling::V));",
            "    mty::Library wilsonLib(\"C9_SM\", \"libs\");",
            "    wilsonLib.addFunction(\"C9\", C9_mu);",
            "    return 0;",
            "}",
            "int main() {",
        };
        for (auto line : source) {
            mod.modifyLine(line);
            mod.addLine(f, line);
        }
        f.close();
        const std::string generated = slurp(out);
        assert(generated.find("HYPERISO_MARTY_ORDER_POLICY_TREE_LEVEL_ONLY") != std::string::npos);
        assert(generated.find("static const std::vector<int> tree_order = {1, 0, 2, 3}") != std::string::npos);
        assert(generated.find("static const std::vector<int> one_loop_order = {1, 0, 3, 2}") != std::string::npos);
        assert(generated.find("static const std::vector<int> tree_order = {0, 3, 1, 2}") != std::string::npos);
        assert(generated.find("static const std::vector<int> one_loop_order = {1, 0, 2, 3}") != std::string::npos);
        assert(generated.find("hyperiso_marty_dimension6_operator(hyperiso_marty_order, model, wil") != std::string::npos);
        assert(generated.find(
            "hyperiso_marty_configured_fermion_order(mty::Order::OneLoop)"
        ) != std::string::npos);
        assert(generated.find("hyperiso_marty_tree_fermion_orders()") == std::string::npos);
        assert(generated.find("hyperiso_marty_candidate.first != CSL_0") == std::string::npos);
        assert(generated.find("hyperiso_marty_configured_fermion_order(mty::Order::TreeLevel)") != std::string::npos);
        assert(generated.find("fermion-order=") != std::string::npos);
        assert(generated.find("hyperiso_marty_order != mty::Order::TreeLevel") != std::string::npos);
    }

    {
        // AUTO must also wrap a TreeLevel-only BSM template.  The complete LO
        // coefficient is evaluated first and a distinct OneLoop order is used
        // only when the LO expression vanishes.
        GeneralModelModifier mod(
            "C5", "ZPrime", "ZPrime", zprime_hdr.string(), std::nullopt,
            false, true, false, true, MartyOrderPolicy::AUTO,
            {1, 0, 2, 3},
            {1, 0, 3, 2}
        );
        fs::path out = root / "c5_tree_template_auto.cpp";
        std::ofstream f(out);
        const std::vector<std::string> source = {
            "#include <iostream>",
            "using namespace sm_input;",
            "int calculate_C5(Model &model, gauge::Type gauge) {",
            "    FeynOptions opts;",
            "    opts.setFermionOrder({1, 0, 2, 3});",
            "    auto wil = model.computeWilsonCoefficients(mty::Order::TreeLevel, process, opts);",
            "    Expr C5 = getWilsonCoefficient(wil, O5);",
            "    mty::Library wilsonLib(\"C5_SM\", \"libs\");",
            "    wilsonLib.addFunction(\"C5\", C5);",
            "    return 0;",
            "}",
            "int main() {",
        };
        for (auto line : source) {
            mod.modifyLine(line);
            mod.addLine(f, line);
        }
        f.close();
        const std::string generated = slurp(out);
        assert(generated.find("hyperiso_marty_use_tree") != std::string::npos);
        assert(generated.find("mty::Order::TreeLevel") != std::string::npos);
        assert(generated.find("mty::Order::OneLoop") != std::string::npos);
        assert(generated.find("static const std::vector<int> tree_order = {1, 0, 2, 3}") != std::string::npos);
        assert(generated.find("static const std::vector<int> one_loop_order = {1, 0, 3, 2}") != std::string::npos);
        assert(generated.find("hyperiso_marty_require_non_sm_diagram_particle(opts)") != std::string::npos);
    }

    {
        // TREE_LEVEL_ONLY must keep a genuine structural zero instead of
        // silently replacing it by a one-loop result.
        GeneralModelModifier mod(
            "C7", "ZPrime", "ZPrime", zprime_hdr.string(), std::nullopt,
            false, true, false, true, MartyOrderPolicy::TREE_LEVEL_ONLY
        );
        fs::path out = root / "c7_tree_only.cpp";
        std::ofstream f(out);
        const std::vector<std::string> source = {
            "#include <iostream>",
            "using namespace sm_input;",
            "int calculate_C7(Model &model, gauge::Type gauge) {",
            "    FeynOptions opts;",
            "    auto wil = model.computeWilsonCoefficients(mty::Order::OneLoop, process, opts);",
            "    Expr C7 = getWilsonCoefficient(wil, O7);",
            "    mty::Library wilsonLib(\"C7_SM\", \"libs\");",
            "    wilsonLib.addFunction(\"C7\", C7);",
            "    return 0;",
            "}",
            "int main() {",
        };
        for (auto line : source) {
            mod.modifyLine(line);
            mod.addLine(f, line);
        }
        f.close();
        const std::string generated = slurp(out);
        assert(generated.find("TreeLevelOnly") != std::string::npos);
        assert(generated.find("if (!hyperiso_marty_use_tree)") == std::string::npos);
        assert(generated.find("orderExternalFermions = true") != std::string::npos);
    }

    {
        // The explicit full-target path remains available for diagnostics and
        // must keep every target diagram.  The normal WilsonBuilder path no
        // longer uses it to infer BSM through a cross-backend subtraction.
        GeneralModelModifier mod(
            "C9", "ZPrime", "ZPrime", zprime_hdr.string(), std::nullopt,
            false, true, true
        );
        fs::path out = root / "target_c9.cpp";
        {
            std::ofstream f(out);
            mod.addLine(f, "#include <iostream>");
            mod.addLine(f, "using namespace sm_input;");
            mod.addLine(f, "int calculate() {");
            mod.addLine(f, "    FeynOptions opts;");
        }
        const std::string generated = slurp(out);
        assert(generated.find("HYPERISO_MARTY_TARGET_SPLIT") != std::string::npos);
        assert(generated.find("hyperiso_marty_require_non_sm_diagram_particle(opts)") == std::string::npos);
    }

    {
        const std::string before = GeneralModelModifier::modelSignature(
            "ZPrime", zprime_hdr.string()
        );
        write_file(zprime_hdr,
            "class ZPrime_Model : public mty::Model { int cache_revision; };\n"
        );
        const std::string after = GeneralModelModifier::modelSignature(
            "ZPrime", zprime_hdr.string()
        );
        assert(before != after && "MARTY cache signatures must track model header content");
    }

    MartyRuntimeConfig::clear_external_install_path();
    std::cout << "UNIT OK\n";
    return 0;
}
