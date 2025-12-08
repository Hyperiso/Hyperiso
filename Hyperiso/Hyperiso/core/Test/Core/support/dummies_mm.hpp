#pragma once
#include <filesystem>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <string>
#include <chrono>

#include "IDataLoader.h"
#include "ISpectrumCalculator.h"
#include "BlockAccessor.h"
#include "CorrelationRepo.h"
#include "Paths.h"
#include "Include.h"

namespace fs = std::filesystem;

struct TestPathsProvider : public IPathsProvider {
    fs::path root;
    explicit TestPathsProvider(fs::path r) : root(std::move(r)) {}

    fs::path assets_root() const override { return root; }

    fs::path default_param_values() const override { return root/"default"/"parameters.json"; }
    fs::path default_obs_values()   const override { return root/"default"/"observables.json"; }
    fs::path default_param_corr()   const override { return root/"default"/"parameters_corr.json"; }
    fs::path default_obs_corr()     const override { return root/"default"/"observables_corr.json"; }

    fs::path user_sm_params()       const override { return root/"input_files"/"parameters"/"sm.yaml"; }
    fs::path user_flavor_params()   const override { return root/"input_files"/"parameters"/"flavor.yaml"; }
    fs::path user_decay_params()    const override { return root/"input_files"/"parameters"/"decay.yaml"; }
    fs::path user_obs_values()      const override { return root/"input_files"/"observables"/"observables.yaml"; }
    fs::path user_param_corr()      const override { return root/"input_files"/"parameters"/"correlations.yaml"; }
    fs::path user_obs_corr()        const override { return root/"input_files"/"observables"/"correlations.yaml"; }

    fs::path spectrum_dir()         const override { return root/"spectrum"; }
    fs::path template_dir_path()         const override { return root/"spectrum"; }
    fs::path param_mapping_dir_path()         const override { return root/"spectrum"; }
};

inline fs::path& sandbox_root_singleton() {
    static fs::path root; 
    return root;
}

inline void touch_file(const fs::path& p, std::string_view content = "{}\n") {
    std::error_code ec{};
    fs::create_directories(p.parent_path(), ec);
    std::ofstream f(p);
    f << content;
}

inline fs::path prepare_assets_for_mm(const fs::path& lha_rel) {
    auto& root = sandbox_root_singleton();
    if (root.empty()) {

        auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        root = fs::temp_directory_path() / ("hyperiso_assets_sandbox_" + std::to_string(stamp));
        fs::create_directories(root);

        fs::create_directories(root/"default");
        fs::create_directories(root/"input_files"/"parameters");
        fs::create_directories(root/"input_files"/"observables");
        fs::create_directories(root/"spectrum");

        touch_file(root/"default"/"parameters.json",      R"({ "MASS": {}, "GAUGE": {} })");
        touch_file(root/"default"/"observables.json",     R"({ "FOBS": {} })");
        touch_file(root/"default"/"parameters_corr.json", R"({})");
        touch_file(root/"default"/"observables_corr.json",R"({})");

        touch_file(root/"input_files"/"parameters"/"sm.yaml",        "sm: {}\n");
        touch_file(root/"input_files"/"parameters"/"flavor.yaml",    "flavor: {}\n");
        touch_file(root/"input_files"/"parameters"/"decay.yaml",     "decay: {}\n");
        touch_file(root/"input_files"/"observables"/"observables.yaml","observables: {}\n");
        touch_file(root/"input_files"/"parameters"/"correlations.yaml","{}\n");
    }

    const fs::path lha_abs = (root / lha_rel).lexically_normal();
    touch_file(lha_abs, "BLOCK MASS\n 25 1.250000e+02\n");

    return root;
}

inline double hash_to_double(const fs::path& p) {
    auto s = p.string();
    size_t h = std::hash<std::string>{}(s);
    return 1.0 + static_cast<double>((h % 9900)) / 100.0;
}

struct DummyBA : public IDataLoader<BlockAccessor> {
    void load(std::shared_ptr<BlockAccessor> dest, fs::path src_file) override {

        double base = hash_to_double(src_file);

        auto mass = std::make_shared<Block>();
        mass->blockname = "MASS";
        mass->store(1, std::make_shared<Parameter>(ParamId("MASS", 1), 0.0047, 0.0, 0.0));
        mass->store(2, std::make_shared<Parameter>(ParamId("MASS", 2), 0.00216, 0.0, 0.0));
        mass->store(3, std::make_shared<Parameter>(ParamId("MASS", 3), 0.0935, 0.0, 0.0));
        mass->store(4, std::make_shared<Parameter>(ParamId("MASS", 4), 1.273, 0.0, 0.0));
        mass->store(25, std::make_shared<Parameter>(ParamId("MASS", 25), 125.0, 0.0, 0.0));
        mass->set_scale(91.1876);

        auto gauge = std::make_shared<Block>();
        gauge->blockname = "GAUGE";
        gauge->store(1, std::make_shared<Parameter>(ParamId("GAUGE", 1), base, 0.0, 0.0));

        auto sminputs = std::make_shared<Block>();
        sminputs->blockname = "SMINPUTS";
        sminputs->store(1, std::make_shared<Parameter>(ParamId("SMINPUTS", 1), 1.27934000e+02, 0.0, 0.0));
        sminputs->store(3, std::make_shared<Parameter>(ParamId("SMINPUTS", 3), 1.17200e-01, 0.0, 0.0));
        sminputs->store(4, std::make_shared<Parameter>(ParamId("SMINPUTS", 4), 9.11876000e+01, 0, 0));
        sminputs->store(5, std::make_shared<Parameter>(ParamId("SMINPUTS", 5), 4.180000000e+00, 0, 0));
        sminputs->store(6, std::make_shared<Parameter>(ParamId("SMINPUTS", 6), 1.72500000e+02, 0, 0));
        auto qcd = std::make_shared<Block>();
        qcd->blockname = "QCD";
        qcd->store(1, std::make_shared<Parameter>(ParamId("QCD", 1), 0.118, 0.0, 0.0));

        auto vckmin = std::make_shared<Block>();
        vckmin->blockname = "VCKMIN";
        vckmin->store(1, std::make_shared<Parameter>(ParamId("VCKMIN", 1), 0.225, 0.0, 0.0)); // lambda
        vckmin->store(2, std::make_shared<Parameter>(ParamId("VCKMIN", 2), 0.82,  0.0, 0.0)); // A
        vckmin->store(3, std::make_shared<Parameter>(ParamId("VCKMIN", 3), 0.135, 0.0, 0.0)); // rho
        vckmin->store(4, std::make_shared<Parameter>(ParamId("VCKMIN", 4), 0.349, 0.0, 0.0)); // eta

        auto fobs = std::make_shared<Block>();
        fobs->blockname = "FOBS";
        fobs->store(1, std::make_shared<Parameter>(ParamId("FOBS", 1), 0.0, 0.0, 0.0));

        auto fw = std::make_shared<Block>();
        fw->blockname = "FWCOEF";
        fw->set_scale(80.0);

        auto few = std::make_shared<Block>();
        few->blockname = "EW_SCALE";
        few->store(1, std::make_shared<Parameter>(ParamId("EW_SCALE", 1), 81., 0.0, 0.0));
        few->set_scale(80.0);

        auto fb = std::make_shared<Block>();
        fb->blockname = "B_SCALE";
        fb->store(1, std::make_shared<Parameter>(ParamId("EW_SCALE", 1), 4.18, 0.0, 0.0));
        fb->set_scale(80.0);

        auto fd = std::make_shared<Block>();
        fd->blockname = "D_SCALE";
        fd->store(1, std::make_shared<Parameter>(ParamId("EW_SCALE", 1), 0.5, 0.0, 0.0));
        fd->set_scale(80.0);

        dest->emplace("MASS",  mass);
        dest->emplace("GAUGE", gauge);
        dest->emplace("SMINPUTS", sminputs);
        dest->emplace("QCD",      qcd);
        dest->emplace("VCKMIN",   vckmin);
        dest->emplace("FOBS",  fobs);
        dest->emplace("FWCOEF",fw);
        dest->emplace("EW_SCALE",few);
        dest->emplace("B_SCALE",fb);
        dest->emplace("D_SCALE",fd);
    }
};

template<typename T>
struct DummyCorr : public IDataLoader<CorrelationMatrixPair<T>> {
    void load(std::shared_ptr<CorrelationMatrixPair<T>> dest, fs::path /*src_file*/) override {
        (void)dest;
    }
};

struct DummySpectrum : public ISpectrumCalculator {
    void calculate_spectrum(fs::path in, fs::path out, Model /*m*/) override {
        std::error_code ec{};
        fs::create_directories(out.parent_path(), ec);

        std::ifstream ifs(in, std::ios::binary);
        std::ofstream ofs(out, std::ios::binary);
        ofs << ifs.rdbuf();
    }
};
