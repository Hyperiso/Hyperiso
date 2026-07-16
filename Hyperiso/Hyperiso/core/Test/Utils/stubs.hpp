#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <utility>

// Stubs pour les tests
#include "IBlockComposer.h"
#include "ICoreAPI.h"
#include "IMartyWilsonProxy.h"
#include "IParameterProxy.h"
#include "IParamSetter.h"
#include "WilsonGroupAdapterConfig.h"
#include "InterpretedParam.h"
#include "GroupDefinition.h"

namespace {

struct DummyBlockComposer : IBlockComposer {
    void compose_block(const std::string&, const std::unordered_map<ParameterType, std::vector<std::string>>&,
                       const DepUpdateFunc&) override {}
    void compose_parameter(const ParamId&, const std::unordered_set<ParamId>&,
                           const DepParamUpdateFunc&) override {}
    void remove_block(const std::string&) override {}
    void update(const std::string&) override {}
    void remove_all_composed_blocks() override {}
};

template<typename T>
struct DummyCoreAPI : ICoreAPI<T> {
    T value;
    DummyCoreAPI(T v) : value(std::move(v)) {}
    T get() override { return value; }
};

template<typename T>
struct DummyMartyProxy : IMartyWilsonProxy<T> {
    static void write_result(const std::string& wilson,
                             const std::string& output_model,
                             double matching_scale) {
        const auto directory = std::filesystem::temp_directory_path() / "pyhyperiso" / "MartyTemp";
        std::filesystem::create_directories(directory);
        std::ofstream output(directory / (output_model + "_wilson.csv"));
        if (!output) {
            throw std::runtime_error("DummyMartyProxy could not create its test CSV");
        }

        output << "Q_match,"
               << wilson << "_real," << wilson << "_img,"
               << wilson << "_SM_SPLIT_real," << wilson << "_SM_SPLIT_img,"
               << wilson << "_BSM_SPLIT_real," << wilson << "_BSM_SPLIT_img,"
               << wilson << "_TOTAL_SPLIT_real," << wilson << "_TOTAL_SPLIT_img\n";
        output << std::setprecision(17) << matching_scale
               << ",0,0,0,0,0,0,0,0\n";
    }

    void calculate(std::string wilson,
                   std::string model,
                   double matching_scale,
                   std::string model_path) override {
        calculate(std::move(wilson), model, model, matching_scale, std::move(model_path), false, false);
    }

    void calculate(std::string wilson,
                   std::string output_model,
                   std::string,
                   double matching_scale,
                   std::string,
                   bool,
                   bool) override {
        write_result(wilson, output_model, matching_scale);
    }

    std::set<std::string> get_special_blocks() override { return {}; }
    std::unordered_set<T> get_dependencies(std::string) override { return {}; }
};

template<typename B, typename I>
struct DummyParamProxy : IParameterProxy<B,I> {
    scalar_t operator()(const B&, const I&) const override { return 0.0; }
    bool exist(const B&, const I&) const override { return false; }
    double get_scale(const B&) const override { return 0.0; }
};

template<typename T>
struct DummyParamSetter : IParamSetter<T> {
    DummyParamSetter(T init) { this->param = init; }
    void set(double) override {}
    void switch_param(T) override {}
};

inline BuildContext make_ctx(Model model,
                             Backend backend,
                             ContributionType contrib,
                             WGroupId gid)
{
    using std::make_shared;

    auto wilson_proxy      = make_shared<DummyParamProxy<std::string, LhaID>>();
    auto block_composer    = make_shared<DummyBlockComposer>();
    auto use_marty_api     = make_shared<DummyCoreAPI<bool>>(false);
    auto marty_model_name  = make_shared<DummyCoreAPI<std::string>>("dummy_model");
    auto marty_model_path  = make_shared<DummyCoreAPI<fs::path>>(fs::path{"/tmp"});
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy =
        make_shared<DummyMartyProxy<InterpretedParam>>();

    WilsonGroupAdapterConfig adapters(
        wilson_proxy,
        block_composer,
        use_marty_api,
        marty_model_name,
        marty_model_path,
        marty_proxy
    );

    return BuildContext{adapters, model, backend, contrib, gid};
}

} // namespace