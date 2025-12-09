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
    void calculate(std::string, std::string, double, std::string, bool = false) override {}
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

// Petit helper pour construire un BuildContext complet
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
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr;

    WilsonGroupAdapterConfig adapters(
        wilson_proxy,
        block_composer,
        use_marty_api,
        marty_model_name,
        marty_model_path,
        marty_proxy
    );

    // ⚠️ Construction agrégée, pas de BuildContext() par défaut
    return BuildContext{adapters, model, backend, contrib, gid};
}

} // namespace