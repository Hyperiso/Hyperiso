#include "HyperisoMaster.h"
#include "ParameterProvider.h"
#include "Include.h"
#include "Logger.h"
#include "CompositeParamCreator.h"
#include "QCDProvider.h"
#include "ParameterSetter.h"
#include "BlockProvider.h"
#include "ParamBlockWriter.h"
#include "JsonParser.h"
#include "YamlParser.h"
#include "LhaParser.h"
#include "FileWriter.h"

static ParamId make_id(ParameterType t, const char* block, int code) {
    return ParamId(t, BlockName(block), LhaID(code));
}

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoConfig config;
    config.model = Model::SM;

    config.flags[ExternalFlag::HAS_WILSON_INPUT] = true;

    HyperisoMaster hi;
    hi.init("lha/testInput.flha", config);

    BlockProvider().log_all_blocks(ParameterType::WILSON);
    
    ParameterProvider sm {ParameterType::SM};
    ParameterProvider wil {ParameterType::WILSON};
    ParameterProvider obs {ParameterType::OBSERVABLE};

    auto obs_p = obs.get_parameter({ParameterType::OBSERVABLE, "FOBS", {521, 2, 3, 321, 13}});
    LOG_INFO(*obs_p);
    auto bin = obs_p->get_bin();
    LOG_INFO("A_FB(B > K mu mu) = [", bin.first, ",", bin.second, "] =", obs_p->get_val());
    

    // 521_11_3_423_-15_16

    auto obs_p2 = obs.get_parameter({ParameterType::OBSERVABLE, "FOBS", {521, 11, 3, 423, -15, 16}});
    LOG_INFO(*obs_p2);
    auto bin2 = obs_p2->get_bin();
    LOG_INFO("A_FB(B > K mu mu) = [", bin2.first, ",", bin2.second, "] =", obs_p2->get_val());

    LOG_INFO(sm("SMINPUTS", 6));

    CompositeParamAdapter cpc;
    std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"SMINPUTS", "MASS"}}};

    auto func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
        QCDProvider qcd;
        auto xt = std::pow(qcd(MassConfig{6, 80, MassType::MSBAR, MassType::POLE}) / src.get_val("MASS",24), 2);
        dep_block->store_or_assign(1, std::make_shared<Parameter>(Parameter({ParameterType::WILSON, "WPARAM", 1}, xt, 0., 0.)));
    };
    cpc.add_block_dependency("WPARAM", src, ParameterType::WILSON, func);

    ParameterSetter ps;
    LOG_INFO("Before: m_W =", sm("MASS", 24), ", x_t =", wil("WPARAM", 1));
    ps.mutate({ParameterType::SM, "MASS", 24}, 100);
    LOG_INFO("After: m_W =", sm("MASS", 24), ", x_t =", wil("WPARAM", 1));

    std::shared_ptr<DBNode> node = std::make_shared<DBNode>();
    std::shared_ptr<BlockAccessor> acc = MemoryManager::GetInstance()->extract_block_accessor();
    ParamBlockWriter().write(node, acc);

    std::cout << "writing to file.." << std::endl;
    hi.switch_lha("lha/si_input.flha", config);
    JSONParser().writeToFile("test.json", node);
    YAMLParser().writeToFile("test.yaml", node);
    LhaParser().writeToFile("test.flha", node);

    FileWriter().write("test2.json");
    FileWriter().write("test2.yaml");
    FileWriter().write("test2.flha");


    try {
        // -------------------------
        // 1) Source block SRC: x
        // -------------------------
        auto src = std::make_shared<Block>();
        src->blockname = "SRC";

        ParamId x_id = make_id(ParameterType::SM, "SRC", 1);
        auto x = std::make_shared<Parameter>(x_id, 1.0, 0.0, 0.0);
        src->store(LhaID(1), x);

        // -------------------------
        // 2) DependentBlock DBLK: a = x + 10
        //    DepUpdateFunc = void(const BlockSrc&, shared_ptr<DependentBlock>)
        // -------------------------
        std::unordered_map<std::string, std::shared_ptr<Block>> blk_sources;
        blk_sources.emplace("SRC", src);

        DepUpdateFunc blk_recalc = DepUpdateFunc(
        [](const BlockSrc& srcs, std::shared_ptr<DependentBlock> self) {
            auto s = srcs.block("SRC");
            double xv = s->retrieve(LhaID(1))->get_val();

            if (!self->contains(LhaID(1))) {
                ParamId a_id(ParameterType::SM, BlockName("DBLK"), LhaID(1));
                auto a = std::make_shared<Parameter>(a_id, xv + 10.0, 0.0, 0.0);
                self->store(LhaID(1), a);
            } else {
                // IMPORTANT: ne pas store() ici (sinon ça n’écrase jamais)
                self->retrieve(LhaID(1))->set_expected(xv + 10.0);
            }
        }
    );

        auto dblk = std::make_shared<DependentBlock>(blk_sources, blk_recalc);
        dblk->blockname = "DBLK";
        dblk->init();
        dblk->update();

        // -------------------------
        // 3) DependentParameter y = 2*a + 1, dépend de DBLK:1
        //    DepParamUpdateFunc = void(const ParamSrc&, shared_ptr<DependentParameter>)
        // -------------------------
        ParamId a_id = make_id(ParameterType::SM, "DBLK", 1);
        ParamId y_id = make_id(ParameterType::SM, "Y", 1);

        std::unordered_map<ParamId, std::shared_ptr<Parameter>> psources;
        psources.emplace(a_id, dblk->retrieve(LhaID(1))); // NB: on pointe vers le param "a" courant

        DepParamUpdateFunc y_recalc = DepParamUpdateFunc(
            [a_id](const ParamSrc& srcs, std::shared_ptr<DependentParameter> self) {
                double av = srcs.get_val(a_id);       // via SourcesView
                self->set_expected(2.0 * av + 1.0);
            }
        );

        auto y = std::make_shared<DependentParameter>(y_id, psources, y_recalc);
        y->init();
        y->update();

        // petit bloc conteneur
        auto out = std::make_shared<Block>();
        out->blockname = "OUT";
        out->store(LhaID(1), y);

        // -------------------------
        // 4) Check initial: x=1 => a=11 => y=23
        // -------------------------
        double y0 = out->retrieve(LhaID(1))->get_val();
        std::cout << "y0 = " << y0 << " (expect 23)\n";

        // -------------------------
        // 5) Change x => 2  => a=12 => y=25
        // -------------------------
        src->assign(LhaID(1), 2.0);
        
        double a1 = dblk->retrieve(LhaID(1))->get_val();
        std::cout << "a after x=2 => " << a1 << " (expect 12)\n";
        double y1 = out->retrieve(LhaID(1))->get_val();
        std::cout << "y1 = " << y1 << " (expect 25)\n";

        if (std::abs(y0 - 23.0) > 1e-9) {
            std::cerr << "[FAIL] y0 != 23\n";
            return 1;
        }
        if (std::abs(y1 - 25.0) > 1e-9) {
            std::cerr << "[FAIL] y1 != 25 -> cascade broken (block->param)\n";
            return 2;
        }

        std::cout << "[OK] DependentBlock + DependentParameter cascade works.\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return 99;
    }
    return 0;
}