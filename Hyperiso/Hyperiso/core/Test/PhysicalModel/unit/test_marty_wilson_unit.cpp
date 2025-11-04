// testMartyWilsonUnit.cpp
#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "MartyWilson.h"                 // classe testée
#include "InterpretedParam.h"            // pour le type utilisé par le proxy
#include "IParameterProxy.h"               // IParameterProxy
#include "IMartyWilsonProxy.h"           // IMartyWilsonProxy<InterpretedParam>
#include "Parameter.h"                   // ParamId, Parameter
#include "Include.h"                     // LhaID, enums
#include "CSVReader.h"                   // (pour signature), pas essentiel ici

namespace fs = std::filesystem;

// --- util ---
static void ensure_parent(const fs::path& p) {
    fs::create_directories(p.parent_path());
}
static std::string slurp(const fs::path& p){
    std::ifstream f(p); std::stringstream ss; ss << f.rdbuf(); return ss.str();
}

// Proxy Marty factice : écrit un CSV (synchrone) à l’emplacement fourni
class FakeMartyProxy : public IMartyWilsonProxy<InterpretedParam> {
public:
    fs::path csv_path;
    // valeurs à écrire (par coef)
    std::unordered_map<std::string, std::pair<double,double>> values{
        {"C7",  { -1.942999e-01, 3.438710e-07 }},
        {"C10", { -1.936392e+00, -3.200705e-02 }},
    };

    // capture des derniers arguments
    mutable std::string last_wilson, last_model;
    mutable double last_Q = 0.0;
    mutable fs::path last_model_path;

    explicit FakeMartyProxy(const fs::path& csv) : csv_path(csv) {}

    void calculate(std::string wilson, std::string model, double Q_match,
                   std::string model_path, bool /*new_params*/ = false) override
    {
        last_wilson = wilson; last_model = model; last_Q = Q_match; last_model_path = model_path;

        // Écrit un CSV avec une seule ligne au Q demandé
        ensure_parent(csv_path);
        std::ofstream out(csv_path);
        // headers: Q_match, <wilson>_real, <wilson>_img
        out << "Q_match," << wilson << "_real," << wilson << "_img\n";
        auto it = values.find(wilson);
        double re = (it!=values.end() ? it->second.first : 1.23);
        double im = (it!=values.end() ? it->second.second: -4.56);

        // volontairement avec des espaces / notation scientifique
        out << Q_match << ",  " << std::scientific << re << " , " << std::scientific << im << "\n";
        out.flush();
    }

    std::set<std::string> get_special_blocks() override {
        // rien de spécial à filtrer → toutes les dépendances remontent
        return {};
    }

    std::unordered_set<InterpretedParam> get_dependencies(std::string wilson) override {
        // Renvoie 2 dépendances (SM + BSM), pour tester le remplissage des sources
        // (wilson ignoré ici)
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

    // On force un chemin CSV test
    const fs::path tmp = fs::temp_directory_path() / "mw_unit" / "SM_wilson.csv";

    // Proxy Marty factice (écrit dans tmp)
    auto mw_proxy = std::make_shared<FakeMartyProxy>(tmp);

    // Prépare la config
    // LhaID pour C7 LO: (305,4422,0,0) — correspond classiquement à C7
    MartyWilsonConfig cfg(
        /*model_name*/ "SM",
        /*coeff_id*/   LhaID(305,4422,0,0),
        /*storage*/    GroupMapper::str(WGroup::B, ScaleType::MATCHING),
        /*model_path*/ "/dev/null",
        /*proxy*/      mw_proxy
    );
    // IMPORTANT: on redirige le csv_path vers /tmp
    cfg.csv_path = tmp.string();

    // Instancie le coefficient
    MartyWilson mw(cfg);

    // Après construction, le compute(LO) a été appelé avec EW_SCALE=1 (dummy)
    // Mais on va le rappeler avec un Q explicite et vérifier le retour
    auto fLO = mw.get_func(QCDOrder::LO);

    // Prépare les sources: uniquement EW_SCALE → lambda ajoute les deps du proxy
    ParamId pid{ParameterType::WILSON, "EW_SCALE", 1};
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> src{
        { pid, P(81.0) }     // Q_match = 81
    };

    scalar_t val = fLO(ParamSrc(src));
    // On s’attend à trouver les valeurs que FakeMartyProxy a écrites pour C7, Q=81
    // (-1.942999e-01, 3.438710e-07)
    assert(std::abs(val.real() - (-1.942999e-01)) < 1e-12);
    assert(std::abs(val.imag() - ( 3.438710e-07)) < 1e-12);

    // Les sources doivent inclure EW_SCALE + nos 2 deps
    auto S = mw.get_sources(QCDOrder::LO);
    assert(S.count(ParamId{ParameterType::WILSON, "EW_SCALE", LhaID(1)}) == 1);
    assert(S.count(ParamId{ParameterType::SM,     "SMINPUTS", LhaID(7,1)}) == 1);
    assert(S.count(ParamId{ParameterType::BSM,    "XBLK",     LhaID(1)})   == 1);

    // Vérifie que calculate a bien été appelé avec (wilson="C7", model="SM", Q=81)
    assert(mw_proxy->last_wilson == "C7");
    assert(mw_proxy->last_model  == "SM");
    assert(std::abs(mw_proxy->last_Q - 81.0) < 1e-12);

    std::cout << "✅ UNIT OK\n";
    return 0;
}
