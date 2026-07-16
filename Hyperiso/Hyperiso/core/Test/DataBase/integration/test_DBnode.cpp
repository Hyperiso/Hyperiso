// test_dbnode_integration.cpp
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>
#include <cmath>

#include "DBNode.h"

static bool contains_all(const std::string& s, const std::vector<std::string>& needles) {
    for (auto& n : needles) if (s.find(n) == std::string::npos) return false;
    return true;
}

int main() {
    std::cout << "== Running INTEGRATION tests for Node ==\n";

    DBNode root;

    // Section "model"
    root.set(std::string("SM"), "model", "type");
    root.set(true, "model", "use_QED");
    root.set(91.1876, "model", "MZ");

    // Section "parameters"
    root.set(0.118, "parameters", "alpha_s");
    root.set(173.1, "parameters", "top_mass");

    // Section "metadata"
    root.set(42, "metadata", "run_id");
    root.set("integration test", "metadata", "description");

    // Observables : final list (valeurs simples stockées sous "" dans chaque élément)
    {
        std::vector<std::shared_ptr<DBNode>> obs_list;
        auto o1 = std::make_shared<DBNode>(); o1->set("BR(B→X_s γ)", "");
        auto o2 = std::make_shared<DBNode>(); o2->set("R(D*)", "");
        obs_list.push_back(o1); obs_list.push_back(o2);
        root.set(obs_list, "observables");
    }

    // Observables détaillées : liste d’objets
    {
        std::vector<std::shared_ptr<DBNode>> det;
        auto n1 = std::make_shared<DBNode>();
        n1->set("R_K", "name");
        n1->set(0.846, "value");
        n1->set(0.05, "unc");
        auto n2 = std::make_shared<DBNode>();
        n2->set("R_K*", "name");
        n2->set(0.90, "value");
        n2->set(0.06, "unc");
        det.push_back(n1); det.push_back(n2);
        root.set(det, "observables_detailed");
    }

    // Validation de quelques valeurs
    assert(std::get<BlockName>(root.get("model", "type")) == "SM");
    assert(std::abs(std::get<double>(root.get("model", "MZ")) - 91.1876) < 1e-9);
    assert(std::get<bool>(root.get("model", "use_QED")) == true);
    assert(std::get<int>(root.get("metadata", "run_id")) == 42);

    // JSON vers std::cout (visuel) et vers stream (testable)
    {
        std::stringstream ss;
        root.printJSONToStream(ss);
        const std::string json = ss.str();

        // Présence des champs clés (sans forcer le format exact des doubles)
        assert(contains_all(json, {
            "\"model\"", "\"type\"", "\"SM\"",
            "\"parameters\"", "\"alpha_s\"", "\"top_mass\"",
            "\"metadata\"", "\"run_id\"", "42", "\"description\"", "\"integration test\"",
            "\"observables\"", "BR(B", "R(D*)",
            "\"observables_detailed\"", "\"name\"", "\"value\"", "\"unc\"", "0.846"
        }));

        // Tolère 0.9 ou 0.90 selon le formatage des doubles
        bool has_RKstar_val = (json.find("0.90") != std::string::npos) || (json.find("0.9") != std::string::npos);
        assert(has_RKstar_val);
    }

    // YAML (capture stdout)
    {
        std::stringstream capture;
        auto* oldbuf = std::cout.rdbuf(capture.rdbuf());
        root.printYAML();
        std::cout.rdbuf(oldbuf);

        const std::string yaml = capture.str();

        // Vérifie les clés/valeurs essentielles (sans imposer "- ")
        std::vector<std::string> must = {
            "model:", "type:", "SM",
            "parameters:", "alpha_s:", "top_mass:",
            "metadata:", "run_id:", "description:",
            "observables:",                 // la liste existe
            "BR(B", "R(D*)",               // éléments de la "final list"
            "observables_detailed:",       // liste d’objets
            "name:", "value:", "unc:"      // clés dans les objets
        };

        // Debug sympa : imprimer le YAML et le token manquant si ça plante
        for (const auto& token : must) {
            if (yaml.find(token) == std::string::npos) {
                std::cerr << "\n--- YAML dump ---\n" << yaml << "\n";
                std::cerr << "Missing token: [" << token << "]\n";
                assert(false && "YAML is missing an expected token");
            }
        }
    }
    std::cout << "\n All DBNode integration tests passed!\n";
    return 0;
}
