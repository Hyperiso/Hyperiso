// test_lha_elements_integration.cpp
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>

#include "LhaBlockPrototype.h"
#include "lha_elements.h"

static std::vector<std::string> split_tabs(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s) {
        if (ch == '\t' || ch == '\n') {
            if (!cur.empty()) out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}
static bool dapprox(double a, double b, double eps=1e-12) {
    return std::abs(a - b) <= eps * (1.0 + std::abs(a) + std::abs(b));
}

int main() {
    std::cout << "[integration] début des tests LHA elements...\n";

    // 1) Création hétérogène via factory
    {
        std::vector<std::shared_ptr<AbstractElement>> elems;

        // FCINFO (string)
        elems.push_back(LhaElementFactory::createElement(FCINFO, {"1", "GeneratorX"}));

        // FMASS (double + scale + scheme)
        elems.push_back(LhaElementFactory::createElement(FMASS, {"5", "4.18", "91.1876", "1"}));

        // HMIX (double + global scale)
        elems.push_back(LhaElementFactory::createElement(HMIX, {"1000.0", "1", "0.10"}));

        // Vérifs types
        size_t nStr = 0, nDbl = 0;
        for (auto& e : elems) {
            if (std::dynamic_pointer_cast<LhaElement<std::string>>(e)) ++nStr;
            if (std::dynamic_pointer_cast<LhaElement<double>>(e)) ++nDbl;
        }
        assert(nStr == 1 && nDbl == 2);

        // Vérifs de quelques champs via toString()
        {
            // FMASS
            auto f = std::dynamic_pointer_cast<LhaElement<double>>(elems[1]);
            auto toks = split_tabs(f->toString());
            // ID + val + scale + scheme
            assert(toks.size() == 4);
            assert(dapprox(std::stod(toks[1]), 4.18));
            assert(dapprox(std::stod(toks[2]), 91.1876));
            assert(std::stoi(toks[3]) == 1);
        }
        {
            // HMIX
            auto h = std::dynamic_pointer_cast<LhaElement<double>>(elems[2]);
            auto toks = split_tabs(h->toString());
            // ID + val + scale
            assert(toks.size() == 3);
            assert(dapprox(std::stod(toks[2]), 1000.0));
        }
    }

    // 2) Cohérence global scale : plusieurs lignes HMIX partagent Q (col 0)
    {
        Prototype P = HMIX;
        std::vector<std::string> line1 = {"246.22", "2", "0.5"};
        std::vector<std::string> line2 = {"246.22", "3", "1.2"};
        LhaElement<double> e1(P, line1);
        LhaElement<double> e2(P, line2);

        assert(dapprox(e1.getScale(), 246.22));
        assert(dapprox(e2.getScale(), 246.22));
        assert(dapprox(e1.getValue(), 0.5));
        assert(dapprox(e2.getValue(), 1.2));
    }

    // 3) Mélange FLHA/LHA avec schémas et échelles (FMASS + FCONST)
    {
        // FMASS: itemCount=4, valueIdx=1, scaleIdx=2, rgIdx=3
        LhaElement<double> m(FMASS, {"13", "0.105658", "2.0", "2"}); // 2 -> DRBAR
        assert(dapprox(m.getValue(), 0.105658));
        assert(dapprox(m.getScale(), 2.0));
        assert(m.getScheme() == static_cast<RenormalizationScheme>(2));

        // FCONST: itemCount=5, valueIdx=2, scaleIdx=3, rgIdx=4
        LhaElement<double> c(FCONST, {"11", "1", "0.00729735", "91.1876", "0"}); // 0 -> POLE
        assert(dapprox(c.getValue(), 0.00729735));
        assert(dapprox(c.getScale(), 91.1876));
        assert(c.getScheme() == static_cast<RenormalizationScheme>(0));

        // Vérif toString() : nb de champs 4 pour les deux
        auto t1 = split_tabs(m.toString());
        auto t2 = split_tabs(c.toString());
        assert(t1.size() == 4);
        assert(t2.size() == 4);
    }

    // 4) Schéma hors-plage (stoi out_of_range) -> exception
    {
        Prototype P = FMASS;
        std::string big = "999999999999999999999999999999";
        bool thrown = false;
        try {
            LhaElement<double> e(P, {"5", "4.2", "91.1876", big});
            (void)e;
        } catch (const std::out_of_range&) {
            thrown = true;
        } catch (...) {}
        assert(thrown);
    }

    std::cout << "All integration tests passed.\n";
    return 0;
}
