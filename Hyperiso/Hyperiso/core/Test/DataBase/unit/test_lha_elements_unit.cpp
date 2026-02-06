#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <sstream>
#include <typeinfo>

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
    std::cout << "[unit] début des tests LHA elements...\n";

    {
        Prototype P = SMINPUTS; 
        std::vector<std::string> line = {"1", "0.118"};
        LhaElement<double> e(P, line);

        assert(dapprox(e.getValue(), 0.118));
        assert(dapprox(e.getScale(), 0.0));
        assert(e.getScheme() == RenormalizationScheme::NONE);

        std::string s = e.toString(); 
        auto toks = split_tabs(s);

        assert(toks.size() == 2);
        assert(dapprox(std::stod(toks[1]), 0.118));

        auto node = e.toDBNode();
        assert(node != nullptr);
    }

    {
        Prototype P = FMASS; 
        std::vector<std::string> line = {"5", "4.2", "91.1876", "1"}; 
        LhaElement<double> e(P, line);

        assert(dapprox(e.getValue(), 4.2));
        assert(dapprox(e.getScale(), 91.1876));

        assert(e.getScheme() == static_cast<RenormalizationScheme>(1));

        auto toks = split_tabs(e.toString());

        assert(toks.size() == 4);
        assert(dapprox(std::stod(toks[1]), 4.2));
        assert(dapprox(std::stod(toks[2]), 91.1876));

        assert(std::stoi(toks[3]) == 1);
    }

    {
        Prototype P = HMIX;
        std::vector<std::string> line = {"1000.0", "1", "0.5"};
        LhaElement<double> e(P, line);

        assert(dapprox(e.getValue(), 0.5));
        assert(dapprox(e.getScale(), 1000.0));
        assert(e.getScheme() == RenormalizationScheme::NONE);

        auto toks = split_tabs(e.toString());

        assert(toks.size() == 3);
        assert(dapprox(std::stod(toks[1]), 0.5));
        assert(dapprox(std::stod(toks[2]), 1000.0));
    }

    {
        Prototype P = FCINFO;
        std::vector<std::string> line = {"1", "SoftSUSY"};
        LhaElement<std::string> e(P, line);

        assert(e.getValue() == "SoftSUSY");
        assert(dapprox(e.getScale(), 0.0));
        assert(e.getScheme() == RenormalizationScheme::NONE);

        auto toks = split_tabs(e.toString());

        assert(toks.size() == 2);
        assert(toks[1] == "SoftSUSY");
    }


    {
        Prototype P = SMINPUTS;
        std::vector<std::string> bad = {"1", "abc"}; 
        bool thrown = false;
        try {
            LhaElement<double> e(P, bad);
            (void)e;
        } catch (const std::runtime_error&) {
            thrown = true;
        }
        assert(thrown);
    }

    {
        Prototype P = FMASS;
        std::vector<std::string> bad = {"5", "4.2", "QX", "1"};
        bool thrown = false;
        try {
            LhaElement<double> e(P, bad);
            (void)e;
        } catch (const std::invalid_argument&) {
            thrown = true;
        } catch (...) {}
        assert(thrown);
    }

    {

        {
            std::vector<std::string> line = {"1", "SPHENO"};
            auto ptr = LhaElementFactory::createElement(FCINFO, line);
            assert(ptr != nullptr);

            auto s = std::dynamic_pointer_cast< LhaElement<std::string> >(ptr);
            assert(s != nullptr);
            assert(s->getValue() == "SPHENO");
        }

        {
            std::vector<std::string> line = {"125.0", "1", "0.12"};
            auto ptr = LhaElementFactory::createElement(HMIX, line);
            assert(ptr != nullptr);
            auto d = std::dynamic_pointer_cast< LhaElement<double> >(ptr);
            assert(d != nullptr);
            assert(dapprox(d->getValue(), 0.12));
            assert(dapprox(d->getScale(), 125.0));
        }

        {
            Prototype SPINFO_ = Prototype{"SPINFO"};
            std::vector<std::string> line = {"1", "Version1.0"};
            auto ptr = LhaElementFactory::createElement(SPINFO_, line);
            assert(ptr != nullptr);
            auto s = std::dynamic_pointer_cast< LhaElement<std::string> >(ptr);
            assert(s != nullptr);
            assert(s->getValue() == "Version1.0");
        }
    }

    std::cout << "All unit tests passed.\n";
    return 0;
}
