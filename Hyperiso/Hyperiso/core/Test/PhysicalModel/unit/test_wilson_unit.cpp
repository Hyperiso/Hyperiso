#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "Wilson.h"
#include "Parameter.h"
#include "IParameterProxy.h"

class StubWC : public WilsonCoefficient {
public:
    explicit StubWC(const std::string& name, const std::string& block)
        : WilsonCoefficient(name, block) {}
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<StubWC>(*this);
    }
};

class SpyProxy : public IParameterProxy<std::string, LhaID> {
public:
    mutable std::string last_block;
    mutable LhaID last_id{0};
    double ret = 42.0;

    scalar_t operator()(const std::string& blk, const LhaID& id) const override {
        last_block = blk;
        last_id    = id;
        return ret;
    }
    bool exist(const std::string&, const LhaID&) const override { return true; }
    double get_scale(const std::string&) const override { return 0.0; }
};

int main() {
    std::cout << "== Wilson BASE UNIT ==\n";

    {
        StubWC w1("C7", "B_MATCH");
        assert(w1.get_base_name() == "C7");
        assert(w1.get_type() == ContributionType::SM);

        StubWC w2("C7_THDM", "B_MATCH");
        assert(w2.get_base_name() == "C7");
        assert(w2.get_type() == ContributionType::BSM);

        StubWC w3("C10_SUSY", "B_MATCH");
        assert(w3.get_base_name() == "C10");
        assert(w3.get_type() == ContributionType::BSM);
    }

    {
        StubWC w("C9", "B_MATCH");
        auto id_lo   = w.id(QCDOrder::LO,   w.get_type());
        auto id_nlo  = w.id(QCDOrder::NLO,  w.get_type());
        auto id_nnlo = w.id(QCDOrder::NNLO, w.get_type());

        assert(w.get_lhaid(QCDOrder::LO)   == id_lo);
        assert(w.get_lhaid(QCDOrder::NLO)  == id_nlo);
        assert(w.get_lhaid(QCDOrder::NNLO) == id_nnlo);

        auto f = w.get_func(QCDOrder::LO);
        std::unordered_map<ParamId, std::shared_ptr<Parameter>> none;
        assert(f(ParamSrc(none)) == 0.0);
    }

    {
        StubWC a("C7", "B_MATCH");
        StubWC b("C7", "B_MATCH");
        assert(a == b);

        a.set_owned(true);
        assert(a != b);
        b.set_owned(true);
        assert(a == b);
    }

    {
        StubWC w("C8", "B_MATCH");
        SpyProxy spy;
        spy.ret = 123.0;

        const auto val = w.get_matching_value("LO", ContributionType::SM,
                                              std::make_shared<SpyProxy>(spy));

                                              assert(std::abs(val.real() - 123.0) < 1e-12);
        assert(std::abs(val.imag()) < 1e-12);

        auto expected_id = w.id(QCDOrder::LO, ContributionType::SM);

        auto sp = std::make_shared<SpyProxy>();
        (void)w.get_matching_value("NLO", ContributionType::SM, sp);
        assert(sp->last_block == w.get_storage_block());
        assert(sp->last_id == w.id(QCDOrder::NLO, ContributionType::SM));

        w.set_storage_block("NEW_BLOCK");
        (void)w.get_matching_value("NNLO", ContributionType::SM, sp);
        assert(sp->last_block == "NEW_BLOCK");
        assert(sp->last_id == w.id(QCDOrder::NNLO, ContributionType::SM));
    }

    std::cout << "UNIT OK\n";
    return 0;
}
