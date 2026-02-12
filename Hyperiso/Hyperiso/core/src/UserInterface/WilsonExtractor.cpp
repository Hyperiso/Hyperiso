#include "WilsonExtractor.h"

std::string wc_key(const std::string& coefName, QCDOrder o, ContributionType c) {
    std::string k = "WC:";
    k += coefName;
    k += ":";
    k += OrderMapper::str(o);
    k += ":";
    k += ContributionTypeMapper::str(c);
    return k;
}


std::vector<std::string> WilsonCoeffExtractor::schema(const OutputSpec& /*spec*/) const {
    std::vector<std::string> cols;
    cols.reserve(coefs_.size() * 6);

    std::vector<WCoefId> sorted(coefs_.begin(), coefs_.end());
    std::sort(sorted.begin(), sorted.end(), [](WCoefId a, WCoefId b){
        return WCoefMapper::str(a) < WCoefMapper::str(b);
    });

    for (auto coef : sorted) {
        auto name = WCoefMapper::str(coef);
        push_cols_for_coef(cols, name);
    }
    return cols;
}

void WilsonCoeffExtractor::extract(std::unordered_map<std::string, Value>& outY,
                const OutputSpec& /*spec*/) const {
    for (auto coef : coefs_) {
        auto name = WCoefMapper::str(coef);
        add_one(outY, coef, name, QCDOrder::LO, ContributionType::SM);
        add_one(outY, coef, name, QCDOrder::LO, ContributionType::BSM);

        if (wbc_.order > QCDOrder::LO) {
            add_one(outY, coef, name, QCDOrder::NLO, ContributionType::SM);
            add_one(outY, coef, name, QCDOrder::NLO, ContributionType::BSM);

            if (wbc_.order > QCDOrder::NLO) {
                add_one(outY, coef, name, QCDOrder::NNLO, ContributionType::SM);
                add_one(outY, coef, name, QCDOrder::NNLO, ContributionType::BSM);
            }
        }
    }
}

void WilsonCoeffExtractor::add_one(std::unordered_map<std::string, Value>& outY,
                WCoefId coef,
                const std::string& name,
                QCDOrder ord,
                ContributionType ct) const
{
    const WGroup  grp = WCoefMapper::group_of(coef);
    const WCoef en = WCoefMapper::enum_of(coef).value();
    double val = wi_.getM(grp, en, ord, ct);
    outY[wc_key(name, ord, ct)] = val;
}

void WilsonCoeffExtractor::push_cols_for_coef(std::vector<std::string>& cols,
                        const std::string& name) const
{
    cols.push_back(wc_key(name, QCDOrder::LO, ContributionType::SM));
    cols.push_back(wc_key(name, QCDOrder::LO, ContributionType::BSM));
    if (wbc_.order > QCDOrder::LO) {
        cols.push_back(wc_key(name, QCDOrder::NLO, ContributionType::SM));
        cols.push_back(wc_key(name, QCDOrder::NLO, ContributionType::BSM));
        if (wbc_.order > QCDOrder::NLO) {
            cols.push_back(wc_key(name, QCDOrder::NNLO, ContributionType::SM));
            cols.push_back(wc_key(name, QCDOrder::NNLO, ContributionType::BSM));
        }
    }
}
