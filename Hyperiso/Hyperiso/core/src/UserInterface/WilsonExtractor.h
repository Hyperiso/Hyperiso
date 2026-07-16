#include "IOutputExtractor.h"
#include "WilsonInterface.h"
#include "Include.h"

std::string wc_key(const std::string& coefName, QCDOrder o, ContributionType c);

class WilsonCoeffExtractor : public IOutputExtractor {
public:
    WilsonCoeffExtractor(WilsonInterface& wi,
                         const WilsonBuildConfig& wbc,
                         const std::unordered_set<WCoefId>& coefs)
    : wi_(wi), wbc_(wbc), coefs_(coefs) {}

    std::vector<std::string> schema(const OutputSpec& /*spec*/) const override;

    void extract(std::unordered_map<std::string, Value>& outY,
                 const OutputSpec& /*spec*/) const override;

private:
    WilsonInterface& wi_;
    const WilsonBuildConfig& wbc_;
    const std::unordered_set<WCoefId>& coefs_;

    void add_one(std::unordered_map<std::string, Value>& outY,
                 WCoefId coef,
                 const std::string& name,
                 QCDOrder ord,
                 ContributionType ct) const;

    void push_cols_for_coef(std::vector<std::string>& cols,
                            const std::string& name) const;
};