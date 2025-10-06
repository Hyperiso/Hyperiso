#ifndef WILSON_SUPER_H
#define WILSON_SUPER_H

#include "Include.h"
#include "Math.h"
#include "Utils.h"
#include "Parameter.h"
#include "IParamAdapter.h"
// #include "Wilson_parameters.h"
// #include "ModelAPI.h"
// #include "HasWilsonAPI.h"
// #include "ParameterProxy.h"

struct MatchingInfo {
    std::unordered_set<ParamId> sources = {};
    std::function<scalar_t(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&)> compute =
        [](const auto&) { return 0.0; };
    LhaID lhaid;

    MatchingInfo()
    : compute([](const auto&) { return 0.0; }) {}

    MatchingInfo(const LhaID& id) : lhaid(id) {}

    MatchingInfo(std::unordered_set<ParamId> src,
                 std::function<scalar_t(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&)> comp,
                 LhaID id)
        : sources(std::move(src)), compute(std::move(comp)), lhaid(std::move(id)) {}

    MatchingInfo(MatchingInfo&& mi) = default;
    MatchingInfo(const MatchingInfo& mi) = default;

    MatchingInfo operator=(const MatchingInfo& mi) {
        this->sources = mi.sources;
        this->compute = mi.compute;
        this->lhaid = mi.lhaid;
        return *this;
    }
};


//TODO : utility function
static bool ends_with(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

class WilsonCoefficient {
public:
    WilsonCoefficient() = default;
    WilsonCoefficient(const std::string& name, const std::string& storage_block) : coeffName(name), storage_block(storage_block) {
        if (ends_with(coeffName, "_THDM") || ends_with(coeffName, "_SUSY")) {
            type = ContributionType::BSM;
        }


        for (auto order : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
            matching_info[order] = MatchingInfo(this->id(order, type));
        }
    }

    WilsonCoefficient(const LhaID &name, const std::string& storage_block, ContributionType ct) : coeffName(name.to_string()), storage_block(storage_block) {
        type = ct;


        for (auto order : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
            matching_info[order] = MatchingInfo(name + LhaID((int)order, (int)ct));
        }
    }


    // TODO : Implement initialization as dependent parameter from MARTY library or from lha
    std::function<scalar_t(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&)> get_func(QCDOrder order);
    std::unordered_set<ParamId> get_sources(QCDOrder order);
    LhaID get_lhaid(QCDOrder order);
    std::string get_base_name() const;
    
    void set_name(std::string name) {this->coeffName = name;}
    void set_owned(bool owned);
    void set_storage_block(std::string block_name);
    void set_contribution_type(ContributionType type);

    complex_t get_matching_value(std::string order, ContributionType cont_type, std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_p) const; 
    std::string get_name() const {return this->coeffName;}
    std::string get_storage_block() const {return this->storage_block;}
    QCDOrder get_max_order() const {return this->max_order;}
    ContributionType get_type() {return this->type;}

    LhaID id(QCDOrder order, ContributionType typ) const;


    bool operator==(const WilsonCoefficient& other) const;
    bool operator!=(const WilsonCoefficient& other) const { return !(*this == other); }

    virtual ~WilsonCoefficient() = default;
    virtual std::shared_ptr<WilsonCoefficient> clone() const = 0;

protected:

    std::string coeffName{};
    QCDOrder max_order;
    ContributionType type {ContributionType::SM};
    bool is_owned {false};
    std::string storage_block;
    std::map<QCDOrder, MatchingInfo> matching_info;
};

// TODO : virtual LO, NLO, NNLO calculation need to be dealt properly



#endif //Wilsonv2