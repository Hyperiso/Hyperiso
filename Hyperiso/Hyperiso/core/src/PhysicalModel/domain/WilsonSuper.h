#ifndef WILSON_SUPER_H
#define WILSON_SUPER_H

#include "Include.h"
#include "Math.h"
#include "Utils.h"
#include "Wilson_parameters.h"
#include "ModelAPI.h"
#include "HasWilsonAPI.h"
#include "ParameterProxy.h"

struct MatchingInfo {
    std::unordered_set<ParamId> sources;
    std::function<double(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&)> compute;
    LhaID lhaid;
};

class WilsonCoefficient {
public:
    WilsonCoefficient() = default;
    WilsonCoefficient(const std::string& name, const std::string& storage_block) : coeffName(name), storage_block(storage_block) {}

    // TODO : Implement initialization as dependent parameter from MARTY library or from lha
    std::function<double(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&)> get_func(QCDOrder order);
    std::unordered_set<ParamId> get_sources(QCDOrder order);
    LhaID get_lhaid(QCDOrder order);
    
    void set_name(std::string name) {this->coeffName = name;}
    void set_owned(bool owned);
    void set_storage_block(std::string block_name);
    void set_contribution_type(ContributionType type);

    complex_t get_matching_value(std::string order) const; 
    std::string get_name() const {return this->coeffName;}
    std::string get_storage_block() const {return this->storage_block;}
    QCDOrder get_max_order() const {return this->max_order;}
    LhaID id(QCDOrder order) const;


    virtual void init_sources() = 0;
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