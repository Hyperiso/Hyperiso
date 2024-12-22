#ifndef __DECAYPARENT_H__
#define __DECAYPARENT_H__

#include <map>
#include <string>
#include "General.h"
#include "WilsonManager.h"
#include "Node.h"

struct WilsonInfo {
    double matching_scale;
    double hadronic_scale;
    Model model;
    QCDOrder order;
    BWilsonBasis basis;
    std::vector<WilsonGroups> wgroups;
};

class DecayParent {

protected:
    std::map<Observables, std::shared_ptr<OperatorNode>> roots;
    WilsonInfo winfo;

public:
    explicit DecayParent() = default;

    std::shared_ptr<CoefficientManager> compute_wilsons();

    scalar_t compute_observable(Observables obs);

    virtual void build_op_tree() = 0;

};

#endif // __DECAYPARENT_H__