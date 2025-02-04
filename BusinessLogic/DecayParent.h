#ifndef __DECAYPARENT_H__
#define __DECAYPARENT_H__

#include <map>
#include <string>
#include "General.h"
#include "WilsonInterface.h"
#include "Node.h"
#include "QCDHelper.h"

struct WilsonInfo {
    double matching_scale;
    double hadronic_scale;
    Model model;
    QCDOrder order;
    BWilsonBasis basis;
    std::vector<WGroup> wgroups;
};

class DecayParent {

protected:
    std::map<Observables, std::shared_ptr<OperatorNode>> roots;
    WilsonInfo winfo;
    std::shared_ptr<WilsonInterface> wilson;
    QCDOrder max_order;

public:
    explicit DecayParent() = default;

    std::shared_ptr<WilsonInterface> get_wilsons(bool force_update=true);
    void set_order(QCDOrder new_order);

    scalar_t compute_observable(Observables obs);

    virtual void build_op_tree() = 0;

};

#endif // __DECAYPARENT_H__