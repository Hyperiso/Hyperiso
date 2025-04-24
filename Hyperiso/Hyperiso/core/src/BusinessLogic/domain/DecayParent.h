#ifndef DECAYPARENT_H
#define DECAYPARENT_H

#include <map>
#include <string>
#include "General.h"
#include "WilsonInterface.h"
#include "Node.h"
#include "ObsUseMarty.h"
#include "WilsonAdapter.h"
#include "ObsWilsonProxy.h"
#include "Math.h"

class DecayParent {

protected:
    std::map<Observables, std::shared_ptr<OperatorNode>> roots;
    QCDOrder max_order;
    ObsWilsonProxy w_proxy;

    QCDOrder check_max_order(QCDOrder order) const {
        if (order > max_order) {
            LOG_WARN("QCD order for decay cannot be higher than", OrderMapper::str(max_order));
            return max_order;
        }
        return order;
    }

public:
    explicit DecayParent() = default;
    void set_order(QCDOrder new_order);

    scalar_t compute_observable(Observables obs);
    size_t get_n_evals(Observables obs);

    virtual void build_op_tree() = 0;

};

#endif // __DECAYPARENT_H__