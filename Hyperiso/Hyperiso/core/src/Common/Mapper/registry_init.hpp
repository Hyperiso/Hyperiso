#pragma once
#include "observable_ids.hpp"
#include "wcoef_ids.hpp"
#include "wgroup_ids.hpp"
#include "parametertype_ids.hpp"
#include "model_ids.hpp"
#include "wilsonbasis_ids.hpp"
#include "contributiontype_ids.hpp"
#include "masstype_ids.hpp"
#include "scaletype_ids.hpp"
#include "uncertaintytype_ids.hpp"
#include "decay_ids.hpp"
#include "decay_graph.hpp"
#include "qcdorder_ids.hpp"

inline void init_decay_graph_builtins(){
    for (const auto& [d, vec] : decay_observable_mapping()){
        DecayId did = DecayMapper::to_id(d);
        for (auto o : vec){
            DecayGraph::instance().link(did, ObservableMapper::to_id(o));
        }
    }
}

inline void init_all_builtins() {
    static bool done=false; if (done) return; done=true;

    ObservableMapper::init_builtins();
    WCoefMapper::init_builtins();
    GroupMapper::init_builtins();
    ParameterTypeMapper::init_builtins();
    ModelMapper::init_builtins();
    WilsonBasisMapper::init_builtins();
    ContributionTypeMapper::init_builtins();
    MassTypeMapper::init_builtins();
    ScaleTypeMapper::init_builtins();
    UncertaintyTypeMapper::init_builtins();
    DecayMapper::init_builtins();
    OrderMapper::init_builtins();

    init_decay_graph_builtins();
}
