#ifndef REGISTRY_INIT_H
#define REGISTRY_INIT_H

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
#include "decay_graph.h"
#include "qcdorder_ids.hpp"

/**
 * @file registry_init.h
 * @brief Central initialization helpers for all mappers and the decay graph.
 *
 * This header provides:
 *   - init_decay_graph_builtins(): populate DecayGraph from static mappings,
 *   - init_all_builtins(): one-shot initialization of all mapper registries
 *     and the decay graph.
 *
 * These functions are intended to be called once at startup or before
 * the first use of the mapping subsystem.
 */

/**
 * @brief Populates DecayGraph with builtin decay–observable relations.
 *
 * This iterates over decay_observable_mapping(), converts each decay and
 * observable enum into their respective ids via DecayMapper and
 * ObservableMapper, and registers the links in the DecayGraph singleton.
 */
inline void init_decay_graph_builtins(){
    for (const auto& [d, vec] : decay_observable_mapping()){
        DecayId did = DecayMapper::to_id(d);
        for (auto o : vec){
            DecayGraph::instance().link(did, ObservableMapper::to_id(o));
        }
    }
}

/**
 * @brief Initializes all builtin mapper registries and the decay graph.
 *
 * This function:
 *   - calls init_builtins() on all mapper types,
 *   - then calls init_decay_graph_builtins(),
 *   - is guarded by a static flag so it runs at most once.
 *
 * It is safe to call multiple times; only the first call performs
 * actual work.
 */
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

#endif