#include "WilsonManager.h"
#include "WilsonBlockNames.h"
#include "WilsonRunningValidation.h"



static int qcd_index(QCDOrder o) {
    switch (o) {
        case QCDOrder::LO:   return 0;
        case QCDOrder::NLO:  return 1;
        case QCDOrder::NNLO: return 2;
        case QCDOrder::NONE: LOG_ERROR("ValueError", "Cannot handle QCDOrder::None as qcd index");
    }
    return 0;
}

static constexpr std::array<QCDOrder, 3> ALL_ORDERS = {
    QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO
};

static std::string hadronic_final_block_name(const std::string& groupName, WilsonBasis basis) {
    return GroupMapper::str(GroupMapper::enum_elt(groupName), ScaleType::HADRONIC, basis);
}

static std::string hadronic_sm_intermediate_block_name(const std::string& groupName, WilsonBasis basis) {
    return hadronic_final_block_name(groupName, basis) + "__SM_INTERMEDIATE";
}

static std::string hadronic_bsm_intermediate_block_name(const std::string& groupName, WilsonBasis basis) {
    return hadronic_final_block_name(groupName, basis) + "__BSM_INTERMEDIATE";
}

static scalar_t get_block_value_or_zero(
    const std::map<LhaID, std::shared_ptr<Parameter>>& items,
    const LhaID& id
) {
    auto it = items.find(id);
    if (it == items.end() || !it->second) {
        return 0.0;
    }
    return it->second->get_val();
}

static bool hard_coded_lo_enabled(const WilsonPortsConfig& pc) {
    return pc.hard_coded_lo && pc.hard_coded_lo->get();
}

static void compose_sm_only_triplet_for_order(
    const std::string& groupName,
    const std::vector<WCoefId>& members,
    QCDOrder order,
    WilsonPortsConfig& ports_config
) {
    WGroupId gid = GroupMapper::enum_elt(groupName);
    const std::string final_block = WilsonBlockNames::matching(gid);

    for (const auto& wcoef_id : members) {
        auto base = WCoefMapper::flha_base(wcoef_id);

        LhaID id_sm (base.first, base.second, qcd_index(order), 0);
        LhaID id_bsm(base.first, base.second, qcd_index(order), 1);
        LhaID id_tot(base.first, base.second, qcd_index(order), 2);

        ParamId pid_sm  { ParameterType::WILSON, final_block, id_sm  };
        ParamId pid_bsm { ParameterType::WILSON, final_block, id_bsm };
        ParamId pid_tot { ParameterType::WILSON, final_block, id_tot };

        ports_config.iblock_c->compose_parameter(
            pid_bsm,
            std::unordered_set<ParamId>{},
            [](const ParamSrc&, std::shared_ptr<DependentParameter> dep) {
                dep->set_expected(0.);
            }
        );

        ports_config.iblock_c->compose_parameter(
            pid_tot,
            std::unordered_set<ParamId>{ pid_sm },
            [pid_sm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                dep->set_expected(src.get_val(pid_sm));
            }
        );
    }
}

static std::vector<QCDOrder> orders_up_to(QCDOrder max) {
    std::vector<QCDOrder> out{QCDOrder::LO};
    if (max >= QCDOrder::NLO)  out.push_back(QCDOrder::NLO);
    if (max >= QCDOrder::NNLO) out.push_back(QCDOrder::NNLO);
    return out;
}

void CoefficientManager::throw_no_group_error(const std::string &groupName) const {
    std::stringstream ss;
    ss << "Coefficient group " << groupName << " not found in manager. Existing groups are:\n";
    for (const auto& group : this->coefficientGroups) {
        ss << "\t- " << group.first << "\n";
    }
    LOG_ERROR("KeyError", ss.str());
}

void CoefficientManager::ensure_final_triplet_defaults_zero(
    const std::string& groupName,
    QCDOrder max_order
) {
    WGroupId gid = GroupMapper::enum_elt(groupName);
    const std::string final_block = WilsonBlockNames::matching(gid);

    const auto& members = this->coefficientGroups.at(groupName)->get_member_ids();

    auto wp = ports_config.wilson_proxy;

    for (auto o : orders_up_to(max_order)) {
        for (const auto& wcoef_id : members) {
            auto base = WCoefMapper::flha_base(wcoef_id);

            LhaID id_sm (base.first, base.second, qcd_index(o), 0);
            LhaID id_bsm(base.first, base.second, qcd_index(o), 1);
            LhaID id_tot(base.first, base.second, qcd_index(o), 2);

            ParamId pid_sm  { ParameterType::WILSON, final_block, id_sm  };
            ParamId pid_bsm { ParameterType::WILSON, final_block, id_bsm };
            ParamId pid_tot { ParameterType::WILSON, final_block, id_tot };

            auto zero = [](const ParamSrc&, std::shared_ptr<DependentParameter> dep) {
                dep->set_expected(0.);
            };

            if (!wp->exist(final_block, id_sm))  ports_config.iblock_c->compose_parameter(pid_sm,  {}, zero);
            if (!wp->exist(final_block, id_bsm)) ports_config.iblock_c->compose_parameter(pid_bsm, {}, zero);
            if (!wp->exist(final_block, id_tot)) ports_config.iblock_c->compose_parameter(pid_tot, {}, zero);
        }
    }
}

// void CoefficientManager::ensure_sm_intermediate_and_copy_to_final(
//     const std::string& groupName,
//     QCDOrder max_order,
//     WilsonPortsConfig& ports_config
// ) {
//     WGroupId gid = GroupMapper::enum_elt(groupName);
//     const std::string sm_block    = WilsonBlockNames::sm_matching(gid);
//     const std::string final_block = WilsonBlockNames::matching(gid);

//     std::shared_ptr<CoefficientGroup> sm_group_ptr;

//     if (ports_config.build_group) {
//         const bool marty_backend = ports_config.use_marty->get();
//         sm_group_ptr = ports_config.build_group(
//             gid,
//             Model::SM,
//             marty_backend,
//             ContributionType::SM,
//             sm_block
//         );
//     } else {
//         sm_group_ptr = this->coefficientGroups.at(groupName)->get_sm_group();
//         sm_group_ptr->set_matching_storage_block(sm_block);
//     }

//     if (!sm_group_ptr) {
//         LOG_ERROR("LogicError", "No SM group found for " + groupName);
//     }

//     if (!this->coefficientGroups.contains(sm_block)) {
//         this->registerCoefficientGroup(sm_block, sm_group_ptr);
//     }

//     // IMPORTANT: init au max => compose LO+NLO(+NNLO) côté SM
//     sm_group_ptr->init(max_order);

//     // Copie SM (par ordre) : sm_block -> final_block
//     auto maybe_g = GroupMapper::enum_of(gid);
//     if (!maybe_g) LOG_ERROR("LogicError", "Bad group id for " + groupName);

//     for (auto o : ALL_ORDERS) {
//         if (qcd_index(o) > qcd_index(max_order)) continue;

//         for (auto wcoef_enum : WCoefMapper::get_group(*maybe_g)) {
//             auto base = WCoefMapper::flha_base(wcoef_enum);

//             LhaID id_sm(base.first, base.second, qcd_index(o), 0);

//             ParamId pid_src { ParameterType::WILSON, sm_block,    id_sm };
//             ParamId pid_dst { ParameterType::WILSON, final_block, id_sm };

//             ports_config.iblock_c->compose_parameter(
//                 pid_dst,
//                 std::unordered_set<ParamId>{ pid_src },
//                 [pid_src](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
//                     dep->set_expected(src.get_val(pid_src));
//                 }
//             );
//         }
//     }
// }

void CoefficientManager::ensure_sm_intermediate_and_copy_to_final(
    const std::string& groupName,
    QCDOrder max_order,
    WilsonPortsConfig& ports_config
) {
    WGroupId gid = GroupMapper::enum_elt(groupName);
    const std::string sm_block    = WilsonBlockNames::sm_matching(gid);
    const std::string final_block = WilsonBlockNames::matching(gid);

    std::shared_ptr<CoefficientGroup> sm_group_ptr;

    const bool marty_backend = ports_config.use_marty->get();
    const bool hard_lo = hard_coded_lo_enabled(ports_config);

    bool allow_hardcoded_sm = (marty_backend && hard_lo && (max_order > QCDOrder::LO));
    if (allow_hardcoded_sm && !ports_config.build_group) {
        LOG_WARN(
            "(CoefficientManager) hard_coded_lo=true but no ports_config.build_group provided; "
            "cannot build a hardcoded SM group. Falling back to LO-only SM."
        );
        allow_hardcoded_sm = false;
    }

    const bool sm_use_marty = allow_hardcoded_sm ? false : marty_backend;

    const QCDOrder sm_max_order = (marty_backend && (max_order > QCDOrder::LO) && !allow_hardcoded_sm)
        ? QCDOrder::LO
        : max_order;

    if (ports_config.build_group) {
        sm_group_ptr = ports_config.build_group(
            gid,
            Model::SM,
            sm_use_marty,
            ContributionType::SM,
            sm_block
        );
    } else {
        sm_group_ptr = this->coefficientGroups.at(groupName)->get_sm_group();
        sm_group_ptr->set_matching_storage_block(sm_block);
    }

    if (!sm_group_ptr) {
        LOG_ERROR("LogicError", "No SM group found for " + groupName);
    }

    if (!this->coefficientGroups.contains(sm_block)) {
        this->registerCoefficientGroup(sm_block, sm_group_ptr);
    }

    sm_group_ptr->init(sm_max_order);

    const auto& members = this->coefficientGroups.at(groupName)->get_member_ids();

    for (auto o : ALL_ORDERS) {
        if (qcd_index(o) > qcd_index(sm_max_order)) continue;

        for (const auto& wcoef_id : members) {
            auto base = WCoefMapper::flha_base(wcoef_id);

            LhaID id_sm(base.first, base.second, qcd_index(o), 0);

            ParamId pid_src { ParameterType::WILSON, sm_block,    id_sm };
            ParamId pid_dst { ParameterType::WILSON, final_block, id_sm };

            ports_config.iblock_c->compose_parameter(
                pid_dst,
                std::unordered_set<ParamId>{ pid_src },
                [pid_src](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                    dep->set_expected(src.get_val(pid_src));
                }
            );
        }
    }
}

std::string CoefficientManager::getModel() {
    return ModelMapper::str(ports_config.model_api->get());
}

void CoefficientManager::init_group_matching(const std::string& groupName, const std::string& order) {
    init_specific_order_group_matching(groupName, order, /*only_total=*/false);

}
void CoefficientManager::ensure_matching_triplet_zeroed(
    const std::string& groupName,
    QCDOrder o
) {
    std::string storage_block = this->coefficientGroups.at(groupName)->get_matching_storage_block();

    for (auto& coeff : *this->coefficientGroups.at(groupName)) {
        auto base = WCoefMapper::enum_elt(coeff.second->get_base_name());

        ParamId pid_SM  { storage_block, WCoefMapper::flha_full(base, o, ContributionType::SM) };
        ParamId pid_BSM { storage_block, WCoefMapper::flha_full(base, o, ContributionType::BSM) };
        ParamId pid_TOT { storage_block, WCoefMapper::flha_full(base, o, ContributionType::TOTAL) };

        auto zero = [] (const ParamSrc&, std::shared_ptr<DependentParameter> dep) {
            dep->set_expected(0.);
        };

        ports_config.iblock_c->compose_parameter(pid_SM,  {}, zero);
        ports_config.iblock_c->compose_parameter(pid_BSM, {}, zero);
        ports_config.iblock_c->compose_parameter(pid_TOT, {}, zero);
    }
}

void CoefficientManager::fill_matching_groups(const std::string& groupName, const std::string& orderStr) {
    QCDOrder max = OrderMapper::enum_elt(orderStr);

    if (max == QCDOrder::LO) {
        ensure_matching_triplet_zeroed(groupName, QCDOrder::NLO);
        ensure_matching_triplet_zeroed(groupName, QCDOrder::NNLO);
    } else if (max == QCDOrder::NLO) {
        ensure_matching_triplet_zeroed(groupName, QCDOrder::NNLO);
    } else {
    }
}


void CoefficientManager::compose_from_fwcoef(
    const std::string& groupName,
    QCDOrder max_order,
    WilsonPortsConfig& ports_config
) {
    WGroupId gid = GroupMapper::enum_elt(groupName);
    const std::string sm_block    = WilsonBlockNames::sm_matching(gid);
    const std::string final_block = WilsonBlockNames::matching(gid);
    const std::string fw_block    = WilsonBlockNames::fwcoef();

    auto wp = ports_config.wilson_proxy;

    const auto& members = this->coefficientGroups.at(groupName)->get_member_ids();

    for (auto o : ALL_ORDERS) {
        if (qcd_index(o) > qcd_index(max_order)) continue;

        for (const auto& wcoef_id : members) {
            auto base = WCoefMapper::flha_base(wcoef_id);

            LhaID id_sm (base.first, base.second, qcd_index(o), 0);
            LhaID id_bsm(base.first, base.second, qcd_index(o), 1);
            LhaID id_tot(base.first, base.second, qcd_index(o), 2);

            ParamId fw_sm  { ParameterType::WILSON, fw_block, id_sm  };
            ParamId fw_bsm { ParameterType::WILSON, fw_block, id_bsm };
            ParamId fw_tot { ParameterType::WILSON, fw_block, id_tot };

            const bool fw_has_sm  = wp->exist(fw_block, id_sm);
            const bool fw_has_bsm = wp->exist(fw_block, id_bsm);
            const bool fw_has_tot = wp->exist(fw_block, id_tot);

            ParamId sm_inter { ParameterType::WILSON, sm_block, id_sm };

            if (fw_has_sm) {
                ports_config.iblock_c->compose_parameter(
                    sm_inter,
                    std::unordered_set<ParamId>{ fw_sm },
                    [fw_sm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                        dep->set_expected(src.get_val(fw_sm));
                    }
                );
            } else if (fw_has_tot && fw_has_bsm) {
                ports_config.iblock_c->compose_parameter(
                    sm_inter,
                    std::unordered_set<ParamId>{ fw_tot, fw_bsm },
                    [fw_tot, fw_bsm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                        dep->set_expected(src.get_val(fw_tot) - src.get_val(fw_bsm));
                    }
                );
            }

            ParamId final_sm  { ParameterType::WILSON, final_block, id_sm  };
            ParamId final_bsm { ParameterType::WILSON, final_block, id_bsm };
            ParamId final_tot { ParameterType::WILSON, final_block, id_tot };

            if (fw_has_bsm) {
                ports_config.iblock_c->compose_parameter(
                    final_bsm,
                    std::unordered_set<ParamId>{ fw_bsm },
                    [fw_bsm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                        dep->set_expected(src.get_val(fw_bsm));
                    }
                );
            } else if (fw_has_tot) {
                ports_config.iblock_c->compose_parameter(
                    final_bsm,
                    std::unordered_set<ParamId>{ fw_tot, final_sm },
                    [fw_tot, final_sm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                        dep->set_expected(src.get_val(fw_tot) - src.get_val(final_sm));
                    }
                );
            } else {
                ports_config.iblock_c->compose_parameter(
                    final_bsm,
                    std::unordered_set<ParamId>{},
                    [](const ParamSrc&, std::shared_ptr<DependentParameter> dep) {
                        dep->set_expected(0.0);
                    }
                );
            }

            if (fw_has_tot) {
                ports_config.iblock_c->compose_parameter(
                    final_tot,
                    std::unordered_set<ParamId>{ fw_tot },
                    [fw_tot](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                        dep->set_expected(src.get_val(fw_tot));
                    }
                );
            } else if (fw_has_bsm) {
                ports_config.iblock_c->compose_parameter(
                    final_tot,
                    std::unordered_set<ParamId>{ final_sm, final_bsm },
                    [final_sm, final_bsm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                        dep->set_expected(src.get_val(final_sm) + src.get_val(final_bsm));
                    }
                );
            } else {
                ports_config.iblock_c->compose_parameter(
                    final_tot,
                    std::unordered_set<ParamId>{ final_sm },
                    [final_sm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                        dep->set_expected(src.get_val(final_sm));
                    }
                );
            }
        }
    }
}

void CoefficientManager::compose_missing_from_calculation(
    const std::string& groupName,
    QCDOrder order,
    WilsonPortsConfig& ports_config
) {
    WGroupId gid = GroupMapper::enum_elt(groupName);
    const std::string final_block = WilsonBlockNames::matching(gid);
    const bool marty_backend = ports_config.use_marty->get();

    const auto& members = this->coefficientGroups.at(groupName)->get_member_ids();

    for (const auto& wcoef_id : members) {
        auto base = WCoefMapper::flha_base(wcoef_id);

        LhaID id_sm (base.first, base.second, qcd_index(order), 0);
        LhaID id_bsm(base.first, base.second, qcd_index(order), 1);
        LhaID id_tot(base.first, base.second, qcd_index(order), 2);

        ParamId pid_sm  { ParameterType::WILSON, final_block, id_sm  };
        ParamId pid_bsm { ParameterType::WILSON, final_block, id_bsm };
        ParamId pid_tot { ParameterType::WILSON, final_block, id_tot };

        if (marty_backend) {
            ports_config.iblock_c->compose_parameter(
                pid_bsm,
                std::unordered_set<ParamId>{ pid_tot, pid_sm },
                [pid_tot, pid_sm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                    dep->set_expected(src.get_val(pid_tot) - src.get_val(pid_sm));
                }
            );
        } else {
            ports_config.iblock_c->compose_parameter(
                pid_tot,
                std::unordered_set<ParamId>{ pid_sm, pid_bsm },
                [pid_sm, pid_bsm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                    dep->set_expected(src.get_val(pid_sm) + src.get_val(pid_bsm));
                }
            );
        }
    }
}

// void CoefficientManager::ensure_sm_model_triplet_in_matching(
//     const std::string& groupName,
//     QCDOrder max_order
// ) {
//     WGroupId gid = GroupMapper::enum_elt(groupName);
//     const std::string final_block = WilsonBlockNames::matching(gid);

//     auto maybe_g = GroupMapper::enum_of(gid);
//     if (!maybe_g) LOG_ERROR("LogicError", "Bad group id for " + groupName);
//     const auto members = WCoefMapper::get_group(*maybe_g);

//     for (auto o : orders_up_to(max_order)) {
//         for (auto wcoef_enum : members) {
//             auto base = WCoefMapper::flha_base(wcoef_enum);

//             LhaID id_sm (base.first, base.second, qcd_index(o), 0);
//             LhaID id_bsm(base.first, base.second, qcd_index(o), 1);
//             LhaID id_tot(base.first, base.second, qcd_index(o), 2);

//             ParamId pid_sm  { ParameterType::WILSON, final_block, id_sm  };
//             ParamId pid_bsm { ParameterType::WILSON, final_block, id_bsm };
//             ParamId pid_tot { ParameterType::WILSON, final_block, id_tot };

//             ports_config.iblock_c->compose_parameter(
//                 pid_bsm,
//                 std::unordered_set<ParamId>{},
//                 [](const ParamSrc&, std::shared_ptr<DependentParameter> dep) {
//                     dep->set_expected(0.);
//                 }
//             );

//             ports_config.iblock_c->compose_parameter(
//                 pid_tot,
//                 std::unordered_set<ParamId>{ pid_sm },
//                 [pid_sm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
//                     dep->set_expected(src.get_val(pid_sm));
//                 }
//             );
//         }
//     }
// }

void CoefficientManager::ensure_sm_model_triplet_in_matching(
    const std::string& groupName,
    QCDOrder max_order
) {
    WGroupId gid = GroupMapper::enum_elt(groupName);
    const std::string final_block = WilsonBlockNames::matching(gid);

    const auto& members = this->coefficientGroups.at(groupName)->get_member_ids();

    for (auto o : orders_up_to(max_order)) {
        for (const auto& wcoef_id : members) {
            auto base = WCoefMapper::flha_base(wcoef_id);

            LhaID id_sm (base.first, base.second, qcd_index(o), 0);
            LhaID id_bsm(base.first, base.second, qcd_index(o), 1);
            LhaID id_tot(base.first, base.second, qcd_index(o), 2);

            ParamId pid_sm  { ParameterType::WILSON, final_block, id_sm  };
            ParamId pid_bsm { ParameterType::WILSON, final_block, id_bsm };
            ParamId pid_tot { ParameterType::WILSON, final_block, id_tot };

            ports_config.iblock_c->compose_parameter(
                pid_bsm,
                std::unordered_set<ParamId>{},
                [](const ParamSrc&, std::shared_ptr<DependentParameter> dep) {
                    dep->set_expected(0.0);
                }
            );

            ports_config.iblock_c->compose_parameter(
                pid_tot,
                std::unordered_set<ParamId>{ pid_sm, pid_bsm },
                [pid_sm, pid_bsm](const ParamSrc& src, std::shared_ptr<DependentParameter> dep) {
                    dep->set_expected(src.get_val(pid_sm) + src.get_val(pid_bsm));
                }
            );
        }
    }
}

// void CoefficientManager::init_specific_order_group_matching(const std::string& groupName,
//                                                             const std::string& orderStr,
//                                                             bool only_total)
// {
//     if (!this->coefficientGroups.contains(groupName)) {
//         throw_no_group_error(groupName);
//     }

//     const bool has_input = ports_config.has_wilson->get();
//     // const bool marty = ports_config.use_marty->get();
//     const bool SM_model = (ports_config.model_api->get() == Model::SM);

//     QCDOrder order = OrderMapper::enum_elt(orderStr);

//     if (!has_input) {
//         if (!only_total) {
//             this->coefficientGroups.at(groupName)->init(order);
//         }

//         if (SM_model) {
//             ensure_sm_model_triplet_in_matching(groupName, order);
//             return;
//         }

//         ensure_sm_intermediate_and_copy_to_final(groupName, order, ports_config);
//         if (!only_total) compose_missing_from_calculation(groupName, order, ports_config);

//     } else {
//         ensure_sm_intermediate_and_copy_to_final(groupName, order, ports_config);

//         compose_from_fwcoef(groupName, order, ports_config);

//     }
// }

void CoefficientManager::init_specific_order_group_matching(const std::string& groupName,
                                                            const std::string& orderStr,
                                                            bool only_total)
{
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    const bool has_input     = ports_config.has_wilson->get();
    const bool marty_backend = ports_config.use_marty->get();
    const bool hard_lo       = hard_coded_lo_enabled(ports_config);
    const bool SM_model      = (ports_config.model_api->get() == Model::SM);

    const QCDOrder requested = OrderMapper::enum_elt(orderStr);

    const bool mixed_orders = marty_backend
        && hard_lo
        && (requested > QCDOrder::LO)
        && static_cast<bool>(ports_config.build_group);

    const QCDOrder marty_calc_order = (marty_backend && (requested > QCDOrder::LO))
        ? QCDOrder::LO
        : requested;

    if (!has_input) {
        if (!only_total) {
            this->coefficientGroups.at(groupName)->init(marty_calc_order);
        }

        if (SM_model) {
            ensure_sm_model_triplet_in_matching(groupName, marty_calc_order);

            for (auto o : ALL_ORDERS) {
                if (qcd_index(o) > qcd_index(requested)) continue;
                if (qcd_index(o) <= qcd_index(marty_calc_order)) continue;
                ensure_matching_triplet_zeroed(groupName, o);
            }
            return;
        }

        ensure_sm_intermediate_and_copy_to_final(groupName, requested, ports_config);

        if (!only_total) {
            if (marty_backend && (requested > QCDOrder::LO)) {
                compose_missing_from_calculation(groupName, QCDOrder::LO, ports_config);

                if (mixed_orders) {
                    // >LO : BSM=0 ; TOTAL=SM
                    for (auto o : ALL_ORDERS) {
                        if (o == QCDOrder::LO) continue;
                        if (qcd_index(o) > qcd_index(requested)) continue;
                        compose_sm_only_triplet_for_order(groupName, this->coefficientGroups.at(groupName)->get_member_ids(), o, ports_config);
                    }
                } else {
                    for (auto o : ALL_ORDERS) {
                        if (o == QCDOrder::LO) continue;
                        if (qcd_index(o) > qcd_index(requested)) continue;
                        ensure_matching_triplet_zeroed(groupName, o);
                    }

                    if (hard_lo && !ports_config.build_group) {
                        LOG_WARN(
                            "(CoefficientManager) hard_coded_lo=true requested but ports_config.build_group is not set; "
                            "cannot compute SM beyond LO. Higher orders are zero-filled."
                        );
                    } else if (!hard_lo) {
                        LOG_WARN(
                            "(CoefficientManager) Marty backend does not support QCD orders beyond LO; "
                            "higher orders are zero-filled. Set hard_coded_lo=true (with build_group hook) to keep SM up to the requested order."
                        );
                    }
                }
            } else {
                compose_missing_from_calculation(groupName, requested, ports_config);
            }
        }

    } else {
        ensure_sm_intermediate_and_copy_to_final(groupName, requested, ports_config);
        compose_from_fwcoef(groupName, requested, ports_config);
    }
}

void CoefficientManager::fill_sources_for_group(const std::string & groupName, const std::string& order, std::unordered_map<ParameterType, std::vector<std::string>>& src, WilsonBasis id) {
    switch (OrderMapper::enum_elt(order))
    {
    case QCDOrder::NNLO:
        for (const auto& [key, vec] : this->coefficientGroups[groupName]->get_sources(QCDOrder::NNLO, id)) {
            auto& targetVec = src[key];
            targetVec.insert(targetVec.end(), vec.begin(), vec.end());
        }
        [[fallthrough]];
        
    case QCDOrder::NLO:
        for (const auto& [key, vec] : this->coefficientGroups[groupName]->get_sources(QCDOrder::NLO, id)) {
            auto& targetVec = src[key];
            targetVec.insert(targetVec.end(), vec.begin(), vec.end());
        }
        [[fallthrough]];

    case QCDOrder::LO:
        for (const auto& [key, vec] : this->coefficientGroups[groupName]->get_sources(QCDOrder::LO, id)) {
            auto& targetVec = src[key];
            targetVec.insert(targetVec.end(), vec.begin(), vec.end());
        }
        break;

    default:
        break;
    }
}

std::pair<WCoefId, std::pair<QCDOrder, ContributionType>> lha_wilson_deserialize(LhaID id) {
    auto parts = id.get_parts();
    auto w_id = std::make_pair<int, int>(parts[0], parts[1]);

    auto maybe = WCoefMapper::from_flha_key(w_id.first, w_id.second);
    if (!maybe) {
        LOG_ERROR("ValueError", "bad lha id for wilson conversion (unknown custom/base key)");
    }
    WCoefId coef = *maybe;
    QCDOrder order = parts[2] ? ((parts[2] -1) ? QCDOrder::NNLO : QCDOrder::NLO) : QCDOrder::LO;
    ContributionType part = parts[3] ? parts[3] -1 ? ContributionType::TOTAL : ContributionType::BSM : ContributionType::SM;

    std::pair<WCoefId, std::pair<QCDOrder, ContributionType>> ret;
    ret = {coef, {order, part}};

    return ret;
}

void CoefficientManager::init_group_hadronic(const std::string& groupName,
                                            const std::string& order,
                                            WilsonBasis basis) {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    std::unordered_map<ParameterType, std::vector<std::string>> running_src = {};
    fill_sources_for_group(groupName, order, running_src, basis);

    const QCDOrder ord = OrderMapper::enum_elt(order);

    const std::map<QCDOrder,
        std::function<std::unordered_map<WCoefId, scalar_t>(
            const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>&,
            const BlockSrc&)>> funcs = {
        {QCDOrder::LO,   this->coefficientGroups[groupName]->get_func(QCDOrder::LO, basis)},
        {QCDOrder::NLO,  this->coefficientGroups[groupName]->get_func(QCDOrder::NLO, basis)},
        {QCDOrder::NNLO, this->coefficientGroups[groupName]->get_func(QCDOrder::NNLO, basis)}
    };

    const std::string matching_block_name = this->coefficientGroups[groupName]->get_matching_storage_block();
    const auto& members = this->coefficientGroups.at(groupName)->get_member_ids();

    const std::string final_block_name = hadronic_final_block_name(groupName, basis);
    const std::string sm_run_block_name = hadronic_sm_intermediate_block_name(groupName, basis);
    const std::string bsm_run_block_name = hadronic_bsm_intermediate_block_name(groupName, basis);

    auto make_partial_running_func =
        [matching_block_name, ord, funcs, members, groupName, basis]
        (ContributionType kept_contrib, const std::string& target_block_name) {

        return [matching_block_name, ord, funcs, kept_contrib, target_block_name, members, groupName, basis]
               (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
            auto raw = src.raw();
            auto it = raw.find(matching_block_name);
            if (it == raw.end() || !it->second) {
                return;
            }

            const auto& matching_coeff = it->second->getItems();

            std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>> matching_one;

            // Defensive defaulting: built-in running functions often use .at(order).at(coef).
            // Missing BSM coefficients, missing higher orders, or custom sparse groups should
            // therefore be represented explicitly by zero instead of leaving the key absent.
            for (auto qcd_order : ALL_ORDERS) {
                if (qcd_index(qcd_order) > qcd_index(ord)) {
                    continue;
                }
                for (const auto& wcoef_id : members) {
                    matching_one[qcd_order][wcoef_id] = 0.0;
                }
            }

            for (const auto& [lha_id, param] : matching_coeff) {
                auto dec = lha_wilson_deserialize(lha_id);
                const WCoefId& wcoef = dec.first;
                const QCDOrder& qcd_order = dec.second.first;
                const ContributionType& contrib = dec.second.second;

                if (qcd_index(qcd_order) > qcd_index(ord) || contrib != kept_contrib || !param) {
                    continue;
                }

                matching_one[qcd_order][wcoef] = param->get_val();
            }

            std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>> res_one;

            auto run_order = [&](QCDOrder qcd_order) {
                WilsonRunningValidation::require_running_function(funcs, qcd_order, groupName, basis);
                WilsonRunningValidation::require_matching_input_complete(matching_one, members, qcd_order, groupName, basis);

                auto result = funcs.at(qcd_order)(matching_one, src);
                WilsonRunningValidation::require_running_result_known_members(result, members, qcd_order, groupName, basis);
                res_one[qcd_order] = std::move(result);
            };

            switch (ord) {
                case QCDOrder::NNLO:
                    run_order(QCDOrder::NNLO);
                    [[fallthrough]];
                case QCDOrder::NLO:
                    run_order(QCDOrder::NLO);
                    [[fallthrough]];
                case QCDOrder::LO:
                    run_order(QCDOrder::LO);
                    break;
                default:
                    break;
            }

            for (const auto& [qcd_order, coef_map] : res_one) {
                for (const auto& [coef_id, coef_val] : coef_map) {
                    const LhaID coef_lha = WCoefMapper::flha_full(coef_id, qcd_order, kept_contrib);
                    const ParamId pid{
                        ParameterType::WILSON,
                        target_block_name,
                        coef_lha
                    };
                    dep_block->store_or_assign(
                        pid.code,
                        std::make_shared<Parameter>(pid, coef_val, 0.0, static_cast<int>(kept_contrib))
                    );
                }
            }
        };
    };

    auto sm_func  = make_partial_running_func(ContributionType::SM,  sm_run_block_name);
    auto bsm_func = make_partial_running_func(ContributionType::BSM, bsm_run_block_name);

    ports_config.iblock_c->compose_block(sm_run_block_name, running_src, sm_func);
    ports_config.iblock_c->compose_block(bsm_run_block_name, running_src, bsm_func);

    std::unordered_map<ParameterType, std::vector<std::string>> combine_src = {
        { ParameterType::WILSON, { sm_run_block_name, bsm_run_block_name } }
    };


    auto combine_func =
        [groupName, basis, ord, final_block_name, sm_run_block_name, bsm_run_block_name, members]
        (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
            auto raw = src.raw();

            auto it_sm  = raw.find(sm_run_block_name);
            auto it_bsm = raw.find(bsm_run_block_name);

            if (it_sm == raw.end() || !it_sm->second) {
                return;
            }
            if (it_bsm == raw.end() || !it_bsm->second) {
                return;
            }

            const auto& sm_items  = it_sm->second->getItems();
            const auto& bsm_items = it_bsm->second->getItems();

            for (auto qcd_order : ALL_ORDERS) {
                if (qcd_index(qcd_order) > qcd_index(ord)) {
                    continue;
                }

                for (const auto& wcoef_id : members) {
                    auto base = WCoefMapper::flha_base(wcoef_id);

                    const LhaID id_sm (base.first, base.second, qcd_index(qcd_order), 0);
                    const LhaID id_bsm(base.first, base.second, qcd_index(qcd_order), 1);
                    const LhaID id_tot(base.first, base.second, qcd_index(qcd_order), 2);

                    const scalar_t sm_val  = get_block_value_or_zero(sm_items,  id_sm);
                    const scalar_t bsm_val = get_block_value_or_zero(bsm_items, id_bsm);
                    const scalar_t tot_val = sm_val + bsm_val;

                    {
                        const ParamId pid{ParameterType::WILSON, final_block_name, id_sm};
                        dep_block->store_or_assign(
                            pid.code,
                            std::make_shared<Parameter>(pid, sm_val, 0.0, static_cast<int>(ContributionType::SM))
                        );
                    }
                    {
                        const ParamId pid{ParameterType::WILSON, final_block_name, id_bsm};
                        dep_block->store_or_assign(
                            pid.code,
                            std::make_shared<Parameter>(pid, bsm_val, 0.0, static_cast<int>(ContributionType::BSM))
                        );
                    }
                    {
                        const ParamId pid{ParameterType::WILSON, final_block_name, id_tot};
                        dep_block->store_or_assign(
                            pid.code,
                            std::make_shared<Parameter>(pid, tot_val, 0.0, static_cast<int>(ContributionType::TOTAL))
                        );
                    }
                }
            }
        };

    ports_config.iblock_c->compose_block(final_block_name, combine_src, combine_func);
}

// void CoefficientManager::init_group_hadronic(const std::string& groupName, const std::string& order, WilsonBasis basis) {

    

//     if (!this->coefficientGroups.contains(groupName)) {
//         throw_no_group_error(groupName);
//     }
//     std::unordered_map<ParameterType, std::vector<std::string>> src = {};
//     fill_sources_for_group(groupName, order, src, basis);
//     QCDOrder ord = OrderMapper::enum_elt(order);
//     std::map<QCDOrder, std::function<std::unordered_map<WCoefId, scalar_t>(const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>&, const BlockSrc&)>> funcs = {
//         {QCDOrder::LO, this->coefficientGroups[groupName]->get_func(QCDOrder::LO, basis)},
//         {QCDOrder::NLO, this->coefficientGroups[groupName]->get_func(QCDOrder::NLO, basis)},
//         {QCDOrder::NNLO, this->coefficientGroups[groupName]->get_func(QCDOrder::NNLO, basis)}
//     };

//     std::string matching_block_name = this->coefficientGroups[groupName]->get_matching_storage_block();
//     auto func = [matching_block_name, ord, funcs, groupName, basis] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
//         auto raw = src.raw();
//         auto it = raw.find(matching_block_name);
//         if (it == raw.end() || !it->second) {
//             return;
//         }

//         auto maybe_dump_coef = [&](const char* tag,
//                             const std::unordered_map<ContributionType,
//                             std::unordered_map<QCDOrder,
//                             std::unordered_map<WCoefId, scalar_t>>>& M)
//     {
//         auto dump_one = [&](ContributionType ct, QCDOrder qo, WCoef wc, const char* name) {
//             auto it1 = M.find(ct);
//             if (it1 == M.end()) return;
//             auto it2 = it1->second.find(qo);
//             if (it2 == it1->second.end()) return;
//             auto it3 = it2->second.find(WCoefMapper::to_id(wc));
//             if (it3 == it2->second.end()) return;
//             std::cout << "[HADDBG] " << tag
//                     << " contrib=" << int(ct)
//                     << " order=" << int(qo)
//                     << " " << name
//                     << " = " << it3->second << "\n";
//         };

//         dump_one(ContributionType::SM,    QCDOrder::LO, WCoef::C9,  "C9");
//         dump_one(ContributionType::BSM,   QCDOrder::LO, WCoef::C9,  "C9");
//         dump_one(ContributionType::TOTAL, QCDOrder::LO, WCoef::C9,  "C9");
//         dump_one(ContributionType::SM,    QCDOrder::LO, WCoef::C10, "C10");
//         dump_one(ContributionType::BSM,   QCDOrder::LO, WCoef::C10, "C10");
//         dump_one(ContributionType::TOTAL, QCDOrder::LO, WCoef::C10, "C10");
//     };
    
//         auto matching_coeff = it->second->getItems();
//         std::unordered_map<ContributionType, std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>> matching_map;
//         for (auto& coef : matching_coeff) {
//             std::pair<WCoefId, std::pair<QCDOrder, ContributionType>> c = lha_wilson_deserialize(coef.first);
//             const WCoefId& wcoef = c.first;
//             const QCDOrder& order = c.second.first;
//             const ContributionType& contrib = c.second.second;
//             matching_map[contrib][order][wcoef] = coef.second->get_val();
//         }

//         maybe_dump_coef("matching_map", matching_map);
//         std::unordered_map<ContributionType, std::unordered_map<QCDOrder,std::unordered_map<WCoefId, scalar_t>>> res;
//         for (auto contri : {ContributionType::SM, ContributionType::BSM, ContributionType::TOTAL}) {
//             switch (ord)
//                 {
//                 case QCDOrder::NNLO:
//                     res[contri][QCDOrder::NNLO] = funcs.at(QCDOrder::NNLO)(matching_map[contri], src);
//                     [[fallthrough]];
                    
//                 case QCDOrder::NLO:
//                     res[contri][QCDOrder::NLO] = funcs.at(QCDOrder::NLO)(matching_map[contri], src);
//                     [[fallthrough]];

//                 case QCDOrder::LO:
//                     res[contri][QCDOrder::LO] = funcs.at(QCDOrder::LO)(matching_map[contri], src);
//                     break;

//                 default:
//                     break;
//                 }
//         }
//         maybe_dump_coef("res", res);
//         for (auto& [c_type, order_map] : res) { // Iterate over the contributions
//             for (auto& [order, coef_map] : order_map) { // Iterate over the orders
//                 for (auto& [coef_id, coef_val] : coef_map) { // Iterate over the coefficients
//                     LhaID coef_lha = WCoefMapper::flha_full(coef_id, order, c_type);
//                     ParamId pid {
//                         ParameterType::WILSON, 
//                         GroupMapper::str(GroupMapper::enum_elt(groupName), ScaleType::HADRONIC, basis), 
//                         coef_lha
//                     };
//                     dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, coef_val, 0., (int)c_type));
//                 }
//             }
//         }
//     };
//     ports_config.iblock_c->compose_block(GroupMapper::str(GroupMapper::enum_elt(groupName), ScaleType::HADRONIC, basis), src, func);
// }

void CoefficientManager::init_group_hadronic_all_bases(const std::string &groupName, const std::string &order) {
    if (OrderMapper::enum_elt(order) < QCDOrder::NNLO) {
        fill_matching_groups(groupName, order);
    }

    for (auto basis : this->coefficientGroups.at(groupName)->get_bases()) {
        this->init_group_hadronic(groupName, order, basis);
    }
}

complex_t CoefficientManager::getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type) {
    // if (!this->coefficientGroups.contains(groupName)) {
    //     throw_no_group_error(groupName);
    // }

    complex_t c = this->coefficientGroups.at(groupName)->get_matching_coefficient(coeffName, order, cont_type);

    return c;
}

complex_t CoefficientManager::getFullMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type) {
    double fact = (*ports_config.wilson_proxy)("WPARAM_MATCH_SM", 1) / (4 * PI);
    size_t max_order = static_cast<size_t>(OrderMapper::enum_elt(order));
    complex_t c {0};
    for (size_t o = 1; o <= max_order; o++) {
        c += this->getMatchingCoefficient(groupName, coeffName, OrderMapper::str(static_cast<QCDOrder>(o)), cont_type) * std::pow(fact, o - 1);
    }
    return c;
}

complex_t CoefficientManager::getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type, WilsonBasis basis) {
    // if (!this->coefficientGroups.contains(groupName)) {
    //     throw_no_group_error(groupName);
    // }

    complex_t c = this->coefficientGroups.at(groupName)->get_running_coefficient(coeffName, order, cont_type, basis);
    return c;
}

complex_t CoefficientManager::getFullRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type, WilsonBasis basis) {
    double fact = (*ports_config.wilson_proxy)("WPARAM_RUN_SM", 1) / (4 * PI);
    size_t max_order = static_cast<size_t>(OrderMapper::enum_elt(order));
    complex_t c {0};
    for (size_t o = 1; o <= max_order; o++) {
        c += this->getRunCoefficient(groupName, coeffName, OrderMapper::str(static_cast<QCDOrder>(o)), cont_type, basis) * std::pow(fact, o - 1);
    }
    return c;
}

void CoefficientManager::registerCoefficientGroup(const std::string& groupName, std::shared_ptr<CoefficientGroup> group) {
    coefficientGroups[groupName] = group;
}

std::shared_ptr<CoefficientGroup> CoefficientManager::getCoefficientGroup(const std::string& groupName) const {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    return this->coefficientGroups.at(groupName);
}

std::map<std::string, std::shared_ptr<CoefficientGroup>> CoefficientManager::getGroups() {
    return this->coefficientGroups;
}

std::unordered_set<WilsonBasis> CoefficientManager::getGroupBases(WGroupId group) {
    return this->coefficientGroups.at(GroupMapper::str(group))->get_bases();
}

void CoefficientManager::printGroupCoefficients(const std::string& groupName) const {
    std::shared_ptr<CoefficientGroup> group = getCoefficientGroup(groupName);
}

CoefficientManager::~CoefficientManager() {
    LOG_TRACE("Call to CoefficientManager destructor");
    if (ports_config.iblock_c) {
        ports_config.iblock_c->remove_all_composed_blocks();
    }

}

void CoefficientManager::update(double mu_W, double mu_h) {
    this->set_matching_scale(mu_W);
    this->set_hadronic_scale(mu_h);
}

std::shared_ptr<CoefficientManager> CoefficientManager::Builder( std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double mu_W, double mu_h, std::string order, WilsonPortsConfig portconfig, std::map<Model, std::shared_ptr<IWilsonParameterHelper>> wilson_param_helpers) {
    
    for (auto& helper : wilson_param_helpers) {
        for (const auto& elem : groups) {
            helper.second->init(2, GroupMapper::enum_elt(elem.first));
        }
    }

    if (groups.empty()) {
        return std::make_shared<CoefficientManager>(portconfig);
    }

    auto manager = std::make_shared<CoefficientManager>(portconfig);
    for (auto& group : groups) {
        LOG_DEBUG("(CoefficientManager) Registering coefficient group", group.first);
        manager->registerCoefficientGroup(group.first, group.second);
    }
    LOG_DEBUG("(CoefficientManager) Setting matching scale");
    manager->set_matching_scale(mu_W);
    LOG_DEBUG("(CoefficientManager) Setting hadronic scale");
    manager->set_hadronic_scale(mu_h);
    for (auto& group: groups) {
        LOG_DEBUG("(CoefficientManager) Initializing group matching", group.first, "at", order);
        manager->init_group_matching(group.first, order);
        LOG_DEBUG("(CoefficientManager) Initializing group hadronic", group.first, "at", order); //TODO : Camilia change, need to be done correctly
        manager->init_group_hadronic_all_bases(group.first, order);
    }
    LOG_DEBUG("(CoefficientManager) Manager successfully initialized");
    return manager;
}

void CoefficientManager::set_hadronic_scale(double mu_h) {
    ports_config.scale_setter_api->switch_param(ScaleType::HADRONIC);
    ports_config.scale_setter_api->set(mu_h);
}

void CoefficientManager::set_matching_scale(double mu_W) {
    ports_config.scale_setter_api->switch_param(ScaleType::MATCHING);
    ports_config.scale_setter_api->set(mu_W);
}
