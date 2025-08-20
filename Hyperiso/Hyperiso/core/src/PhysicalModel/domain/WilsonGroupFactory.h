#ifndef __WILSONGROUPFACTORY_H__
#define __WILSONGROUPFACTORY_H__

#include "Include.h"
#include "WilsonGroup.h"
#include "BWilsonGroup.h"
#include "MesonMixingWilsonGroup.h"
#include "Wilson_SUSY_super.h"
#include "Wilson_THDM_super.h"

class WilsonGroupFactory {
public:
    static std::shared_ptr<CoefficientGroup> create_coefficient_group(WGroup group, Model model) {
        if (UseMarty().get()) {
            model = Model::SM;
        }

        switch (model) {
            case Model::SM:
            switch (group) {
                case WGroup::B:
                    return std::make_shared<BCoefficientGroup>();
                case WGroup::BPrime:
                    return std::make_shared<BPrimeCoefficientGroup>();
                case WGroup::BScalar:
                    return std::make_shared<BScalarCoefficientGroup>();
                case WGroup::BCC:
                    return std::make_shared<BclnuCoefficientGroup>();
                case WGroup::MESON_MIXING:
                    return std::make_shared<MesonMixingCoefficientGroup>();
                default:
                    LOG_ERROR("Invalid Argument", "Unknown group type", GroupMapper::str(group));
            }
            case Model::SUSY:
            switch (group) {
                case WGroup::B:
                    return std::make_shared<BCoefficientGroup_susy>();
                case WGroup::BPrime:
                    return std::make_shared<BPrimeCoefficientGroup_susy>();
                case WGroup::BScalar:
                    return std::make_shared<BScalarCoefficientGroup_susy>();
                case WGroup::BCC:
                    return std::make_shared<BclnuCoefficientGroup_SUSY>();
                case WGroup::MESON_MIXING:
                    LOG_ERROR("Not Implemented", "Coefficient group", GroupMapper::str(group), "has not yet been implmented in SUSY");
                default:
                    LOG_ERROR("Invalid Argument", "Unknown group type", GroupMapper::str(group));
            }
            case Model::THDM:
            switch (group) {
                case WGroup::B:
                    return std::make_shared<BCoefficientGroup_THDM>();
                case WGroup::BPrime:
                    return std::make_shared<BPrimeCoefficientGroup_THDM>();
                case WGroup::BScalar:
                    return std::make_shared<BScalarCoefficientGroup_THDM>();
                case WGroup::BCC:
                    return std::make_shared<BclnuCoefficientGroup_THDM>();
                case WGroup::MESON_MIXING:
                    LOG_ERROR("Not Implemented", "Coefficient group", GroupMapper::str(group), "has not yet been implmented in THDM");
                default:
                    LOG_ERROR("Invalid Argument", "Unknown group type", GroupMapper::str(group));
            }
            default:
            LOG_ERROR("Invalid Argument", "Unknown model", ModelMapper::str(model));
        }
    }
};

#endif // __WILSONGROUPFACTORY_H__
