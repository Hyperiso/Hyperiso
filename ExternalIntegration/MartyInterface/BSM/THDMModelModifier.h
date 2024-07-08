#ifndef THDM_MODEL_MODIFIER_H
#define THDM_MODEL_MODIFIER_H

#include "ModelModifier.h"

class THDMModelModifier : public ModelModifier {
public:
    void modifyLine(std::string& line) override;
};

#endif // THDM_MODEL_MODIFIER_H
