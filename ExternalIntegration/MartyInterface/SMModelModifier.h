#ifndef SM_MODEL_MODIFIER_H
#define SM_MODEL_MODIFIER_H

#include <unordered_map>
#include "ModelModifier.h"

class SMModelModifier : public ModelModifier {
public:
    SMModelModifier(std::string wilson) {this->wilson = wilson;}
    void modifyLine(std::string& line) override;
};

#endif // SM_MODEL_MODIFIER_H
