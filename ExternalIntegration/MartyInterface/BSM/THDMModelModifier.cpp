#include "THDMModelModifier.h"

void THDMModelModifier::modifyLine(std::string& line) {
    if (line.find("//MODEL_SPECIFIC_CODE") != std::string::npos) {
        line.replace(line.find("//MODEL_SPECIFIC_CODE"), 20, "THDM specific code");
    }
}
