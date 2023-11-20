#ifndef ISAJET_INTERFACE_H
#define ISAJET_INTERFACE_H

#include <string>
#include "ModelParameters.h"
#include "ModelOutput.h"

class IsajetInterface {
public:
    IsajetInterface(const std::string& executablePath);
    ModelOutput runIsajet(const ModelParameters& params);

private:
    std::string isajetExecutable;
    void writeConfigFile(const ModelParameters& params);
    ModelOutput parseOutputFile(const std::string& filename);
};

#endif // ISAJET_INTERFACE_H
