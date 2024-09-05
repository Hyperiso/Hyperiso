#ifndef MARTY_INTERFACE_H
#define MARTY_INTERFACE_H

#include "CodeGenerator.h"
#include "GppCompilerStrategy.h"
#include "SMModelModifier.h"
// #include "THDMModelModifier.h"

class MartyInterface {
public:
    void generate(std::string wilson);
    void compile(std::string wilson);

private:
    bool generated = false;
    bool compiled = false;
};

#endif // MARTY_INTERFACE_H
