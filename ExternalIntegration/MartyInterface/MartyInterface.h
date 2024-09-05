#ifndef MARTY_INTERFACE_H
#define MARTY_INTERFACE_H

#include "CodeGenerator.h"
#include "GppCompilerStrategy.h"
#include "SMModelModifier.h"
// #include "THDMModelModifier.h"

class MartyInterface {
public:
    void generate();
    void run();
};

#endif // MARTY_INTERFACE_H
