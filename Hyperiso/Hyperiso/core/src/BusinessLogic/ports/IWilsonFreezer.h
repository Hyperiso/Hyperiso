#ifndef __IWILSONFREEZER_H__
#define __IWILSONFREEZER_H__

#include "Freezer.h"

template<typename T>
class IWilsonFreezer {
public:
    virtual void freeze(T) = 0;
    virtual void unfreeze(T) = 0;
};

#endif // __IWILSONFREEZER_H__
