#ifndef __SCALESETTER_H__
#define __SCALESETTER_H__

#include "IParamSetter.h"

class MatchingScaleSetter : public IParamSetter {
public:
    void set(double value) override;
};

class HadronicScaleSetter : public IParamSetter {
public:
    void set(double value) override;
};

#endif // __SCALESETTER_H__
