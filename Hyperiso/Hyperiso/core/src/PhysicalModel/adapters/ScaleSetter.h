#ifndef __SCALESETTER_H__
#define __SCALESETTER_H__

#include "IParamSetter.h"
#include "Include.h"

class ScaleSetter : public IParamSetter {
public:
    ScaleSetter(ScaleType scale_type) : scale_type(scale_type) {}

    void set(double value) override;

private:
    ScaleType scale_type;
};

#endif // __SCALESETTER_H__
