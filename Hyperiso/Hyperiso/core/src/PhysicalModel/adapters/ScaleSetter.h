#ifndef __SCALESETTER_H__
#define __SCALESETTER_H__

#include "IParamSetter.h"
#include "Include.h"

class ScaleSetter : public IParamSetter<ScaleType> {
public:
    ScaleSetter(ScaleType scale_type) {this->param = scale_type;}

    void switch_param(ScaleType scale_type) override {this->param = scale_type;}

    void set(double value) override;

};

#endif // __SCALESETTER_H__
