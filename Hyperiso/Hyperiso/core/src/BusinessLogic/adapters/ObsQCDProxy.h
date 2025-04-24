#ifndef OBS_QCD_PROXY_H
#define OBS_QCD_PROXY_H

#include "QCDProvider.h"

class ObsQCDProxy {
public:
    double operator()(AlphasConfig config);
    double operator()(MassConfig config);

    QCDConstants *get_constants();

};

#endif