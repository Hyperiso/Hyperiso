#ifndef __IQCDPROVIDER_H__
#define __IQCDPROVIDER_H__

#include "AbstractConfig.h"
#include "QCDHelper.h"

class IQCDProvider {
public:
    virtual ~IQCDProvider() = default;

    virtual QCDConstants* get_constants() = 0;
};

#endif // __IQCDPROVIDER_H__
