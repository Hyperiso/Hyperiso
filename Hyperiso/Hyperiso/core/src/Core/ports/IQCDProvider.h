#ifndef IQCDPROVIDER_H
#define IQCDPROVIDER_H

#include "AbstractConfig.h"
#include "QCDHelper.h"

/**
 * @class IQCDProvider
 * @ingroup DataProvidersModule
 * @brief Abstract interface for accessing QCD constants.
 */
class IQCDProvider {
public:
    virtual ~IQCDProvider() = default;

    /**
     * @brief Virtual method to retrieve QCD constants.
     * @return Pointer to a QCDConstants structure.
     */
    virtual QCDConstants* get_constants() = 0;
};

#endif // IQCDPROVIDER_H
