#ifndef QCDADAPTER_H
#define QCDADAPTER_H

#include "IDataProvider.h"
#include "IQCDProvider.h"
#include "QCDHelper.h"
#include "Configs.h"

/**
 * @class QCDProvider
 * @ingroup DataProvidersModule
 * @brief Provides strong coupling constant and MSbar mass values at various energy scales.
 */
class QCDProvider : public IDataProvider<QCDProvider>, public IQCDProvider {
public:

    /**
     * @brief Computes the strong coupling constant \f$\alpha_s\f$ at a given scale.
     * @param config AlphasConfig object containing scale and mass type settings.
     * @return Value of \f$\alpha_s\f$ at the specified scale.
     */
    double operator()(AlphasConfig);

    /**
     * @brief Computes the MS-bar running mass at a given scale for a quark.
     * @param config MassConfig object containing PDG ID and scale information.
     * @return The MS-bar mass at the specified scale.
     */
    double operator()(MassConfig);

    /**
     * @brief Retrieves the constants related to QCD calculations.
     * @return Pointer to QCDConstants structure.
     */
    QCDConstants* get_constants() override;
};


#endif // QCDADAPTER_H
