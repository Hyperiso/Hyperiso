#include "ICoreAPI.h"
#include "HyperisoMaster.h"

/**
 * @file SMFromHypProxy.h
 * @brief Core API exposing whether Hyperiso should provide the SM part when using Marty.
 *
 * This header defines @ref SMFromHypProxy, a concrete @ref ICoreAPI<bool> implementation.
 * It queries Hyperiso's active configuration flag @ref ExternalFlag::HYP_AS_SM_MARTY
 * through @ref HyperisoMaster.
 *
 * The flag is intended to control workflows where Marty is used for Wilson coefficient
 * generation but the SM contribution (or SM defaults) should be taken from Hyperiso's
 * built-in/hard-coded SM implementation instead of relying on Marty for SM beyond LO.
 *
 * @see HyperisoMaster
 * @see HyperisoConfig
 * @see ExternalFlag
 */

/**
 * @class SMFromHypProxy
 * @ingroup MonitoringSystemModule
 * @brief Returns whether Hyperiso is configured to provide the SM part in the Marty workflow.
 *
 * This class is a lightweight adapter over @ref HyperisoMaster:
 * @code
 *   SMFromHypProxy smFromHyp;
 *   if (smFromHyp.get()) {
 *       // Use Hyperiso/builtin SM part (e.g. to keep SM up to requested QCD order)
 *   }
 * @endcode
 *
 * Internally, it delegates to:
 * @code
 *   HyperisoMaster().check_flag(ExternalFlag::HYP_AS_SM_MARTY);
 * @endcode
 *
 * @note This reads the flag from the active HyperisoConfig (stored in the framework cache).
 */
class SMFromHypProxy : public ICoreAPI<bool> {
public:
    /**
     * @brief Default constructor.
     */
    explicit SMFromHypProxy() = default;

    /**
     * @brief Returns true if @ref ExternalFlag::HYP_AS_SM_MARTY is active.
     *
     * @return true if Hyperiso should be used as SM provider in the Marty workflow, false otherwise.
     */
    inline bool get() override { return HyperisoMaster().check_flag(ExternalFlag::HYP_AS_SM_MARTY);}

};