#ifndef WILSONFREEZER_H
#define WILSONFREEZER_H

#include "IWilsonFreezer.h"
#include "Freezer.h"
#include "Include.h"
#include "ObsWilsonProxy.h"
#include "ObsWilsonBuilder.h"

/**
 * @file WilsonFreezer.h
 * @brief Concrete freezer for Wilson coefficient blocks (matching + hadronic).
 *
 * This implementation freezes/unfreezes:
 *  - the matching-scale block of a Wilson group,
 *  - all hadronic-scale blocks for every available @ref WilsonBasis of that group.
 *
 * The set of bases is queried from the observable Wilson proxy.
 *
 * Typical usage:
 * @code
 *   auto builder = std::make_shared<ObsWilsonBuilder>(...);
 *   WilsonFreezer freezer(builder);
 *   freezer.freeze(GroupMapper::to_id(WGroup::B));
 *   ...
 *   freezer.unfreeze(GroupMapper::to_id(WGroup::B));
 * @endcode
 *
 * @see IWilsonFreezer
 * @see Freezer
 * @see ObsWilsonProxy
 * @see GroupMapper
 */
class WilsonFreezer : public IWilsonFreezer<WGroupId> {
public:
    /**
     * @brief Construct a freezer using an observable Wilson builder.
     * @param wil_builder Builder used to obtain the Wilson proxy and discover bases.
     */
    WilsonFreezer(const std::shared_ptr<IObsWilsonBuilder> &wil_builder) {
        this->w_proxy = wil_builder->get_proxy();
    }

    /**
     * @brief Freezes matching and hadronic blocks for the given group.
     */
    void freeze(WGroupId group) override;

    /**
     * @brief Unfreezes matching and hadronic blocks for the given group.
     */
    void unfreeze(WGroupId group) override;

private:
    std::shared_ptr<IObsWilsonProxy> w_proxy;
};

#endif // WILSONFREEZER_H
