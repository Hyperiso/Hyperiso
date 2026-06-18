#ifndef STATISTIC_CACHE_PRINTER_H
#define STATISTIC_CACHE_PRINTER_H

#include <iosfwd>

#include "StatisticManager.h"

/**
 * @class StatisticCachePrinter
 * @brief Formatting helper for @ref StatCache diagnostics.
 *
 * Keeping cache formatting outside StatisticManager makes the manager focus on
 * state construction and computation.  Printing remains opt-in through
 * StatisticConfig and can be expanded without growing the manager itself.
 */
class StatisticCachePrinter {
public:
    /** Print a compact diagnostic dump of the statistic cache. */
    static void print(const StatCache& cache, std::ostream& os);
};

#endif
