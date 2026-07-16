#ifndef MNCONTOUREXTRACTOR_H
#define MNCONTOUREXTRACTOR_H

#include "IContourExtractor.h"

/**
 * @file MnContourExtractor.h
 * @brief Minuit-backed implementation of the contour extraction interface.
 *
 * @see IContourExtractor
 */

/**
 * @class MnContourExtractor
 * @brief Extracts a 2D contour using the configured Minuit fit backend.
 *
 * The extractor first minimizes the scalar field in the two displayed
 * coordinates, then delegates contour construction to the backend contour
 * routine. The target level is passed as the Minuit contour `up` value.
 */
class MnContourExtractor: public IContourExtractor {
public:
    /**
     * @brief Computes a Minuit contour for the requested scalar field level.
     *
     * @param field Scalar field evaluated as \f$f(x,y)\f$.
     * @param cr Contour extraction request containing target level and axis
     *           parameter definitions.
     * @return Contour returned by the Minuit backend.
     */
    Contour extract(const ScalarField2D& field, const ContourRequest& cr) override;
};

#endif // MNCONTOUREXTRACTOR_H
