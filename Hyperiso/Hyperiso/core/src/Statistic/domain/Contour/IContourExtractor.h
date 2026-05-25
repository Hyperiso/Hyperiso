#ifndef ICONTOUREXTRACTOR_H
#define ICONTOUREXTRACTOR_H

#include <vector>
#include <set>
#include <memory>
#include <functional>

#include "Math.h"

/**
 * @file IContourExtractor.h
 * @brief Common data structures and interface for 2D contour extraction.
 *
 * A contour extractor receives a scalar field \f$f(x,y)\f$ and returns the
 * points approximating the level set \f$f(x,y)=c\f$ requested by
 * @ref ContourRequest::level.
 */

/**
 * @struct ContourRequest
 * @brief Input configuration for a two-dimensional contour extraction.
 */
struct ContourRequest {
    double level;   ///< Target scalar-field level to extract.

    /**
     * @brief Extraction domain as {xmin, xmax, ymin, ymax}.
     */
    std::array<double, 4> bounds;

    /**
     * @brief Parameter definitions for the two scanned axes.
     *
     * The definitions provide central values, step hints, and bounds to
     * extractor implementations that rely on optimizer-backed contouring.
     */
    std::array<fit_app::ParameterDefinition, 2> p_defs;

    std::size_t resolution = 40;    ///< Requested grid or sampling resolution.
};

/**
 * @struct Contour
 * @brief Output of a contour extraction algorithm.
 */
struct Contour {
    std::set<Path> paths;   ///< Extracted contour paths.
    double level;           ///< Level actually targeted by the extraction.
    bool success = false;   ///< Whether extraction produced a valid contour.
};

/**
 * @class IContourExtractor
 * @brief Abstract interface for algorithms extracting 2D contour paths.
 *
 * Concrete implementations may use marching squares, Minuit contours, adaptive
 * sampling, or any other method. The interface deliberately exposes only the
 * scalar field and request metadata so that algorithms remain interchangeable.
 */
class IContourExtractor {
public:
    virtual ~IContourExtractor() = default;

    /**
     * @brief Extracts contour paths from a scalar field.
     *
     * @param field Scalar function evaluated as \f$f(x,y)\f$.
     * @param cr Contour extraction request.
     * @return Extracted contour object, including paths and success flag.
     *
     * @throws std::exception if the concrete extraction algorithm encounters a
     *         non-recoverable error.
     */
    virtual Contour extract(const ScalarField2D& field, const ContourRequest& cr) = 0;
};


#endif // ICONTOUREXTRACTOR_H
