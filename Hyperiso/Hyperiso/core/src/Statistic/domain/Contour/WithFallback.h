// #ifndef __WITHFALLBACK_H__
// #define __WITHFALLBACK_H__

// #include "IContourExtractor.h"
// #include "Logger.h"

// class WithFallback : public IContourExtractor {
// public:
//     WithFallback(std::shared_ptr<IContourExtractor> primary, std::shared_ptr<IContourExtractor> fallback) : primary(std::move(primary)), fallback(std::move(fallback)) {}

//     Contour extract(const ScalarField2D& field, const ContourRequest& cr) override {
//         auto res = primary->extract(field, cr);
//         if (!res.success) {
//             LOG_WARN("Primary contour extraction method failed, falling back.");
//             res = fallback->extract(field, cr);
//         }
//         return res;
//     }

// private:
//     std::shared_ptr<IContourExtractor> primary;
//     std::shared_ptr<IContourExtractor> fallback;
// };

// #endif // __WITHFALLBACK_H__
//for segfault

#ifndef WITHFALLBACK_H
#define WITHFALLBACK_H

#include <exception>

#include "IContourExtractor.h"
#include "Logger.h"

/**
 * @file WithFallback.h
 * @brief Decorator adding fallback behavior to contour extraction.
 *
 * @see IContourExtractor
 */

/**
 * @class WithFallback
 * @brief Tries a primary contour extractor and falls back to another extractor on failure.
 *
 * The wrapper is intended to make contour production more robust when a fast or
 * preferred extraction method can occasionally fail. The fallback is used when
 * the primary extractor throws an exception or returns a contour marked as
 * unsuccessful.
 */
class WithFallback : public IContourExtractor {
public:
    /**
     * @brief Constructs a fallback contour extractor wrapper.
     *
     * @param primary Extractor attempted first.
     * @param fallback Extractor used if the primary fails or returns an
     *                 unsuccessful contour.
     */
    WithFallback(std::shared_ptr<IContourExtractor> primary, std::shared_ptr<IContourExtractor> fallback)
        : primary(std::move(primary)), fallback(std::move(fallback)) {}

    /**
     * @brief Extracts a contour using the primary extractor, then fallback if needed.
     *
     * @param field Scalar field from which to extract the contour.
     * @param cr Contour extraction request.
     * @return Contour produced by the primary extractor when successful,
     *         otherwise by the fallback extractor.
     *
     * @throws std::exception propagated from the fallback extractor if fallback
     *         extraction also fails.
     */
    Contour extract(const ScalarField2D& field, const ContourRequest& cr) override {
        Contour res;
        bool primary_returned = false;

        try {
            res = primary->extract(field, cr);
            primary_returned = true;
        } catch (const std::exception& e) {
            LOG_WARN("Primary contour extraction threw, falling back: ", e.what());
        } catch (...) {
            LOG_WARN("Primary contour extraction threw unknown exception, falling back.");
        }

        if (!primary_returned || !res.success) {
            LOG_WARN("Primary contour extraction method failed or returned no contour, falling back.");
            res = fallback->extract(field, cr);
        }

        return res;
    }

private:
    std::shared_ptr<IContourExtractor> primary;     ///< Preferred extractor attempted first.
    std::shared_ptr<IContourExtractor> fallback;    ///< Secondary extractor used when the primary fails.
};

#endif // WITHFALLBACK_H
