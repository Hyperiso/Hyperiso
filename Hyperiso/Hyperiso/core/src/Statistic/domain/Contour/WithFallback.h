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
#ifndef __WITHFALLBACK_H__
#define __WITHFALLBACK_H__

#include "IContourExtractor.h"
#include "Logger.h"

#include <exception>

class WithFallback : public IContourExtractor {
public:
    WithFallback(std::shared_ptr<IContourExtractor> primary, std::shared_ptr<IContourExtractor> fallback)
        : primary(std::move(primary)), fallback(std::move(fallback)) {}

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
    std::shared_ptr<IContourExtractor> primary;
    std::shared_ptr<IContourExtractor> fallback;
};

#endif // __WITHFALLBACK_H__
