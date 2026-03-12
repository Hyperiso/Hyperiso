#ifndef __WITHFALLBACK_H__
#define __WITHFALLBACK_H__

#include "IContourExtractor.h"
#include "Logger.h"

class WithFallback : public IContourExtractor {
public:
    WithFallback(std::shared_ptr<IContourExtractor> primary, std::shared_ptr<IContourExtractor> fallback) : primary(std::move(primary)), fallback(std::move(fallback)) {}

    void extract(std::shared_ptr<Contour> contour) override {
        primary->extract(contour);
        if (!contour->success) {
            LOG_WARN("Primary contour extraction method failed, falling back.");
            fallback->extract(contour);
        }
    }

private:
    std::shared_ptr<IContourExtractor> primary;
    std::shared_ptr<IContourExtractor> fallback;
};

#endif // __WITHFALLBACK_H__
