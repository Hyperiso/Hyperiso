#include "CopulaFactory.h"

std::unique_ptr<ICopula> make(unsigned int seed, CopulaConfig cfg) {
    return std::visit([seed](auto&& c) -> std::unique_ptr<ICopula> {
        using T = std::decay_t<decltype(c)>;
        if constexpr (std::is_same_v<T, GaussianCopulaConfig>)
            return std::make_unique<GaussianCopula>(seed, c.R);
        else if constexpr (std::is_same_v<T, StudentTCopulaConfig>)
            return std::make_unique<StudentTCopula>(seed, c.R, c.nu);
    }, cfg);
}

std::unique_ptr<ICopula> CopulaFactory::create(CopulaType name,
                                               CopulaConfig config,
                                               unsigned int seed) {
    switch (name) {
        case CopulaType::GAUSSIAN:
            return make(seed, config);
        case CopulaType::STUDENT_T:
            return make(seed, config);
        default:
            throw std::invalid_argument("Unknown copula.");
    }
}