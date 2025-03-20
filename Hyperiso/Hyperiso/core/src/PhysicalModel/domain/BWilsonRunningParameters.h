#ifndef __BWILSONRUNNINGPARAMETERS_H__
#define __BWILSONRUNNINGPARAMETERS_H__

#include <array>

struct BWilsonRunningParameters {
    static constexpr int array_size {10};

    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m00;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m10;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m11;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m20;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m21;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> m22;

    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> l00;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> l01;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> l10;
    static const std::array<std::array<std::array<double, array_size>, array_size>, array_size> l11;

    static constexpr std::array<double, array_size> ai = {14.0 / 23.0, 16.0 / 23.0, 6.0 / 23.0, -12.0 / 23.0, 0.408619, -0.422989, -0.899395, 0.145649, -1.0, -1.0};
    static constexpr std::array<double, array_size> ai2 = {6./23., -12./23., 0.4086, -0.4230, -0.8994, 0.1456, 16./23., 14./23., 11./23., 29./23.}; 
};

#endif // __BWILSONRUNNINGPARAMETERS_H__
