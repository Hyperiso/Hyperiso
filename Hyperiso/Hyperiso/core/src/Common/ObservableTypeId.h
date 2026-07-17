#ifndef OBSERVABLE_TYPE_ID_H
#define OBSERVABLE_TYPE_ID_H

#include <stdexcept>

/**
 * @file ObservableTypeId.h
 * @brief Helpers for project-defined FLHA observable-type identifiers.
 */

/**
 * @brief Polarization category encoded by the project-specific FLHA pattern
 *        @c 92ij.
 *
 * The decimal digit @c i identifies the polarization category:
 *
 * - 1: fermion polarization;
 * - 2: longitudinal vector polarization;
 * - 3: transverse vector polarization.
 */
enum class PolarizationKind : long {
    Fermion = 1,
    LongitudinalVector = 2,
    TransverseVector = 3,
};

/** @brief FLHA observable type for the transverse fraction @f$F_T@f$. */
inline constexpr long kTransverseFractionTypeId = 932;

/** @brief FLHA observable type for the angular coefficient @f$\alpha_K@f$. */
inline constexpr long kAlphaKTypeId = 933;

/** @brief Legacy v1.0.0--v1.0.1 type used for charged-lepton polarization. */
inline constexpr long kLegacyFermionPolarizationTypeId = 92015;

/** @brief Legacy v1.0.0--v1.0.1 type used for longitudinal vector polarization. */
inline constexpr long kLegacyLongitudinalVectorPolarizationTypeId = 921423;

/**
 * @brief Build a project-defined FLHA polarization type id using @c 92ij.
 *
 * The decimal digit @c j is the one-based position of the particle of interest
 * in the ordered list of daughter particles stored in the FLHA observable id.
 * For example, in @f$B\to D\tau\nu@f$, the tau is daughter 2, hence its
 * fermion-polarization type id is @c 9212.  In
 * @f$B\to D^*\tau\nu@f$, the @f$D^*@f$ is daughter 1, hence its longitudinal
 * vector-polarization type id is @c 9221.
 *
 * @param kind Polarization category encoded in digit @c i.
 * @param daughter_position One-based daughter position encoded in digit @c j.
 * @return The four-digit type id @c 92ij.
 * @throws std::invalid_argument if @p kind is outside 1--3 or if
 *         @p daughter_position is outside 1--9.
 */
constexpr long make_polarization_type_id(
    PolarizationKind kind,
    unsigned int daughter_position
) {
    const long kind_digit = static_cast<long>(kind);
    if (kind_digit < 1 || kind_digit > 3) {
        throw std::invalid_argument("polarization kind must be encoded by i=1,2,3");
    }
    if (daughter_position < 1 || daughter_position > 9) {
        throw std::invalid_argument("polarized daughter position must be in [1,9]");
    }
    return 9200L + 10L * kind_digit + static_cast<long>(daughter_position);
}

static_assert(make_polarization_type_id(PolarizationKind::Fermion, 2) == 9212);
static_assert(make_polarization_type_id(PolarizationKind::LongitudinalVector, 1) == 9221);
static_assert(make_polarization_type_id(PolarizationKind::TransverseVector, 1) == 9231);

#endif // OBSERVABLE_TYPE_ID_H
