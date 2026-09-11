#ifndef HYPERISO_MARTY_NUMERICAL_POLICY_H
#define HYPERISO_MARTY_NUMERICAL_POLICY_H

/**
 * @brief Internal numerical prescriptions for MARTY-generated Wilson libraries.
 *
 * These are deliberately code-level constants rather than user configuration:
 * changing a regulator prescription is a developer/validation choice.  Keep
 * the ordinary matching regulator small; the O(1) value is reserved for the
 * opt-in raw massless-photon diagnostic used by C9/CP9 (and split CP10).
 */
namespace MartyNumericalPolicy {
inline constexpr double kDefaultRegProp = 1e-10;
inline constexpr double kPhotonDiagnosticRegProp = 1.0;
}

#endif
