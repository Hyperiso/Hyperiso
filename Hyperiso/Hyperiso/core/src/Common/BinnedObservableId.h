#ifndef BINNEDOBSERVABLEID_H
#define BINNEDOBSERVABLEID_H

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "observable_ids.hpp"

/**
 * @file BinnedObservableId.h
 * @brief Utilities for identifying observables together with numerical bins.
 *
 * This header defines the @ref BinnedObservableId type and the helper
 * functions used to serialize bin boundaries into an integer-only FLHA/LhaID
 * representation. The encoding is designed to avoid storing raw floating-point
 * values directly inside an LhaID, while still allowing robust reconstruction
 * of the original bin edges up to @ref kMaxBinFracDigits decimal places.
 */

/**
 * @brief Normalizes signed zero to positive zero.
 *
 * Floating-point values can distinguish +0.0 and -0.0 at the bit level. For
 * observable bin boundaries, both values should be treated identically. This
 * helper maps both representations to +0.0 and leaves all other values
 * unchanged.
 *
 * @param x Floating-point value to normalize.
 * @return +0.0 if @p x is either +0.0 or -0.0, otherwise @p x unchanged.
 */
static inline double norm_zero(double x) noexcept {
    return (x == 0.0) ? 0.0 : x;
}

/**
 * @brief Returns the bit representation of a double after zero normalization.
 *
 * This helper is used when ordering and hashing bin boundaries. By normalizing
 * zero before bit-casting, +0.0 and -0.0 receive the same representation.
 *
 * @param x Floating-point value to convert.
 * @return Unsigned integer containing the IEEE-754 bit representation of @p x
 *         after zero normalization.
 */
static inline std::uint64_t bits_norm_zero(double x) noexcept {
    x = norm_zero(x);
    return std::bit_cast<std::uint64_t>(x);
}

/**
 * @struct EncodedBin
 * @brief Integer representation of a decimal bin boundary.
 *
 * A bin edge is encoded as three integer components:
 *   - an integer part,
 *   - a fractional part,
 *   - the number of decimal digits used for the fractional part.
 *
 * For example, the value 12.345 can be represented as:
 * @code
 * EncodedBin{12, 345, 3}
 * @endcode
 *
 * Negative values store the sign either in @ref int_part when the absolute
 * value is greater than or equal to 1, or in @ref frac_part for values between
 * -1 and 0.
 */
struct EncodedBin {
    long int_part;     /**< Integer part of the encoded value. */
    long frac_part;    /**< Fractional digits encoded as an integer. */
    long frac_digits;  /**< Number of decimal digits stored in @ref frac_part. */
};

/**
 * @brief Position at which bin components are inserted in a binned FLHA id.
 *
 * The bin encoding is inserted into the unbinned LhaID parts vector starting at
 * this index.
 */
static constexpr std::size_t kBinnedFlhaInsertPos = 2;

/**
 * @brief Number of integer components used to encode one complete bin.
 *
 * A bin is composed of two edges. Each edge uses three integer components:
 * integer part, fractional part, and number of fractional digits.
 */
static constexpr std::size_t kBinnedFlhaPartCount = 6;

/**
 * @brief Maximum number of decimal digits retained for a bin edge.
 *
 * Bin values are rounded to this number of fractional decimal digits before
 * being encoded into integer components.
 */
static constexpr long kMaxBinFracDigits = 12;

/**
 * @brief Computes 10 raised to a non-negative integer power.
 *
 * This small integer helper is used when reconstructing decimal values from
 * encoded fractional components.
 *
 * @param n Exponent. Expected to be non-negative.
 * @return 10^@p n as a long integer.
 */
static inline long pow10_long(long n) noexcept {
    long v = 1;
    while (n-- > 0) {
        v *= 10;
    }
    return v;
}

/**
 * @brief Removes insignificant trailing zeroes from a decimal string.
 *
 * If the string contains a decimal point, trailing zeroes are removed from the
 * fractional part. If this leaves the decimal point at the end, the decimal
 * point is removed as well.
 *
 * Examples:
 * @code
 * trim_decimal_string("12.3400") == "12.34"
 * trim_decimal_string("12.0000") == "12"
 * trim_decimal_string("0.0000")  == "0"
 * @endcode
 *
 * @param s Decimal string to trim.
 * @return Trimmed decimal string.
 */
static inline std::string trim_decimal_string(std::string s) {
    auto pos = s.find('.');
    if (pos == std::string::npos) {
        return s;
    }

    while (!s.empty() && s.back() == '0') {
        s.pop_back();
    }

    if (!s.empty() && s.back() == '.') {
        s.pop_back();
    }

    if (s.empty()) {
        return "0";
    }

    return s;
}

/**
 * @brief Encodes a bin boundary into integer FLHA-compatible components.
 *
 * The input value is first required to be finite, then normalized so that
 * signed zeroes are treated identically. The value is rounded to
 * @ref kMaxBinFracDigits decimal places and converted into an @ref EncodedBin.
 *
 * The sign convention is chosen so that values whose absolute value is at
 * least 1 carry the sign on @ref EncodedBin::int_part, while negative values
 * between -1 and 0 carry the sign on @ref EncodedBin::frac_part.
 *
 * @param x Bin boundary value, expressed in GeV.
 * @return Integer encoding of @p x.
 *
 * @throws std::runtime_error if @p x is not finite.
 * @throws std::invalid_argument or std::out_of_range if an intermediate string
 *         component cannot be converted by std::stol.
 */
static inline EncodedBin encode_bin_gev(double x) {
    if (!std::isfinite(x)) {
        throw std::runtime_error("Bin value must be finite");
    }

    x = norm_zero(x);

    constexpr long double kScale = 1.0e12L;

    long double xr = static_cast<long double>(x);
    xr = std::round(xr * kScale) / kScale;

    const bool neg = std::signbit(static_cast<double>(xr));
    const long double ax = std::fabs(xr);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(kMaxBinFracDigits) << ax;

    std::string s = trim_decimal_string(oss.str());

    const auto pos = s.find('.');
    std::string int_str  = (pos == std::string::npos) ? s : s.substr(0, pos);
    std::string frac_str = (pos == std::string::npos) ? "" : s.substr(pos + 1);

    if (int_str.empty()) {
        int_str = "0";
    }

    const long int_part = std::stol(int_str);
    const long frac_part = frac_str.empty() ? 0L : std::stol(frac_str);
    const long frac_digits = static_cast<long>(frac_str.size());

    if (neg) {
        if (int_part != 0) {
            return EncodedBin{-int_part, frac_part, frac_digits};
        }
        if (frac_part != 0) {
            return EncodedBin{0, -frac_part, frac_digits};
        }
    }

    return EncodedBin{int_part, frac_part, frac_digits};
}

/**
 * @brief Decodes integer FLHA-compatible bin components into a double value.
 *
 * This function is the inverse of @ref encode_bin_gev for valid encoded
 * values. The sign is reconstructed from either the integer part or the
 * fractional part.
 *
 * @param int_part Integer part of the encoded bin edge.
 * @param frac_part Fractional digits encoded as an integer.
 * @param frac_digits Number of decimal digits stored in @p frac_part.
 * @return Decoded bin boundary value, expressed in GeV.
 *
 * @throws std::runtime_error if @p frac_digits is outside the supported range.
 * @throws std::runtime_error if @p frac_part cannot fit inside
 *         @p frac_digits decimal digits.
 */
static inline double decode_bin_gev(long int_part, long frac_part, long frac_digits) {
    if (frac_digits < 0 || frac_digits > kMaxBinFracDigits) {
        throw std::runtime_error("decode_bin_gev: invalid frac_digits");
    }

    const long scale = pow10_long(frac_digits);
    if (std::labs(frac_part) >= scale && frac_digits > 0) {
        throw std::runtime_error("decode_bin_gev: frac_part does not fit frac_digits");
    }

    const bool neg = (int_part < 0) || (frac_part < 0);
    const long abs_int  = std::labs(int_part);
    const long abs_frac = std::labs(frac_part);

    long double out = static_cast<long double>(abs_int);

    if (frac_digits > 0) {
        out += static_cast<long double>(abs_frac) / static_cast<long double>(scale);
    }

    if (neg) {
        out = -out;
    }

    return norm_zero(static_cast<double>(out));
}

/**
 * @struct BinnedObservableId
 * @brief Identifies an observable together with a numerical bin.
 *
 * A BinnedObservableId combines an @ref ObservableId with a pair of bin
 * boundaries, typically interpreted as the lower and upper edges of the bin.
 * This is useful for observables whose experimental value is defined in a
 * specific kinematic interval.
 *
 * The bin is stored as a pair of doubles:
 *   - @ref p.first is the lower bin edge,
 *   - @ref p.second is the upper bin edge.
 *
 * The type provides equality comparison, strict weak ordering, hashing,
 * conversion to a robust FLHA/LhaID representation, reconstruction from that
 * representation, and a human-readable string representation.
 */
struct BinnedObservableId {
    ObservableId s;              /**< Identifier of the underlying observable. */
    std::pair<double, double> p; /**< Bin boundaries, usually {lower, upper}. */

    /**
     * @brief Default constructor.
     *
     * Leaves the observable identifier and bin boundaries default-initialized.
     */
    BinnedObservableId() = default;

    /**
     * @brief Constructs an unbinned observable identifier.
     *
     * The observable is initialized from @p id and the bin boundaries are set
     * to {0.0, 0.0}, which is used here as the default unbinned interval.
     *
     * @param id Identifier of the observable.
     */
    BinnedObservableId(ObservableId id) : s(id), p({0.0, 0.0}) {}

    /**
     * @brief Constructs a binned observable identifier.
     *
     * @param id Identifier of the observable.
     * @param bin Pair containing the lower and upper bin boundaries.
     */
    BinnedObservableId(ObservableId id, std::pair<double, double> bin) : s(id), p(bin) {}

    /**
     * @brief Constructs an unbinned observable identifier.
     *
     * The observable is initialized from @p id and the bin boundaries are set
     * to {0.0, 0.0}, which is used here as the default unbinned interval.
     *
     * @param id Identifier of the observable.
     */
    BinnedObservableId(Observables id) : s(ObservableMapper::to_id(id)), p({0.0, 0.0}) {}

    /**
     * @brief Constructs a binned observable identifier.
     *
     * @param id Identifier of the observable.
     * @param bin Pair containing the lower and upper bin boundaries.
     */
    BinnedObservableId(Observables id, std::pair<double, double> bin) : s(ObservableMapper::to_id(id)), p(bin) {}

    /**
     * @brief Equality comparison operator.
     *
     * Two BinnedObservableId objects are equal if they refer to the same
     * observable and have identical normalized bin boundaries. Positive and
     * negative zero are treated as the same value.
     *
     * @param o Object to compare with.
     * @return true if both identifiers represent the same binned observable,
     *         false otherwise.
     */
    bool operator==(BinnedObservableId const& o) const noexcept {
        return s == o.s
            && norm_zero(p.first)  == norm_zero(o.p.first)
            && norm_zero(p.second) == norm_zero(o.p.second);
    }

    /**
     * @brief Strict weak ordering for binned observable identifiers.
     *
     * Ordering is lexicographic on the observable identifier and the bitwise
     * representation of the normalized bin edges. This allows the type to be
     * used as a key in ordered containers such as std::map or std::set.
     *
     * @param o Object to compare with.
     * @return true if this object is ordered before @p o.
     */
    bool operator<(BinnedObservableId const& o) const noexcept {
        const std::uint64_t a = bits_norm_zero(p.first);
        const std::uint64_t b = bits_norm_zero(p.second);
        const std::uint64_t c = bits_norm_zero(o.p.first);
        const std::uint64_t d = bits_norm_zero(o.p.second);

        return std::tie(s, a, b) < std::tie(o.s, c, d);
    }

    /**
     * @brief Converts this binned observable identifier to its FLHA encoding.
     *
     * The underlying observable is first converted to its unbinned FLHA
     * identifier using @ref ObservableMapper::flha_of. The lower and upper bin
     * boundaries are then encoded as six integer components and inserted into
     * the resulting LhaID at @ref kBinnedFlhaInsertPos.
     *
     * @return LhaID containing both the observable identifier and the encoded
     *         bin boundaries.
     *
     * @throws std::runtime_error if no FLHA mapping is known for the
     *         underlying observable.
     * @throws std::runtime_error if one of the bin boundaries cannot be
     *         encoded, for example because it is not finite.
     */
    LhaID flha() const {
        auto flha_opt = ObservableMapper::flha_of(this->s);
        if (!flha_opt.has_value()) {
            throw std::runtime_error("ObservableId to flha mapping unknown");
        }

        LhaID unbinned_id = flha_opt.value();

        const EncodedBin low  = encode_bin_gev(this->p.first);
        const EncodedBin high = encode_bin_gev(this->p.second);

        std::vector<long> bin_parts{
            low.int_part,
            low.frac_part,
            low.frac_digits,
            high.int_part,
            high.frac_part,
            high.frac_digits
        };

        auto parts = unbinned_id.get_parts();
        parts.insert(parts.begin() + kBinnedFlhaInsertPos, bin_parts.begin(), bin_parts.end());

        return LhaID(parts);
    }

    /**
     * @brief Reconstructs a binned observable identifier from an FLHA encoding.
     *
     * This function expects an LhaID containing the encoded bin information at
     * @ref kBinnedFlhaInsertPos. The six bin-related components are extracted,
     * decoded into lower and upper bin boundaries, and removed before
     * converting the remaining unbinned identifier back to an @ref ObservableId.
     *
     * @param id FLHA/LhaID representation of a binned observable.
     * @return Reconstructed BinnedObservableId.
     *
     * @throws std::runtime_error if @p id does not contain enough components
     *         for the encoded bin.
     * @throws std::runtime_error if the unbinned FLHA identifier cannot be
     *         mapped back to an ObservableId.
     * @throws std::runtime_error if the encoded bin components are invalid.
     */
    static BinnedObservableId from_flha(LhaID const& id) {
        auto parts = id.get_parts();

        if (parts.size() < kBinnedFlhaInsertPos + kBinnedFlhaPartCount) {
            throw std::runtime_error("from_flha: LhaID has not enough parts to contain robust binning");
        }

        const long low_int      = parts.at(kBinnedFlhaInsertPos + 0);
        const long low_frac     = parts.at(kBinnedFlhaInsertPos + 1);
        const long low_ndigits  = parts.at(kBinnedFlhaInsertPos + 2);
        const long high_int     = parts.at(kBinnedFlhaInsertPos + 3);
        const long high_frac    = parts.at(kBinnedFlhaInsertPos + 4);
        const long high_ndigits = parts.at(kBinnedFlhaInsertPos + 5);

        parts.erase(parts.begin() + kBinnedFlhaInsertPos,
                    parts.begin() + kBinnedFlhaInsertPos + kBinnedFlhaPartCount);
        LhaID unbinned_id(parts);

        auto obs_opt = ObservableMapper::from_flha(unbinned_id);
        if (!obs_opt.has_value()) {
            throw std::runtime_error("from_flha: flha to ObservableId mapping unknown");
        }

        BinnedObservableId out;
        out.s = obs_opt.value();
        out.p = {
            decode_bin_gev(low_int, low_frac, low_ndigits),
            decode_bin_gev(high_int, high_frac, high_ndigits)
        };

        return out;
    }

    /**
     * @brief Returns a human-readable representation of the binned observable.
     *
     * The returned string contains the observable string representation
     * followed by the bin boundaries.
     *
     * @return String of the form "observable [low, high]".
     */
    std::string str() const {
        std::stringstream ss;
        ss << s.str() << " [" << p.first << ", " << p.second << "]";
        return ss.str();
    }
};

/**
 * @brief Hash functor specialization for BinnedObservableId.
 *
 * This specialization allows BinnedObservableId to be used as a key in
 * std::unordered_map or std::unordered_set. The hash combines the hash of the
 * underlying observable identifier with the bitwise representation of the
 * normalized bin boundaries.
 */
template<>
struct std::hash<BinnedObservableId> {
    /**
     * @brief Computes the hash of a binned observable identifier.
     *
     * Positive and negative zero are normalized before hashing the bin
     * boundaries, consistently with equality comparison.
     *
     * @param id Binned observable identifier to hash.
     * @return Hash value for @p id.
     */
    std::size_t operator()(BinnedObservableId const& id) const noexcept {
        std::size_t h = std::hash<ObservableId>{}(id.s);

        auto mix = [](std::size_t& seed, std::size_t v) {
            seed ^= v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        };

        mix(h, std::hash<std::uint64_t>{}(bits_norm_zero(id.p.first)));
        mix(h, std::hash<std::uint64_t>{}(bits_norm_zero(id.p.second)));

        return h;
    }
};

#endif // BINNEDOBSERVABLEID_H
