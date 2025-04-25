#ifndef __IFREEZER_H__
#define __IFREEZER_H__

template<typename T, typename U, typename V>
class IFreezer {
public:
    virtual ~IFreezer() = default;

    /**
     * @brief Freezes the object, preventing further updates.
     */
    static void freeze(const T&, const U&) = delete;

    /**
     * @brief Freezes the object, preventing further updates.
     */
    static void freeze(const V&) = delete;

    /**
     * @brief Freezes the object, preventing further updates.
     */
    static void unfreeze(const T&, const U&) = delete;

    /**
     * @brief Freezes the object, preventing further updates.
     */
    static void unfreeze(const V&) = delete;
};


#endif // __IFREEZER_H__
