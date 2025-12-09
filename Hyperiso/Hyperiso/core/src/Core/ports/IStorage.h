#ifndef ISTORAGE_H
#define ISTORAGE_H

#include <memory>

/**
 * @brief Generic storage interface for key–value containers.
 *
 * This template defines the minimal operations required for a storage-like
 * container:
 *  - insertion (@ref store),
 *  - replacement (@ref assign),
 *  - store-or-assign,
 *  - membership test,
 *  - retrieval,
 *  - removal.
 *
 * It is used by @ref Block as a base abstraction.
 *
 * @tparam Key  Type of the key (e.g. LhaID).
 * @tparam Data Stored data type (e.g. Parameter).
 */
template <typename Key, typename Data>
class IStorage {
public:
    virtual ~IStorage() = default;

    /**
     * @brief Stores a new element.
     *
     * Behavior when the key already exists depends on the implementation
     * (it may overwrite, ignore, or log an error).
     *
     * @param key  Key associated with the data.
     * @param data Shared pointer to the data object.
     */
    virtual void store(const Key& key, std::shared_ptr<Data> data) = 0;

    /**
     * @brief Assigns/replaces the data for an existing key.
     *
     * Implementations typically expect the key to exist and may error
     * or log if it does not.
     *
     * @param key  Key to update.
     * @param data Shared pointer to the new data object.
     */
    virtual void assign(const Key& key, std::shared_ptr<Data> data) = 0;

    /**
     * @brief Stores or assigns depending on key presence.
     *
     * Typically:
     *  - if key exists → behaves like @ref assign(),
     *  - otherwise     → behaves like @ref store().
     *
     * @param key  Key to store or update.
     * @param data Shared pointer to the data object.
     */
    virtual void store_or_assign(const Key& key, std::shared_ptr<Data> data) = 0;

    /**
     * @brief Checks if a given key is present in the storage.
     * @param key Key to test.
     * @return True if the key exists, false otherwise.
     */
    virtual bool contains(const Key& key) const = 0;

    /**
     * @brief Retrieves the data associated with a key.
     *
     * Semantics on non-existing keys are up to the implementation
     * (may throw, log, or return nullptr).
     *
     * @param key Key to retrieve.
     * @return Shared pointer to the stored data.
     */
    virtual std::shared_ptr<Data> retrieve(const Key& key) = 0;

    /**
     * @brief Removes the entry associated with a key.
     *
     * If the key does not exist, the behavior depends on the implementation
     * (may be a no-op, log, or throw).
     *
     * @param key Key to remove.
     */
    virtual void remove(const Key& key) = 0;
};


#endif // ISTORAGE_H
