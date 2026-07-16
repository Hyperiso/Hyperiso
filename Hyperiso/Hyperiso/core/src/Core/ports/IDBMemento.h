#ifndef IDBMEMENTO_H
#define IDBMEMENTO_H

#include <cstddef>
#include <memory>

#include "BlockAccessor.h"

/**
 * @file IDBMemento.h
 * @brief Interface for snapshotting and restoring BlockAccessor states.
 *
 * This header declares the IDBMemento interface, which provides a generic API
 * for:
 *  - taking snapshots of a BlockAccessor (parameter blocks),
 *  - restoring one or more steps back in the snapshot history,
 *  - inspecting stored snapshots.
 *
 * Concrete implementations (e.g. DBMemento) decide how snapshots are stored
 * in memory (stack, ring buffer, etc.).
 */

/**
 * @class IDBMemento
 * @brief Abstract interface for memento-like management of BlockAccessor states.
 *
 * This interface provides a minimal set of operations to:
 *  - capture the current state of a BlockAccessor (takeSnapshot),
 *  - roll back a given number of steps (restore),
 *  - inspect the content of stored snapshots (print_snapshot_content),
 *  - query the number of available snapshots (stack_size).
 *
 * It is intended to be used by higher-level components such as MemoryManager
 * to manage reversible changes to parameter blocks.
 */
class IDBMemento {
public:
    /// Virtual destructor for proper polymorphic cleanup.
    virtual ~IDBMemento() = default;

    /**
     * @brief Takes a snapshot of the given BlockAccessor.
     *
     * The implementation decides whether to deep-copy the data or only
     * capture references. In typical usage (e.g. DBMemento), a deep copy
     * is taken so that later modifications do not affect the stored state.
     *
     * @param blocks Shared pointer to the BlockAccessor to snapshot.
     */
    virtual void takeSnapshot(std::shared_ptr<BlockAccessor> blocks) = 0;

    /**
     * @brief Restores the state n_steps snapshots back.
     *
     * Implementations typically:
     *  - pop n_steps snapshots from an internal stack,
     *  - overwrite the current BlockAccessor state with the last popped one.
     *
     * If n_steps exceeds the number of stored snapshots, behavior is
     * implementation-defined (often a no-op or clamped restore).
     *
     * @param n_steps Number of snapshots to roll back (default 1).
     */
    virtual void restore(std::size_t n_steps = 1) = 0;

    /**
     * @brief Prints the content of a stored snapshot.
     *
     * Mainly intended for debugging. The exact format is implementation-
     * dependent. If n exceeds the number of snapshots, behavior is
     * implementation-defined (often using the most recent snapshot).
     *
     * @param n Index of the snapshot to print (0 can mean most recent).
     */
    virtual void print_snapshot_content(std::size_t n = 0) = 0;

    /**
     * @brief Returns the number of snapshots currently stored.
     *
     * @return Size of the internal snapshot stack.
     */
    virtual std::size_t stack_size() = 0;
};

#endif // IDBMEMENTO_H
