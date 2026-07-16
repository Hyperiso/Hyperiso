#ifndef DBMEMENTO_H
#define DBMEMENTO_H

#include <stack>
#include <memory>

#include "IDBMemento.h"
#include "BlockAccessor.h"

/**
 * @file DBMemento.h
 * @brief Concrete memento for BlockAccessor snapshots.
 *
 * This header defines DBMemento, a stack-based implementation of the
 * IDBMemento interface. It stores successive copies of BlockAccessor
 * instances, allowing the system to:
 *  - save intermediate states (takeSnapshot),
 *  - roll back a configurable number of steps (restore),
 *  - introspect snapshots for debugging (print_snapshot_content).
 */

 /**
 * @class DBMemento
 * @brief Stack-based implementation of IDBMemento for BlockAccessor.
 *
 * DBMemento maintains an internal stack of BlockAccessor snapshots:
 * each call to takeSnapshot() pushes a new copy onto the stack, while
 * restore() pops one or more entries and overwrites the target
 * BlockAccessor with the stored state.
 *
 * This is used by MemoryManager to:
 *  - save the state after loading default inputs,
 *  - restore prior states when reloading user input or switching LHA files.
 *
 * @see IDBMemento
 * @see BlockAccessor
 * @see MemoryManager
 */
class DBMemento : public IDBMemento {
public:
    /**
     * @brief Takes a snapshot of the given BlockAccessor.
     *
     * The snapshot is stored internally and can be restored later with
     * restore().
     *
     * @param blocks Shared pointer to the BlockAccessor to snapshot.
     */
    void takeSnapshot(std::shared_ptr<BlockAccessor> blocks);

    /**
     * @brief Restores the state n_steps snapshots back.
     *
     * Pops up to @p n_steps snapshots from the internal stack and uses the
     * last popped snapshot to overwrite the current BlockAccessor state.
     *
     * @param n_steps Number of snapshots to roll back (default 1).
     */
    void restore(size_t n_steps = 1);

    /**
     * @brief Prints the content of a stored snapshot.
     *
     * Mainly for debugging. If @p n is 0, it usually refers to the most
     * recent snapshot (implementation dependent).
     *
     * @param n Index of the snapshot to print (default 0).
     */
    void print_snapshot_content(size_t n = 0);

    /**
     * @brief Returns the number of snapshots currently stored.
     *
     * @return Size of the snapshot stack.
     */
    size_t stack_size();

private:
    /**
     * @brief Helper to overwrite a receiver BlockAccessor with the source one.
     *
     * This is typically implemented as a deep copy of all blocks and
     * parameters.
     *
     * @param receiver Target BlockAccessor reference to overwrite.
     * @param source   Shared pointer to the source BlockAccessor snapshot.
     */
    static void Overwrite(std::shared_ptr<BlockAccessor>& receiver, std::shared_ptr<BlockAccessor> source);

    /// Internal stack of snapshots.
    std::stack<std::shared_ptr<BlockAccessor>> snapshots;
};


#endif // DBMEMENTO_H
