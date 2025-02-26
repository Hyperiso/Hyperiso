#include "DBMemento.h"

void DBMemento::takeSnapshot(std::shared_ptr<BlockAccessor> blocks) {
    DBMemento::snapshots.push(blocks);
}

void DBMemento::restore(size_t n) {
    auto head = DBMemento::snapshots.top();
    DBMemento::snapshots.pop();

    for (size_t i = 0; i < n - 1; i++)
        DBMemento::snapshots.pop();

    auto source = DBMemento::snapshots.top();
    DBMemento::snapshots.pop();
    DBMemento::Overwrite(head, source);
    DBMemento::snapshots.push(head);
}

void DBMemento::Overwrite(std::shared_ptr<BlockAccessor> &reciever, std::shared_ptr<BlockAccessor> source) {
    for (const auto& block_name : reciever->get_block_names()) {
        if (!source->has_block(block_name)) {
            reciever->remove_block(block_name);
        } else {
            for (const auto& id : reciever->get_block(block_name)->getAllIDs()) {
                if (!source->exist(block_name, id)) {
                    reciever->remove_item(block_name, id);
                } else {
                    reciever->setParameter(block_name, id, source->getParameter(block_name, id));
                }
            }
        }
    }
}
