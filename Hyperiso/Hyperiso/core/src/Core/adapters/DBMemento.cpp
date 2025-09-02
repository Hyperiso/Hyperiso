#include "DBMemento.h"

void DBMemento::takeSnapshot(std::shared_ptr<BlockAccessor> blocks) {
    snapshots.push(blocks);
}

void DBMemento::restore(size_t n) {
    auto head = snapshots.top();
    snapshots.pop();

    for (size_t i = 0; i < n - 1; i++)
        snapshots.pop();

    auto source = snapshots.top();
    snapshots.pop();
    DBMemento::Overwrite(head, source);
    snapshots.push(head);
}

void DBMemento::print_snapshot_content(size_t n) {
    if (n > snapshots.size() - 1) {
        LOG_ERROR("LogicError", "Please specify an existing snapshot (only", snapshots.size(), "stored, asked", n);
    }

    std::stack<std::shared_ptr<BlockAccessor>> temp;
    for (size_t i = 0; i < n; i++) {
        temp.push(snapshots.top());
        snapshots.pop();
    }

    LOG_INFO(snapshots.top());

    for (size_t i = 0; i < n; i++) {
        snapshots.push(temp.top());
        temp.pop();
    }
}

size_t DBMemento::stack_size() {
    return snapshots.size();
}

void DBMemento::Overwrite(std::shared_ptr<BlockAccessor> &reciever, std::shared_ptr<BlockAccessor> source) {
    for (const auto& block_name : reciever->get_block_names()) {
        if (!source->contains(block_name)) {
            reciever->erase(block_name);
            continue;
        } 
        for (const auto& id : reciever->at(block_name)->getAllIDs()) {
            if (!source->has_param(block_name, id)) {
                reciever->remove_item(block_name, id);
            } else {
                reciever->remove_item(block_name, id);
                reciever->setParameter(block_name, id, std::move(source->getParameter(block_name, id)));
            }
        }
    }
}
