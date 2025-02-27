#ifndef __DBMEMENTO_H__
#define __DBMEMENTO_H__

#include <stack>
#include <memory>
#include "BlockAccessor.h"

class DBMemento {
public:
    void takeSnapshot(std::shared_ptr<BlockAccessor> blocks);
    void restore(size_t n_steps = 1);

private:
    static void Overwrite(std::shared_ptr<BlockAccessor>& reciever, std::shared_ptr<BlockAccessor> source);

    inline static std::stack<std::shared_ptr<BlockAccessor>> snapshots;
};


#endif // __DBMEMENTO_H__
