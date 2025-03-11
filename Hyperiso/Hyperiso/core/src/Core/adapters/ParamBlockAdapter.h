#ifndef __BLOCKSCREATOR_H__
#define __BLOCKSCREATOR_H__

#include <memory>
#include "IDataLoader.h"
#include "BlockAccessor.h"
#include "NodeProviderFactory.h"

class ParamBlockLoader : public IDataLoader<BlockAccessor> {
public:
    void load(std::shared_ptr<BlockAccessor> dest, fs::path src_file) override;
};

#endif // __BLOCKSCREATOR_H__
