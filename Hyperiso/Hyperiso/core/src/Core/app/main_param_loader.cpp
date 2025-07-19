#include "ParamBlockLoader.h"

int main(){

    ParamBlockLoader loader;
    ba = std::make_shared<BlockAccessor>();
    loader.load(ba,"Assets/default/lha/testInput.flha");
    LOG_INFO("Loaded blocks:", ba->size());
    loader.save("Assets/default/lha/testOutput.flha", ba);
    LOG_INFO("Saved blocks to Assets/default/lha/testOutput.flha");
    
    return 0;
}