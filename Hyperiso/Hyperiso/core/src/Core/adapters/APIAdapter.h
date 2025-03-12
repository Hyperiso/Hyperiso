#ifndef __APIADAPTER_H__
#define __APIADAPTER_H__

#include "IMonitor.h"
#include "IPathProvider.h"
#include "IDataMonitor.h"
#include "MemoryManager.h"

enum class APIPath {
    LHA_PATH,
};

class APIAdapter : public IMonitor<ExternalFlag>, IPathProvider<APIPath>, IDataMonitor {
public:
    bool check_flag(ExternalFlag flag);
    fs::path get_path(APIPath path_name);
    std::unordered_set<std::string> get_all_blocks();
    std::unordered_set<std::string> get_blocks_list(ParameterType param_type = ParameterType::SM);
    std::map<LhaID, double> get_block_infos(const std::string& block, ParameterType param_type = ParameterType::SM);
    std::vector<ParameterType> get_type_of_block(const std::string& block);
};


#endif // __APIADAPTER_H__
