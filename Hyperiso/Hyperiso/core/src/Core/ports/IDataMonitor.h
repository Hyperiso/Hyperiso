#ifndef __IDATAMONITOR_H__
#define __IDATAMONITOR_H__

#include "General.h"
#include <unordered_set>
#include <map>
#include <vector>

class IDataMonitor {
public:
    virtual ~IDataMonitor() = default;

    virtual std::unordered_set<std::string> get_all_blocks() = 0;
    virtual std::unordered_set<std::string> get_blocks_list(ParameterType param_type = ParameterType::SM) = 0;
    virtual std::map<LhaID, double> get_block_infos(const std::string& block, ParameterType param_type = ParameterType::SM) = 0;
    virtual std::vector<ParameterType> get_type_of_block(const std::string& block) = 0;
};


#endif // __IDATAMONITOR_H__
