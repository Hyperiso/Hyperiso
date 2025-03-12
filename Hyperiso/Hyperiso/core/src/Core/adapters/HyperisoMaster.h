#ifndef __CONFIGPROVIDER_H__
#define __CONFIGPROVIDER_H__

#include "IMonitor.h"

class HyperisoMaster : public IMonitor<ExternalFlag> {
public:
    void init(const std::string &lhaFile, Config config);
    bool check_flag(ExternalFlag flag);
    Model get_model();
};

#endif // __CONFIGPROVIDER_H__
