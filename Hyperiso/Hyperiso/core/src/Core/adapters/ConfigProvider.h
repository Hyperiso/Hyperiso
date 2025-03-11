#ifndef __CONFIGPROVIDER_H__
#define __CONFIGPROVIDER_H__

#include "IMonitor.h"

class ConfigProvider : public IMonitor<ExternalFlag> {
public:
    bool check_flag(ExternalFlag flag);
    Model get_model();
};

#endif // __CONFIGPROVIDER_H__
