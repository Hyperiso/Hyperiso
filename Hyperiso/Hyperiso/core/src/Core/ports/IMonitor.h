#ifndef __IMONITOR_H__
#define __IMONITOR_H__

#include "MemoryManager.h"

template<typename FlagType>
class IMonitor {
public:
    virtual ~IMonitor() = default;

    virtual bool check_flag(FlagType flag) = 0;
};

#endif // __IMONITOR_H__
