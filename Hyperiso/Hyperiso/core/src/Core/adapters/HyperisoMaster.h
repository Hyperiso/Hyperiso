#ifndef __HYPERISOMASTER_H__
#define __HYPERISOMASTER_H__

#include "MemoryManager.h"
#include "IMonitor.h"

class HyperisoMaster : public IMonitor {
private:

public:
    void init();
    Config getConfig();
};


#endif // __HYPERISOMASTER_H__
