#ifndef __MARTYMODELNAMEAPI_H__
#define __MARTYMODELNAMEAPI_H__

#include "MartyAdapter.h"
#include "Config.h"
#include "ICoreAPI.h"

class MartyModelNameAPI : public ICoreAPI<std::string> {
public:
    inline std::string get() override { return MartyAdapter().get_marty_model_name(); }
}; 

#endif // __MARTYMODELNAMEAPI_H__
