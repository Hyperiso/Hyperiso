#ifndef __MARTYMODELPATHAPI_H__
#define __MARTYMODELPATHAPI_H__

#include "MartyAdapter.h"
#include "Config.h"
#include "ICoreAPI.h"

class MartyModelPathAPI : public ICoreAPI<fs::path> {
public:
    inline fs::path get() override { return MartyAdapter().get_path(MartyPath::MODEL_FILE); }
}; 

#endif // __MARTYMODELNAMEAPI_H__
