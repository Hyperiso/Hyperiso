#ifndef __PARAMETERROUTER_H__
#define __PARAMETERROUTER_H__

#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#include "General.h"
#include "MemoryManager.h"
#include "Utils.h"

struct ParameterBlockRepartition {
    static inline const std::map<ParameterType, std::unordered_set<std::string>> BLOCKS {
        {ParameterType::SM, {"SMINPUTS", "MASS", "VCKMIN", "UPMNSIN", "GAUGE", "RECKM", "IMCKM", "UPMNS", "IMUPMNS"}},
        {ParameterType::SUSY, {"MASS", "HMIX", "ALPHA", "MSOFT", "NMIX", "UMIX", "VMIX", "A0MIX", "H0MIX", "STOPMIX", "SBOTMIX", "STAUMIX", "AU", "AD", "AE", "YU", "YD", "YE"}},
        {ParameterType::THDM, {"MASS", "ALPHA", "UCOUPL", "DCOUPL", "LCOUPL"}},
        {ParameterType::FLAVOR, {"FMASS", "FLIFE", "FCONST", "FCONSTRATIO", "FBAG", "FPARAM"}},
        {ParameterType::WILSON, {"FWCOEF", "IMFWCOEF"}},
        {ParameterType::DECAY, {"B_Ks", "B_ll", "B_Xs", "B_Dlnu", "B_Dslnu"}},
        {ParameterType::OBSERVABLE, {"FOBS", "FOBSERR", "FOBSSM", "FOBSSMERR", "FDIPOLE"}},
        {ParameterType::PASSTHROUGH, {"MODSEL", "SPINFO", "FMODSEL", "FCINFO", "MINPAR", "EXTPAR"}}
    };

    static std::vector<std::string> filter_custom_blocks(const std::vector<std::string>& source);
};

struct ParametersAccessRights {
    static inline const std::map<std::string, std::unordered_set<long>> SM_RIGHTS {
        {"MASS", {1, 2, 3, 4, 11, 12, 13, 14, 15, 16, 21, 22, 24, 25}}, 
        {"GAUGE", {1, 2, 3}},
    };

    static inline const std::map<std::string, std::unordered_set<long>> THDM_RIGHTS {
        {"MASS", {25, 35, 36, 37}}, 
        {"GAUGE", {{}}},
    };

    static inline const std::map<std::string, std::unordered_set<long>> SUSY_RIGHTS {
        {"MASS", {25, 35, 36, 37, 
                  1000001, 1000002, 1000003, 1000004, 1000005, 1000006, 1000011, 1000012, 1000013, 1000014, 1000015, 1000016, 
                  2000001, 2000002, 2000003, 2000004, 2000005, 2000006, 2000011, 2000013, 2000015, 
                  1000021, 1000022, 1000023, 1000024, 1000025, 1000035, 1000037, 1000039}}, 
        {"GAUGE", {{}}},
    };
};

class ParamRouter {
public:
    static ParameterType GetType(std::string block, LhaID id);
    static std::unordered_set<std::string> GetOwnedBlocks(ParameterType ptype);
   
};

#endif // __PARAMETERROUTER_H__
