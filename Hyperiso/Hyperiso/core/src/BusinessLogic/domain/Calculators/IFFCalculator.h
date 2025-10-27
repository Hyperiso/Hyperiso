#ifndef __IFFCALCULATOR_H__
#define __IFFCALCULATOR_H__

template<typename FFType>
class IFFCalculator {
public:
    virtual ~IFFCalculator() = default;	
    virtual double get(FFType, double) = 0;
};

#endif // __IFFCALCULATOR_H__
