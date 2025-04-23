#ifndef IOBS_CORE_API_H
#define IOBS_CORE_API_H

template <typename T>
class IObsCoreAPI{
public:
    inline virtual T get() = 0;
};

#endif