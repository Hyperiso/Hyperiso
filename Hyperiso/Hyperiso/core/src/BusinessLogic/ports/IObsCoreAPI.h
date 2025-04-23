#ifndef ICORE_API_H
#define ICORE_API_H

template <typename T>
class IObsCoreAPI{
public:
    inline virtual T get() = 0;
};

#endif