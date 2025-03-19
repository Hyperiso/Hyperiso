#ifndef __ISTORAGE_H__
#define __ISTORAGE_H__

template <typename Key, typename Data>
class IStorage {
public:
    virtual ~IStorage() = default;

    virtual void store(const Key& key, Data&& data) = 0;
    virtual void remove(const Key& key) = 0;
    virtual Data& retrieve(const Key& key) = 0;
    virtual bool contains(const Key& key) const = 0;
    virtual void update(const Key& key, Parameter&& param) = 0;
};


#endif // __ISTORAGE_H__
