#ifndef __ISTORAGE_H__
#define __ISTORAGE_H__

template <typename Key, typename Data>
class IStorage {
public:
    virtual ~IStorage() = default;

    virtual void store(const Key& key, std::shared_ptr<Data> data) = 0;
    virtual void assign(const Key& key, std::shared_ptr<Data> data) = 0;
    virtual void store_or_assign(const Key& key, std::shared_ptr<Data> data) = 0;
    virtual bool contains(const Key& key) const = 0;
    virtual std::shared_ptr<Data> retrieve(const Key& key) = 0;
    virtual void remove(const Key& key) = 0;
};


#endif // __ISTORAGE_H__
