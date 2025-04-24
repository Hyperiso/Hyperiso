#ifndef __IPARAMMODIFIER_H__
#define __IPARAMMODIFIER_H__

template <typename T, typename U, typename V>
class IDataMutator {
public:
    virtual ~IDataMutator() = default;

    virtual void mutate(const T&, U) = 0;
    virtual void change_mode(const T&, V) = 0;
};

#endif // __IPARAMMODIFIER_H__
