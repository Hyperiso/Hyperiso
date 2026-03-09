#ifndef __IDEPENDENCYMOUNTER_H__
#define __IDEPENDENCYMOUNTER_H__

template <typename Type_T, typename Block_T, typename Id_T> 
class IDependencyPruner {
public:
    virtual ~IDependencyPruner() = default;

    virtual void reattach_block(Type_T, const Block_T&) = 0;
    virtual void detach_block(Type_T, const Block_T&) = 0;

    virtual void reattach_parameter(Type_T, const Block_T&, const Id_T&) = 0;
    virtual void detach_parameter(Type_T, const Block_T&, const Id_T&) = 0;
};

#endif // __IDEPENDENCYMOUNTER_H__
