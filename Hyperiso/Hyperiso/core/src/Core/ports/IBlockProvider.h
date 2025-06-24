#ifndef IBLOCK_PROVIDER_H
#define IBLOCK_PROVIDER_H

template<typename T, typename U>
class IBlockProvider {
public:
    virtual bool exists(U blockname, T) = 0;
    virtual void log_all_blocks(T type) = 0;
    virtual void log_block(T type, U blockname) = 0;

    
};

#endif