#ifndef HYPERISO_LHA_ELEMENTS_H
#define HYPERISO_LHA_ELEMENTS_H

#include <string>
#include <optional>
#include <vector>
#include <memory>
#include <map>
#include "Logger.h"

class LhaBlock;

enum class RenormalizationScheme {
    POLE, MSBAR, DRBAR, ONE_S, KIN, INVARIANT, MOM, SMOM, NONE
};

class AbstractElement {
protected:
    const std::string id;
    const LhaBlock* block;

public:
    inline explicit AbstractElement(LhaBlock* block, const std::string& id) : block(block), id(id) {}

    inline std::string getId() const { return this->id; }
    virtual std::string toString() const = 0;
};

template <typename T>
class LhaElement : public AbstractElement {
private:
    T value;
    std::optional<RenormalizationScheme> rScheme;
    std::optional<double> Q;

    std::string encodeId(LhaBlock* block, const std::vector<std::string>& line);

public:
    LhaElement(LhaBlock* block, const std::vector<std::string>& line);

    inline LhaBlock* getBlock() const { return this->block; }
    inline T getValue() const { return this->value; }
    inline RenormalizationScheme getScheme() const { return this->rScheme.has_value() ? this->rScheme.value() : RenormalizationScheme::NONE; }
    inline double getScale() const { return this->Q.has_value() ? this->Q.value() : 0; }

    std::string toString() const override;
};

class LhaElementFactory {
public:
    static std::unique_ptr<AbstractElement> createElement(LhaBlock* block, const std::vector<std::string>& line);
};

#endif // HYPERISO_LHA_ELEMENTS_H
