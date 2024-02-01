#ifndef HYPERISO_LHA_ELEMENTS_H
#define HYPERISO_LHA_ELEMENTS_H

#include <string>
#include <optional>
#include <vector>
#include <memory>
#include <map>
#include "../Core/Logger.h"

enum class RenormalizationScheme {
    POLE, MSBAR, DRBAR, ONE_S, KIN, INVARIANT, MOM, SMOM, NONE
};

const std::map<std::string, int> SCALE_BLOCKS = {{"FMASS", 3}, {"FCONST", 4}, {"FCONSTRATIO", 6}, {"FBAG", 4}, {"FWCOEF", 0}, {"IMFWCOEF", 0}, {"FOBS", 3}, {"FOBSSM", 3}, {"FOBSERR", 3}};
const std::map<std::string, int> SCHEME_BLOCKS = {{"FMASS", 2}, {"FCONST", 3}, {"FCONSTRATIO", 5}, {"FBAG", 3}};
const std::map<std::string, int> VALUE_POS = {{"FCONST", 2}, {"FCONSTRATIO", 4}, {"FBAG", 2}, {"FWCOEF", 5}, {"IMFWCOEF", 5}, {"FOBS", 2}, {"FOBSSM", 2}, {"FOBSERR", 2}, {"FDIPOLE", 3}};

class AbstractElement {
protected:
    const std::string id;

public:
    inline explicit AbstractElement(const std::string& id) : id(id) {}

    inline std::string getId() const { return this->id; }
    virtual std::string toString() const = 0;
};

template <typename T>
class LhaElement : public AbstractElement {
private:
    const std::string block;
    T value;
    std::optional<RenormalizationScheme> rScheme;
    std::optional<double> Q;

    std::string encodeId(const std::string& block, const std::vector<std::string>& line);

public:
    LhaElement(const std::string& block, const std::vector<std::string>& line);

    inline std::string getBlock() const { return this->block; }
    inline T getValue() const { return this->value; }
    inline RenormalizationScheme getScheme() const { return this->rScheme.has_value() ? this->rScheme.value() : RenormalizationScheme::NONE; }
    inline double getScale() const { return this->Q.has_value() ? this->Q.value() : 0; }

    std::string toString() const override;
};

class LhaElementFactory {
public:
    static std::unique_ptr<AbstractElement> createElement(const std::string& blockName, const std::vector<std::string>& line);
};

#endif // HYPERISO_LHA_ELEMENTS_H
