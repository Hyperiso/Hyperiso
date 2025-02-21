#ifndef HYPERISO_LHA_READER_H
#define HYPERISO_LHA_READER_H

#include <string>
#include <map>
#include <memory>
#include <filesystem>
#include <climits>
#include "lha_blocks.h"
#include "lha_elements.h"
#include "lha_parser.h"

/**
 * @class LhaReader
 * @brief Reads and manages LHA blocks and elements from LHA or FLHA files.
 */
class LhaReader {
private:
    std::vector<Prototype> blockPrototypes;                     /**< List of block prototypes used for parsing. */
    std::map<std::string, std::shared_ptr<LhaBlock>> blocks;    /**< Map of block names to LhaBlock instances. */
    std::map<std::string, std::vector<std::vector<std::string>>> rawBlocks;
    bool isFLHA = false;                                        /**< Flag indicating if the file is in FLHA format. */
    std::filesystem::path lhaFile;                              /**< Path to the LHA file being read. */

    /**
     * @brief Parses the source string into blocks, optionally including comments.
     * @param comments If `true`, includes comments in the parsing process.
     */
    void parse_tokens(std::vector<Token> tokens, bool comments = false);

    /**
     * @brief Adds a new block to the reader from parsed lines.
     * @param id Block identifier.
     * @param lines Vector of lines containing the block's data.
     */
    void addBlock(const std::string& id, const std::vector<std::vector<std::string>>& lines);

public:
    /**
     * @brief Constructs an LhaReader with a specified file path.
     * @param path Path to the LHA file to read.
     */
    LhaReader(std::string_view path);
    
    /**
     * @brief Checks if a block exists in the reader.
     * @param id Block identifier.
     * @return `true` if the block exists, otherwise `false`.
     */
    bool hasBlock(const std::string& id) const;

    /**
     * @brief Checks if an element exists in the reader.
     * @param block_id Block identifier.
     * @param elt_id Block identifier.
     * @return `true` if the block exists, otherwise `false`.
     */
    bool hasElement(const std::string& block_id, const LhaID& elt_id) const;

    bool hasElement(const std::string& block_id, const std::vector<int>& elt_id) const { return hasElement(block_id, LhaID(elt_id)); }

    /**
     * @brief Reads all blocks from the LHA file.
     */
    void readAll();

    /**
     * @brief Finds a prototype by block name.
     * @param name Name of the block.
     * @return `Prototype` instance if found; otherwise, an empty prototype.
     */
    Prototype findPrototype(std::string name) const;

    /**
     * @brief Gets the file path of the LHA file.
     * @return String representing the path to the LHA file.
     */
    std::string getLhaPath() const;

    /**
     * @brief Updates the LHA file path and reloads all blocks.
     * @param newLha New file path for the LHA file.
     */
    void update(std::string_view newLha);

    /**
     * @brief Adds a new block type prototype to the reader.
     * @param blockName Name of the block.
     * @param itemCount Number of items in the block.
     * @param valueIdx Index of the value column.
     * @param scaleIdx Index of the scale column.
     * @param rgIdx Index of the renormalization group column.
     * @param globalScale Flag indicating if the block uses a global scale.
     */
    inline void addBlockType(std::string blockName, int itemCount=2, int valueIdx=1, int scaleIdx=-1, int rgIdx=-1, bool globalScale=false) {
        std::transform(blockName.begin(), blockName.end(), blockName.begin(), ::toupper);  // Make sure block name is uppercase 
        this->blockPrototypes.emplace_back(Prototype{blockName, itemCount, valueIdx, scaleIdx, rgIdx, globalScale});
    }

    
    template <typename T>
    inline void extractFromBlock(std::string blockName, std::vector<T*>& vars) {
        LhaBlock* block = this->getBlock(blockName);
        if (block) {
            for (size_t id=0; id!=vars.size(); ++id) {
                auto e = block->get(id + 1);
                *(vars.at(id)) = e ? static_cast<LhaElement<T>*>(e)->getValue() : T {};
            }
        }
    }

    template <typename T>
    inline void extractFromBlock(std::string blockName, std::vector<T>& vars) {
        LhaBlock* block = this->getBlock(blockName);
        if (block) {
            for (size_t id=0; id < vars.size(); ++id) {
                auto e = block->get(id + 1);
                vars.at(id) = e ? static_cast<LhaElement<T>*>(e)->getValue() : T {};
            }
        }
    }

    template <typename T>
    inline void extractFromBlock(std::string blockName, std::vector<T>& vars, const std::vector<int>& ids) {
        LhaBlock* block = this->getBlock(blockName);
        if (block) {
            for (size_t i=0; i < vars.size(); ++i) {
                auto e = block->get(ids.at(i));
                vars.at(i) = e ? static_cast<LhaElement<T>*>(e)->getValue() : T {};
            }
        }
    }

    template <typename T>
    inline void extractFromBlock(std::string blockName, std::vector<T>& vars, std::vector<LhaID>& ids) {
        LhaBlock* block = this->getBlock(blockName);
        if (block) {
            for (size_t i=0; i < vars.size(); ++i) {
                auto e = block->get(ids.at(i));
                vars[i] = e ? static_cast<LhaElement<T>*>(e)->getValue() : T {};
            }
        }
    }

    /**
     * @brief Retrieves a block by its identifier.
     * @param id Block identifier.
     * @return Pointer to the `LhaBlock` if it exists; otherwise, `nullptr`.
     */
    inline LhaBlock* getBlock(std::string id) const {
        std::transform(id.begin(), id.end(), id.begin(), ::toupper);
        return this->hasBlock(id) ? blocks.at(id).get() : nullptr;
    }

    /**
     * @brief Retrieves a value from a specified block and element.
     * @tparam T Type of the value to retrieve.
     * @param blockName Name of the block.
     * @param eltId Identifier of the element.
     * @return Value of the specified type.
     */
    template <typename T>
    inline T getValue(const std::string& blockName, const LhaID& eltId) {
        auto block = this->getBlock(blockName);
        if (block) {
            auto elt = block->get(eltId);
            if (elt) {
                return static_cast<LhaElement<T>*>(elt)->getValue();
            } else {
                LOG_ERROR("LHAReader", "Trying to access undefined element", eltId, "in block", blockName); 
            }   
        } else {
            LOG_ERROR("LHAReader", "Trying to access undefined block", blockName); 
        }
        
    }

    template <typename T>
    inline T getValue(const std::string& blockName, const std::vector<int>& eltId) {
        return getValue<T>(blockName, LhaID(eltId));
    }

    /**
     * @brief Gets the total number of blocks in the reader.
     * @return Number of blocks.
     */
    inline int getBlockCount() const {
        return blocks.size();
    }

    /**
     * @brief Retrieves the names of all blocks.
     * @return Vector of block names.
     */
    inline std::vector<std::string> getBlocksNames() const {
        std::vector<std::string> temp;
        for (auto &elem : blocks) {
            temp.push_back(elem.first);
        }
        return temp;
    }

    /**
     * @brief Converts the reader's content to a string representation.
     * @return String representing all blocks and their contents.
     */
    std::string toString() const;
};


#endif // HYPERISO_LHA_READER_H