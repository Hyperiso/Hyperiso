#ifndef HYPERISO_LHA_READER_H
#define HYPERISO_LHA_READER_H

#include <string>
#include <map>
#include <memory>
#include <filesystem>
#include <climits>
#include <regex>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "IParser.h"
#include "lha_blocks.h"
#include "lha_elements.h"
#include "Logger.h"

/**
 * @file LhaParser.h
 * @brief Parser for LHA / FLHA-style files into LhaBlock structures and Node trees.
 *
 * LhaParser decodes a textual LHA (or FLHA) file into:
 *   - a set of LhaBlock objects (one per block),
 *   - then a Node tree suitable for use by the rest of the framework.
 *
 * It uses:
 *   - a regex-based tokenizer,
 *   - a simple state machine to group tokens by block,
 *   - a set of Prototype descriptors to interpret each block's layout.
 */

/**
 * @enum TokenType
 * @brief Token categories used during LHA lexical analysis.
 *
 * The order of the enum values is tightly coupled to the capturing groups
 * of the regular expression @ref analyzer_rx. Changing one requires
 * updating the other accordingly.
 */
enum class TokenType {
    FLOAT,  /**< Floating point number. */
    INTEGER,/**< Integer number. */
    BLOCK,  /**< Block keyword. */
    DECAY,  /**< Decay keyword. */
    NEWLINE,/**< Newline character. */
    SKIP,   /**< Whitespace character. */
    COMMENT,/**< Comment line. */
    WORD,   /**< Word token. */
    OTHER   /**< Other types of tokens. */
};

/**
 * @var analyzer_rx
 * @brief Regular expression used to tokenize LHA/FLHA file content.
 *
 * The regex is structured into capturing groups that map directly to
 * values of TokenType:
 *   - FLOAT
 *   - INTEGER
 *   - BLOCK
 *   - DECAY
 *   - NEWLINE
 *   - SKIP
 *   - COMMENT
 *   - WORD
 *   - OTHER
 *
 * The lexer walks over all matches and selects the first matching
 * capturing group to decide the token type.
 */
const std::regex analyzer_rx(
    R"x(((?:[+-])?(?:\d+\.\d*|\.\d+)(?:[eEdD][+-]\d+)?|(\d+(?:[eEdD][+-]\d+)?))|(?:[+-]?\d+(?!\.))|(block)|(decay)|(\n)|([ \t]+)|(#.*)|([\w\=\.]+)|([^#]*))x",
    std::regex_constants::icase
); 

/**
 * @struct Token
 * @brief Single lexical token produced from an LHA file.
 *
 * A token is annotated with:
 *   - its type (TokenType),
 *   - its literal text,
 *   - its (row, col) logical position in the original input.
 */
struct Token {
    TokenType type;     /**< Type of the token. */
    std::string value;  /**< Text value of the token. */
    int row;            /**< Row number where the token is located. */
    int col;            /**< Column number where the token is located. */
};

/**
 * @class LhaParser
 * @brief Parser for LHA/FLHA files, producing blocks and a Node representation.
 *
 * LhaParser implements the IParser interface using the LHA/FLHA conventions:
 *   - BLOCK sections are recognized via the "Block" keyword,
 *   - DECAY sections are recognized via the "Decay" keyword (and currently skipped),
 *   - each block is associated to a Prototype (block layout),
 *   - each line in a block is converted into an LhaElement<T>,
 *   - all blocks are then converted into a top-level Node tree.
 *
 * The set of known prototypes can be configured via set_prototypes()
 * (e.g. using LHA_BLOCKS, SLHA_BLOCKS, FLHA_BLOCKS).
 */
class LhaParser : public IParser {
private:
    std::unordered_set<Prototype> blockPrototypes;                     /**< Set of block prototypes used for parsing. */

    /**
     * @brief Populate the internal prototype set from some LHA-related path.
     *
     * Currently unused in the provided implementation but reserved for
     * future extensions where prototypes could be discovered or configured
     * from an external directory or file.
     *
     * @param lha_path Path to an LHA configuration or reference file.
     */
    void fill_prototypes(std::string_view lha_path);

    /**
     * @brief Lexical analysis: splits the source string into tokens.
     *
     * Uses @ref analyzer_rx to iterate over the input and emit tokens,
     * tracking (row, col) coordinates. Whitespace tokens are skipped,
     * and NEWLINE tokens are used to advance the row counter.
     *
     * @param src Raw LHA/FLHA file content.
     * @return Vector of tokens in order of appearance.
     */

    std::vector<Token> tokenize(const std::string& src) const;

    /**
     * @brief High-level token parser: groups tokens into blocks and lines.
     *
     * The parser:
     *   - detects "Block" and "Decay" headers,
     *   - matches block names against known Prototypes,
     *   - skips unknown blocks and decay blocks (with log warnings),
     *   - handles optional global scale Q= in the header,
     *   - accumulates tokens into rows of strings for each block.
     *
     * @param tokens   Token stream to parse.
     * @param comments If true, comment tokens are kept as lines; otherwise they are ignored.
     * @return A map from block name to raw line tokens.
     */
    std::map<BlockName, std::vector<std::vector<std::string>>> parse_tokens(std::vector<Token> tokens, bool comments = false) const;

    /**
     * @brief Builds and inserts an LhaBlock from raw line data.
     *
     * This helper:
     *   - finds the appropriate Prototype for the block name,
     *   - constructs an LhaBlock,
     *   - calls LhaBlock::readData(lines),
     *   - normalizes the block name to upper-case,
     *   - stores the block in the provided map.
     *
     * @param blocks Map of blocks to populate.
     * @param id     Name/identifier of the block.
     * @param lines  Raw tokenized lines belonging to this block.
     */
    void addBlock(std::map<BlockName, std::shared_ptr<LhaBlock>>& blocks, const BlockName& id, const std::vector<std::vector<std::string>>& lines) const;

    /**
     * @brief Finds the Prototype associated with a given block name.
     *
     * The comparison is case-insensitive (BlockName is case-insensitive
     * internally). If no matching prototype is found, an "empty" Prototype
     * with blockName == "" is returned.
     *
     * @param name Block name to search for.
     * @return Prototype describing the block layout, or a default/empty
     *         Prototype if none is found.
     */
    Prototype findPrototype(BlockName name) const;

    /**
     * @brief Converts a set of parsed blocks into a top-level Node tree.
     *
     * For each block:
     *   - blockPtr->toDBNode() is called to obtain the block-level Node,
     *   - the internal group under the block name is copied to @p root,
     *   - if the block node exposes a "scale" entry, a top-level
     *     "<BLOCKNAME>.scale" is also added as a numeric field.
     *
     * @param blocks Map from block name to LhaBlock instances.
     * @return Shared pointer to the root Node representing the whole LHA file.
     */
    std::shared_ptr<Node> toDBNode(std::map<BlockName, std::shared_ptr<LhaBlock>> blocks) const;

public:
    /**
     * @brief Parses an LHA/FLHA file from disk and returns its Node representation.
     *
     * This implementation:
     *   - loads the file into memory,
     *   - delegates to parse(src).
     *
     * @param input_file Path to the LHA/FLHA file to read.
     * @return Shared pointer to the root Node representation.
     *
     * @throws std::runtime_error if the file cannot be opened or parsing fails.
     */
    std::shared_ptr<Node> readFromFile(const std::string& input_file) const override;

    /**
     * @brief Parses LHA/FLHA content from an in-memory string.
     *
     * The pipeline is:
     *   - tokenize(src),
     *   - parse_tokens(...),
     *   - construct LhaBlocks for all recognized blocks,
     *   - convert them to a Node via toDBNode().
     *
     * @param src Raw LHA/FLHA text.
     * @return Shared pointer to the root Node representation.
     *
     * @throws std::runtime_error on lexical or structural errors.
     */
    std::shared_ptr<Node> parse(const std::string& src) const override;

    /**
     * @brief Serializes a Node structure back to an LHA-like file.
     *
     * Currently marked as TODO in the implementation. The intended design
     * is to invert the toDBNode/LhaBlock/LhaElement pipeline in order to
     * reconstruct a valid LHA/FLHA text file from a Node tree.
     *
     * @param filename Target file name.
     * @param root     Root Node to serialize.
     */
    void writeToFile(const std::string& filename, const std::shared_ptr<Node>& root) const;

    /**
     * @brief Sets the set of block prototypes used by the parser.
     *
     * Typical usage is to pass one of the predefined sets:
     *   - LHA_BLOCKS,
     *   - SLHA_BLOCKS,
     *   - FLHA_BLOCKS,
     * or their union, depending on the type of file being parsed.
     *
     * @param prototypes Set of Prototype objects describing known blocks.
     */
    void set_prototypes(const std::unordered_set<Prototype>& prototypes);
};

#endif // HYPERISO_LHA_READER_H