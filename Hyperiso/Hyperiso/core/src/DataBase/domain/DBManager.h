#ifndef DBMANAGER_H
#define DBMANAGER_H

#include <memory>
#include <filesystem>

#include "DBNode.h"
#include "ParserFactory.h"
#include "LhaBlockPrototype.h"

namespace fs = std::filesystem;

/**
 * @file DBManager.h
 * @brief High-level interface to load and save Node trees from/to various file formats.
 *
 * DBManager provides a unified entry point for:
 *   - detecting the appropriate parser based on the file extension,
 *   - running basic file validity checks (existence, non-emptiness, encoding),
 *   - configuring LHA/FLHA block prototypes when needed,
 *   - sanitizing the resulting Node tree.
 *
 * Supported formats are determined by ParserFactory::Type and the EXTENSIONS table:
 *   - JSON (.json)
 *   - YAML (.yml, .yaml)
 *   - LHA / SLHA / FLHA (.lha, .slha, .flha)
 */

/**
 * @class DBManager
 * @brief Front-end for reading and writing structured data as Node trees.
 *
 * A DBManager instance encapsulates:
 *   - a static mapping between parser types and supported file extensions,
 *   - a static set of LHA/SLHA/FLHA block prototypes used by LhaParser,
 *   - helper methods to sanitize files and Node trees.
 *
 * Typical usage:
 * @code
 *   DBManager db;
 *   auto root = db.read_from_file("input.flha");
 *   // ... work with root ...
 *   db.write_to_file("output.yaml", root);
 * @endcode
 */
class DBManager {
private:
    /**
     * @brief Mapping between parser type and supported filename extensions.
     *
     * This table is used by deduce_parser_type() to choose the correct
     * parser implementation from ParserFactory.
     */
    static inline const std::map<ParserFactory::Type, std::unordered_set<std::string>> EXTENSIONS {
        {ParserFactory::Type::JSON, {".json"}},
        {ParserFactory::Type::YAML, {".yml", ".yaml"}},
        {ParserFactory::Type::LHA, {".lha", ".slha", ".flha"}},
    }; 

    /**
     * @brief Set of LHA/SLHA/FLHA block prototypes used when parsing LHA files.
     *
     * This set is:
     *   - automatically populated with default prototypes based on the file
     *     extension (LHA_BLOCKS, SLHA_BLOCKS, FLHA_BLOCKS),
     *   - extensible at runtime via add_lha_prototype().
     */
    static inline std::unordered_set<Prototype> lha_prototypes {};

    /**
     * @brief Determines which parser type should be used for a given file.
     *
     * The decision is based solely on the file extension and the EXTENSIONS
     * lookup table. If no matching extension is found, a fatal error is
     * logged via LOG_ERROR.
     *
     * @param file_path Path to the file whose type should be deduced.
     * @return ParserFactory::Type corresponding to the file extension.
     */
    ParserFactory::Type deduce_parser_type(fs::path file_path);

    /**
     * @brief Ensures that default LHA block prototypes are loaded.
     *
     * This function:
     *   - always inserts all entries from LHA_BLOCKS,
     *   - for ".slha" or ".lha" also inserts SLHA_BLOCKS,
     *   - for ".flha" inserts FLHA_BLOCKS,
     *   - does not overwrite existing prototypes; it only inserts missing ones.
     *
     * @param file_path Path to the file being parsed (used for its extension).
     */
    void add_default_lha_prototypes(fs::path file_path);

    /**
     * @brief Performs basic validation on an input or output file.
     *
     * The checks include:
     *   - existence,
     *   - non-empty content,
     *   - ability to open the file,
     *   - absence of suspicious control characters (aside from \n, \r, \t).
     *
     * @param file_path Path of the file to inspect.
     *
     * @throws std::runtime_error if one of the checks fails.
     */
    void sanitize_file(const fs::path& file_path);

    /**
     * @brief Performs basic validity checks on a Node tree.
     *
     * This helper ensures:
     *   - the root Node is non-null,
     *   - the root contains at least one child,
     *   - all keys in @p required_keys are present at the root level.
     *
     * @param root          Root Node to check.
     * @param required_keys Optional list of keys that must be present.
     *
     * @throws std::runtime_error if the tree is null, empty, or missing
     *         any required key.
     */
    void sanitize_tree(const std::shared_ptr<Node>& root, const std::vector<BlockName>& required_keys = {});

public:
    /**
     * @brief Reads a structured file and returns its Node representation.
     *
     * The pipeline is:
     *   - sanitize_file(file_path),
     *   - deduce_parser_type(file_path),
     *   - create an appropriate parser via ParserFactory,
     *   - for LHA files: add default LHA prototypes and configure LhaParser,
     *   - parser->readFromFile(file_path),
     *   - sanitize_tree(root).
     *
     * @param file_path Path to the file to be read.
     * @return Shared pointer to the root Node.
     *
     * @throws std::runtime_error for file errors, parser errors, or tree
     *         validation errors.
     */
    std::shared_ptr<Node> read_from_file(fs::path file_path);

    /**
     * @brief Writes a Node tree to disk using the appropriate format.
     *
     * The pipeline is:
     *   - sanitize_tree(root, {"SMINPUTS"}) to ensure minimal consistency,
     *   - deduce_parser_type(file_path),
     *   - create an appropriate parser,
     *   - parser->writeToFile(file_path, root),
     *   - sanitize_file(file_path) to validate the written file.
     *
     * Note: requiring "SMINPUTS" is domain-specific and encodes the
     * expectation that any physical input set at least carries this block.
     *
     * @param file_path Path where the file should be written.
     * @param root      Root Node to serialize.
     *
     * @throws std::runtime_error if validation fails or writing fails.
     */
    void write_to_file(fs::path file_path, std::shared_ptr<Node> root);

    /**
     * @brief Adds or overrides an LHA block prototype used by LhaParser.
     *
     * This function lets the user extend or customize the interpretation
     * of LHA/SLHA/FLHA blocks at runtime. The block name is normalized
     * to upper-case before insertion.
     *
     * If a prototype with the same block name already exists:
     *   - if it is identical to @p new_prototype, a warning is logged and
     *     the call is ignored,
     *   - otherwise, a fatal error is logged (conflicting layouts).
     *
     * @param blockName   Name of the block (case-insensitive, will be uppercased).
     * @param itemCount   Total number of columns in the block (default: 2).
     * @param valueIdx    Index of the value column (default: 1).
     * @param scaleIdx    Index of the scale column, or -1 if none (default: -1).
     * @param rgIdx       Index of the renormalization scheme column, or -1 (default: -1).
     * @param binIdx      Index of the bin lower-bound column, or -1 if unbinned (default: -1).
     * @param globalScale True if the block uses a global Q= scale in the header.
     */
    void add_lha_prototype(BlockName blockName, size_t itemCount=2, size_t valueIdx=1, int scaleIdx=-1, int rgIdx=-1, int binIdx=-1, bool globalScale=false);

    friend class NodeProviderFactory;
};


#endif // DBMANAGER_H
