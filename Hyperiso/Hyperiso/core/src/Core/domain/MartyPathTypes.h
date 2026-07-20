#ifndef HYPERISO_MARTY_PATH_TYPES_H
#define HYPERISO_MARTY_PATH_TYPES_H

/**
 * @enum MartyPath
 * @brief Enumerates MARTY-related filesystem resources.
 *
 * MODEL_FILE: path to the MARTY model header file.
 * SM_MODEL_FILE: packaged Standard Model MARTY header used by the SM backend.
 * TEMPLATE_DIR: read-only MARTY template directory.
 * PARAM_MAPPING_DIR: read-only directory containing mapping JSON files.
 * SM_MAPPING_FILE: read-only Hyperiso/SM mapping file.
 * BSM_MAPPING_FILE: optional user-provided BSM mapping file.
 * MARTY_TEMP_DIR: writable generated-code/cache directory.
 */
enum class MartyPath {
    MODEL_FILE,
    SM_MODEL_FILE,
    TEMPLATE_DIR,
    PARAM_MAPPING_DIR,
    SM_MAPPING_FILE,
    BSM_MAPPING_FILE,
    MARTY_TEMP_DIR
};

#endif // HYPERISO_MARTY_PATH_TYPES_H
