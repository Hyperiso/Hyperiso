find_package(Doxygen)

if (DOXYGEN_FOUND)
    message(STATUS "WAOUW")
    add_custom_target(
        docs
        ${DOXYGEN_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Docs  
    )
endif()