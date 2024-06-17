function(target_set_warnings TARGET ENABLE ENABLED_AS_ERRORS)
    if (NOT ${ENABLE})
        message(STATUS "Warnings disabled for: ${TARGET}")
        return()
    endif()

    set(CLANG_WARNINGS
        -Wall
        -Wextra
        -Wpedantic)
    set(GCC_WARNINGS ${CLANG_WARNINGS})

    if(${ENABLED_AS_ERRORS})
        set(CLANG_WARNINGS ${CLANG_WARNINGS} -Werror)
        set(GCC_WARNINGS ${GCC_WARNINGS} -Werror)
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "CLANG")
        set(WARNINGS ${CLANG_WARNINGS})
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GCC")
        set(WARNINGS ${GCC_WARNINGS})
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(WARNINGS ${GCC_WARNINGS})
    else()
    endif()

    target_compile_options(${TARGET} PRIVATE ${WARNINGS})
    message(STATUS ${WARNINGS})
endfunction()