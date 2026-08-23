cmake_minimum_required(VERSION 3.20)
project(next_rdp_core_parity LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(rdp_dna5_exact STATIC
    vendor/dna5/MathFuncsDll.cpp
)
target_compile_definitions(rdp_dna5_exact PRIVATE MATHFUNCSDLL_EXPORTS)
target_include_directories(rdp_dna5_exact
    PUBLIC vendor/dna5
    PRIVATE include
)
target_compile_options(rdp_dna5_exact PRIVATE
    -Wno-unknown-pragmas
    -Wno-deprecated-volatile
)

add_executable(rdp-core
    src/main.cpp
)
target_include_directories(rdp-core PRIVATE include vendor/dna5 tests)
target_link_libraries(rdp-core PRIVATE rdp_dna5_exact)

if(EMSCRIPTEN)
    set_target_properties(rdp-core PROPERTIES SUFFIX ".js")
    target_link_options(rdp-core PRIVATE
        "-sENVIRONMENT=node"
        "-sALLOW_MEMORY_GROWTH=1"
        "-sEXIT_RUNTIME=1"
        "-sNODERAWFS=1"
    )
endif()
