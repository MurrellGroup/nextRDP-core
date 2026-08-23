cmake_minimum_required(VERSION 3.20)
project(next_rdp_core_parity LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(rdp_dna5_exact STATIC
    src/msvcrt_compat.cpp
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
    src/distance_state.cpp
    src/main.cpp
    src/scan_state.cpp
    src/tree_state.cpp
    src/xover_state.cpp
)
target_include_directories(rdp-core PRIVATE include vendor/dna5 tests)
target_link_libraries(rdp-core PRIVATE rdp_dna5_exact)

if(EMSCRIPTEN)
    # The pinned 32-bit build evaluates source float expressions through the
    # x87 extended-precision stack. This preserves that compiler ABI while
    # leaving the supplied routine bodies byte-identical.
    target_compile_options(rdp_dna5_exact PRIVATE
        "-ffp-eval-method=extended"
    )
    # VB6 and the legacy DNA.dll source use the same 32-bit x87 evaluation
    # model in the pinned oracle. The orchestration port contains scalar
    # bodies from both, so it must use that ABI as well.
    target_compile_options(rdp-core PRIVATE
        "-ffp-eval-method=extended"
    )
    set_target_properties(rdp-core PROPERTIES SUFFIX ".js")
    target_link_options(rdp-core PRIVATE
        "-sENVIRONMENT=node"
        "-sALLOW_MEMORY_GROWTH=1"
        "-sEXIT_RUNTIME=1"
        "-sNODERAWFS=1"
    )
endif()
