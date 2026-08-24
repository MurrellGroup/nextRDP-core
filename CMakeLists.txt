cmake_minimum_required(VERSION 3.20)
project(next_rdp_core LANGUAGES CXX)

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

add_library(next_rdp_core STATIC
    src/analysis.cpp
    src/distance_state.cpp
    src/identification_state.cpp
    src/legacy_method_state.cpp
    src/mutation_state.cpp
    src/rescan_schedule.cpp
    src/round_state.cpp
    src/scan_state.cpp
    src/tree_state.cpp
    src/three_seq_state.cpp
    src/xover_state.cpp
)
target_include_directories(next_rdp_core
    PUBLIC src
    PRIVATE include vendor/dna5
)
target_link_libraries(next_rdp_core PUBLIC rdp_dna5_exact)

add_executable(rdp-core src/main.cpp)
target_include_directories(rdp-core PRIVATE include vendor/dna5 tests)
target_link_libraries(rdp-core PRIVATE next_rdp_core)

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
    target_compile_options(next_rdp_core PRIVATE
        "-ffp-eval-method=extended"
    )
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

    add_executable(next-rdp-core-web src/web_api.cpp)
    target_link_libraries(next-rdp-core-web PRIVATE next_rdp_core)
    target_compile_options(next-rdp-core-web PRIVATE
        "-ffp-eval-method=extended"
    )
    set_target_properties(next-rdp-core-web PROPERTIES SUFFIX ".mjs")
    target_link_options(next-rdp-core-web PRIVATE
        "--no-entry"
        "-sENVIRONMENT=web,worker,node"
        "-sMODULARIZE=1"
        "-sEXPORT_ES6=1"
        "-sEXPORT_NAME=createNextRdpCore"
        "-sALLOW_MEMORY_GROWTH=1"
        "-sEXPORTED_FUNCTIONS=['_malloc','_free','_next_rdp_analyze','_next_rdp_result_json','_next_rdp_error','_next_rdp_version']"
        "-sEXPORTED_RUNTIME_METHODS=['UTF8ToString','HEAPU8']"
    )
endif()
