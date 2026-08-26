cmake_minimum_required(VERSION 3.20)
project(next_rdp_core LANGUAGES CXX)

option(NEXT_RDP_ENABLE_OPENMP
    "Enable the source DLL's OpenMP runtime (not output-preserving for every legacy routine)" OFF)
option(NEXT_RDP_ENABLE_MANUAL_THREADS
    "Enable deterministic per-triplet pthread fan-out when OpenMP is unavailable" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(rdp_dna5_exact STATIC
    src/msvcrt_compat.cpp
    src/rdp_nj_legacy.cpp
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

# The original DNA5 DLL parallelizes AlistGC2/AlistMC3/AlistChi over
# independent triplet rows.  Use the platform runtime when requested; the
# compatibility omp.h remains serial for builds without a usable runtime.
if(NEXT_RDP_ENABLE_OPENMP AND NOT EMSCRIPTEN)
    find_package(OpenMP REQUIRED)
    target_compile_definitions(rdp_dna5_exact PRIVATE NEXT_RDP_USE_REAL_OPENMP=1)
    target_link_libraries(rdp_dna5_exact PUBLIC OpenMP::OpenMP_CXX)
endif()

# The original RDP tree path was built by the 32-bit MinGW source build with
# optimized x87 evaluation. Keep this isolated to the literal Clearcut/NJ
# translation unit; applying it to the complete DNA5 DLL changes upstream
# matrix generation.
set_source_files_properties(src/rdp_nj_legacy.cpp PROPERTIES
    COMPILE_OPTIONS "-O2"
)
if(NOT EMSCRIPTEN)
    set_property(SOURCE src/rdp_nj_legacy.cpp APPEND PROPERTY COMPILE_OPTIONS
        "-mfpmath=387" "-mno-sse2"
    )
endif()

add_library(next_rdp_core STATIC
    src/analysis.cpp
    src/burt_state.cpp
    src/distance_state.cpp
    src/identification_state.cpp
    src/legacy_method_state.cpp
    src/legacy_optional_bridge.cpp
    src/legacy_optional/bootscan.cpp
    src/legacy_optional/siscan.cpp
    src/mutation_state.cpp
    src/rescan_schedule.cpp
    src/round_state.cpp
    src/scan_state.cpp
    src/selection_state.cpp
    src/tree_state.cpp
    src/three_seq_state.cpp
    src/xover_state.cpp
)
target_include_directories(next_rdp_core
    PUBLIC src
    PRIVATE include vendor/dna5
)
target_link_libraries(next_rdp_core PUBLIC rdp_dna5_exact)
if(NEXT_RDP_ENABLE_OPENMP AND NOT EMSCRIPTEN)
    target_compile_definitions(next_rdp_core PRIVATE NEXT_RDP_USE_REAL_OPENMP=1)
endif()
if(NEXT_RDP_ENABLE_MANUAL_THREADS)
    find_package(Threads REQUIRED)
    target_compile_definitions(next_rdp_core PRIVATE NEXT_RDP_METHOD_MANUAL_THREADS=1)
    target_link_libraries(next_rdp_core PUBLIC Threads::Threads)
endif()

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
    if(NEXT_RDP_ENABLE_MANUAL_THREADS)
        target_compile_options(rdp_dna5_exact PRIVATE "-pthread")
        target_compile_options(next_rdp_core PRIVATE "-pthread")
        target_compile_options(rdp-core PRIVATE "-pthread")
    endif()
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
        $<$<BOOL:${NEXT_RDP_ENABLE_MANUAL_THREADS}>:-pthread>
        $<$<BOOL:${NEXT_RDP_ENABLE_MANUAL_THREADS}>:-sUSE_PTHREADS=1>
        $<$<BOOL:${NEXT_RDP_ENABLE_MANUAL_THREADS}>:-sPTHREAD_POOL_SIZE=8>
        "-sENVIRONMENT=node"
        "-sALLOW_MEMORY_GROWTH=1"
        "-sEXIT_RUNTIME=1"
        "-sNODERAWFS=1"
    )

    add_executable(next-rdp-core-web src/web_api.cpp)
    target_link_libraries(next-rdp-core-web PRIVATE next_rdp_core)
    target_compile_options(next-rdp-core-web PRIVATE
        "-ffp-eval-method=extended"
        $<$<BOOL:${NEXT_RDP_ENABLE_MANUAL_THREADS}>:-pthread>
    )
    set_target_properties(next-rdp-core-web PROPERTIES SUFFIX ".mjs")
    target_link_options(next-rdp-core-web PRIVATE
        $<$<BOOL:${NEXT_RDP_ENABLE_MANUAL_THREADS}>:-pthread>
        $<$<BOOL:${NEXT_RDP_ENABLE_MANUAL_THREADS}>:-sUSE_PTHREADS=1>
        $<$<BOOL:${NEXT_RDP_ENABLE_MANUAL_THREADS}>:-sPTHREAD_POOL_SIZE=8>
        "--no-entry"
        "-sENVIRONMENT=web,worker,node"
        "-sMODULARIZE=1"
        "-sEXPORT_ES6=1"
        "-sEXPORT_NAME=createNextRdpCore"
        "-sALLOW_MEMORY_GROWTH=1"
        "-sSTACK_SIZE=67108864"
        "-sINITIAL_MEMORY=134217728"
        "-sMAXIMUM_MEMORY=4294967296"
        "-sEXPORTED_FUNCTIONS=['_malloc','_free','_next_rdp_analyze','_next_rdp_result_json','_next_rdp_error','_next_rdp_version','_rdp_create','_rdp_destroy','_rdp_version','_rdp_set_worker_threads','_rdp_load_alignment','_rdp_get_summary_json','_rdp_scan_begin','_rdp_scan_batch','_rdp_reconcile','_rdp_cancel','_rdp_get_progress_json','_rdp_get_results_json','_rdp_get_signal_plot_json','_rdp_get_event_alignment_json','_rdp_get_event_trees_json','_rdp_get_event_phylpro_json','_rdp_set_review_state','_rdp_set_event_review_state','_rdp_update_event','_rdp_update_event_group','_rdp_reconcile_after','_rdp_export_csv','_rdp_export_enabled_sequences_fasta','_rdp_export_masked_or_disabled_sequences_fasta','_rdp_export_recombinant_sequences_removed_fasta','_rdp_export_recombinant_columns_removed_fasta','_rdp_export_recombination_free_fasta','_rdp_export_fragmented_fasta','_rdp_export_project_json','_rdp_get_error']"
        "-sEXPORTED_RUNTIME_METHODS=['UTF8ToString','HEAPU8']"
    )
endif()
