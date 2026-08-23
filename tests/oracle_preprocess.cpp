#include "preprocess_fixture.hpp"

#include <windows.h>

#include <iostream>

template <typename Function>
Function load_function(HMODULE library, const char* name) {
    return reinterpret_cast<Function>(GetProcAddress(library, name));
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: oracle-preprocess <DNA5.dll>\n";
        return 2;
    }

    HMODULE library = LoadLibraryA(argv[1]);
    if (library == nullptr) {
        std::cerr << "could not load DNA5.dll\n";
        return 1;
    }

    PreprocessApi api{
        load_function<decltype(PreprocessApi::make_a_list_p2)>(library, "MakeAListP2"),
        load_function<decltype(PreprocessApi::count_nucs)>(library, "CountNucs"),
        load_function<decltype(PreprocessApi::recode_nucs)>(library, "RecodeNucs"),
        load_function<decltype(PreprocessApi::do_recode_p)>(library, "DoRecodeP"),
        load_function<decltype(PreprocessApi::make_compress_seq_p)>(library, "MakeCompressSeqP"),
    };

    if (api.make_a_list_p2 == nullptr || api.count_nucs == nullptr ||
        api.recode_nucs == nullptr || api.do_recode_p == nullptr ||
        api.make_compress_seq_p == nullptr) {
        std::cerr << "one or more DNA5 exports were not found\n";
        FreeLibrary(library);
        return 1;
    }

    const int result = run_preprocess_fixture(api, std::cout);
    FreeLibrary(library);
    return result;
}

