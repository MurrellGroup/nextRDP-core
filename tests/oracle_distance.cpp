#include "distance_fixture.hpp"

#include <windows.h>

#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: oracle-distance <DNA5.dll>\n";
        return 2;
    }

    HMODULE library = LoadLibraryA(argv[1]);
    if (library == nullptr) {
        std::cerr << "could not load DNA5.dll\n";
        return 1;
    }

    const auto super_dist_p = reinterpret_cast<SuperDistPFunction>(
        GetProcAddress(library, "SuperDistP"));
    if (super_dist_p == nullptr) {
        std::cerr << "SuperDistP export was not found\n";
        FreeLibrary(library);
        return 1;
    }

    const int result = run_distance_fixture(super_dist_p, std::cout);
    FreeLibrary(library);
    return result;
}

