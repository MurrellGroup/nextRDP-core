#pragma once

#include <cstdint>
#include <string>
#include <vector>

#if defined(_WIN32)
#define RDP_SCAN_CALL __stdcall
#else
#define RDP_SCAN_CALL
#endif

struct Dna5ScanPreprocessApi {
    int(RDP_SCAN_CALL* make_analysis_list)(float, int, short*, int, short*);
    int(RDP_SCAN_CALL* count_nucleotides)(
        int, int, int, short*, unsigned char*, int, int*);
    int(RDP_SCAN_CALL* rank_nucleotides)(
        int, int, int, int*, int, unsigned char*);
    int(RDP_SCAN_CALL* apply_recoding)(
        int, int, int, short*, int, unsigned char*, unsigned char*, int,
        unsigned char*);
    int(RDP_SCAN_CALL* compress_sequences)(
        int, int, unsigned char*, int, unsigned char*, int, int,
        unsigned char*);
};

struct RdpScanState {
    int next_no = -1;
    int sequence_length = 0;
    int analysis_list_last = -1;
    int compressed_sequence_ub = -1;
    std::vector<short> sequence_data;
    std::vector<short> analysis_list;
    std::vector<unsigned char> compressed_sequence;
};

RdpScanState build_rdp_scan_state_from_fasta(
    const std::string& fasta_path, const Dna5ScanPreprocessApi& api);

