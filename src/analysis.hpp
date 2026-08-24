#pragma once

#include "scan_state.hpp"
#include "xover_state.hpp"

#include <string>

struct RdpInitialAnalysisOptions {
    bool circular = true;
    double p_value_cutoff = 0.05;
    int window_sites = 30;
};

struct RdpInitialAnalysisResult {
    RdpScanState alignment;
    RdpRawEventState events;
};

RdpInitialAnalysisResult run_rdp_initial_analysis_from_fasta_text(
    const std::string& fasta_text,
    const RdpInitialAnalysisOptions& options = {});

RdpInitialAnalysisResult run_rdp_initial_analysis_from_fasta_file(
    const std::string& fasta_path,
    const RdpInitialAnalysisOptions& options = {});
