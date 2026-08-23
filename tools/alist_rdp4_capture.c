/*
 * One-shot Windows oracle capture for DNA5!AlistRDP4.
 *
 * This DLL replaces DNA5.dll next to a sandboxed RDP5CL.exe.  Every export
 * except AlistRDP4 is forwarded by the generated .def file.  AlistRDP4 is
 * captured before and after the call and delegated to DNA5_original.dll.
 * There is deliberately no C runtime I/O here: the proxy uses kernel32 so it
 * can be built as a small, locally loadable 32-bit DLL.
 */

typedef void *HANDLE;
typedef void *HMODULE;
typedef void *FARPROC;
typedef unsigned long DWORD;
typedef unsigned int u32;
typedef int BOOL;

#define STDCALL __attribute__((stdcall))
#define DLLIMPORT __declspec(dllimport)
#define GENERIC_WRITE 0x40000000UL
#define FILE_SHARE_READ 0x00000001UL
#define FILE_SHARE_WRITE 0x00000002UL
#define CREATE_ALWAYS 2UL
#define OPEN_EXISTING 3UL
#define FILE_ATTRIBUTE_NORMAL 0x00000080UL
#define FILE_END 2UL
#define INVALID_HANDLE_VALUE ((HANDLE)(long)-1)

DLLIMPORT HMODULE STDCALL LoadLibraryA(const char *name);
DLLIMPORT FARPROC STDCALL GetProcAddress(HMODULE module, const char *name);
DLLIMPORT HANDLE STDCALL CreateFileA(
    const char *name, DWORD access, DWORD share, void *security,
    DWORD creation, DWORD attributes, HANDLE template_file);
DLLIMPORT DWORD STDCALL SetFilePointer(
    HANDLE file, long distance, long *distance_high, DWORD method);
DLLIMPORT BOOL STDCALL WriteFile(
    HANDLE file, const void *buffer, DWORD bytes, DWORD *written,
    void *overlapped);
DLLIMPORT BOOL STDCALL CloseHandle(HANDLE handle);

typedef int(STDCALL *AlistRDP4Fn)(
    int, double *, short *, int, int, int, int, double, unsigned char *, int,
    int, int, double, int, int, int, int, float *, int, float *, int, int,
    unsigned char *, short *, int, short, unsigned char *, int, int, int,
    double *, int, double *, double *);
typedef int(STDCALL *FindSubSeqPB3Fn)(
    int *, int, int, int, int, int, int, int, int, unsigned char *, int,
    char *, unsigned char *);
typedef int(STDCALL *XOHomologyPFn)(short, int, int, short, char *, int *);
typedef int(STDCALL *FindNextPFn)(int, int, int, int, int, int, int, int *);
typedef int(STDCALL *DefineEventP2Fn)(
    int, int, int, int, int, int, int, int, int, int, int, int, int, int,
    int *, int *, int *, int *, int *, char *, int *);
typedef double(STDCALL *ProbCalcP2Fn)(double *, int, int, int, double, int);
typedef int(STDCALL *FindSubSeqPB4Fn)(
    int *, int, int, int, int, int, int, int, int, unsigned char *, int,
    char *, int *, int *, unsigned char *);
typedef int(STDCALL *FindFirstCOPFn)(int, int, int, int, int, int *);
typedef double(STDCALL *ProbCalcPFn)(double *, int, int, double, int);
typedef int(STDCALL *CleanXOSNWFn)(int, int, int, char *);

typedef struct XOVERDEFINE {
  unsigned char OutsideFlag;
  unsigned char MissIdentifyFlag;
  unsigned char ProgramFlag;
  unsigned char SBPFlag;
  unsigned char Accept;
  short MajorP;
  short MinorP;
  short Daughter;
  int Beginning;
  int Ending;
  int LHolder;
  int Eventnumber;
  float PermPVal;
  int BeginP;
  int EndP;
  double Probability;
  double DHolder;
} XOVERDEFINE;

typedef int(STDCALL *MakeTestPVsFn)(
    int, unsigned char *, int, int, int, short *, XOVERDEFINE *, double *);
typedef int(STDCALL *FindBestRecSignalP2Fn)(
    char, int, int, int, double *, char *, int *, short *, double *);
typedef int(STDCALL *UFDistFn)(
    int, int, int, int, float *, float *, float *, float *, int *, int,
    short *);
typedef double(STDCALL *SuperDistP2Fn)(
    int, int, int, int, int, int, int, int, int, double *, float *, float *,
    float *, short *, int *, short *, short *, short *, short *, short *,
    short *, short *, char *, char *, char *, char *, char *, char *, char *,
    char *, char *, char *, char *);
typedef SuperDistP2Fn SuperDistPFn;
typedef int(STDCALL *CheckMatrixPFn)(
    int *, int *, int, int, int, int, unsigned char *, int, float *, int,
    float *, int, float *, float *, int *, int *);
typedef int(STDCALL *MakeNJTreesP2Fn)(
    int, int, int, int *, unsigned char *, unsigned char *, int, int, int,
    int, int *, int *, int, float *, int, float *, int, float *, int, float *,
    int *, char *, char *, float *, float *);
typedef int(STDCALL *MarkOutsidesFn)(
    int, unsigned char *, int, int, short *, XOVERDEFINE *);
typedef int(STDCALL *MakeSDMP2Fn)(
    int, int, int *, int *, int *, int *, unsigned char *, short *, double *,
    double *);
typedef int(STDCALL *FillRmatFn)(
    int, int, int, int, int, int, int, double *, double *, unsigned char *);

#pragma pack(push, 1)
struct CaptureHeader {
  char magic[8];
  u32 version;
  int store_lpv_ub;
  int list_length;
  int start;
  int end;
  int next_no;
  double sub_threshold;
  int circular;
  int mc_correction;
  int mc_flag;
  double lowest_probability;
  int target_x;
  int sequence_length;
  int short_output;
  int distance_ub;
  int tree_distance_ub;
  int fss_rdp_ub;
  int compressed_sequence_ub;
  int xover_window;
  short xover_window_x;
  int probability_file_flag;
  int probability_one_ub;
  int probability_two_ub;
  int fact_three_ub;
};

struct SectionHeader {
  u32 id;
  u32 bytes;
};

struct FindSubSeqPB3Header {
  char magic[8];
  u32 version;
  int fss_ub;
  int xover_window;
  int compressed_sequence_ub;
  int sequence_length;
  int next_no;
  int seq1;
  int seq2;
  int seq3;
  int xover_sequence_ub;
};

struct XOHomologyPHeader {
  char magic[8];
  u32 version;
  short inlyer;
  int sequence_length;
  int xover_length;
  short xover_window;
};

struct FindNextPHeader {
  char magic[8];
  u32 version;
  int homology_ub;
  int start;
  int high;
  int med;
  int low;
  int xover_length;
  int xover_window;
};

struct DefineEventP2Header {
  char magic[8];
  u32 version;
  int homology_ub;
  int short_output;
  int long_winded;
  int med;
  int high;
  int low;
  int target_x;
  int circular;
  int xx;
  int xover_window;
  int sequence_length;
  int xover_length;
  int sequence_daughter;
  int sequence_minor;
};

struct ProbCalcP2Header {
  char magic[8];
  u32 version;
  int fact_three_ub;
  int xover_length;
  int number_in_common;
  double individual_probability;
  int informative_length;
};

struct FindSubSeqPB4Header {
  char magic[8];
  u32 version;
  int fss_ub;
  int xover_window;
  int compressed_sequence_ub;
  int sequence_length;
  int next_no;
  int seq1;
  int seq2;
  int seq3;
  int xover_sequence_ub;
};

struct FindFirstCOPHeader {
  char magic[8];
  u32 version;
  int x;
  int med;
  int high;
  int xover_length;
  int homology_ub;
};

struct ProbCalcPHeader {
  char magic[8];
  u32 version;
  int xover_length;
  int number_in_common;
  double individual_probability;
  int informative_length;
};

struct CleanXOSNWHeader {
  char magic[8];
  u32 version;
  int xover_length;
  int xover_window;
  int xover_sequence_ub;
};

struct MakeTestPVsHeader {
  char magic[8];
  u32 version;
  int done_sequence_ub;
  int next_no;
  int xover_rows_ub;
  int xover_slots_ub;
  int xover_struct_bytes;
};

struct FindBestRecSignalP2Header {
  char magic[8];
  u32 version;
  char done_target;
  int next_no;
  int probability_rows_ub;
  int probability_columns_ub;
};

struct UFDistHeader {
  char magic[8];
  u32 version;
  int sequence_length;
  int begin;
  int end;
  int pair_matrix_ub;
  int sequence_data_ub;
};

struct SuperDistP2Header {
  char magic[8];
  u32 version;
  int x;
  int next_no;
  int ub14;
  int ub04;
  int ub13;
  int ub03;
  int ub12;
  int ub02;
  int ub11;
};

struct CheckMatrixPHeader {
  char magic[8];
  u32 version;
  int next_no;
  int sco;
  int minimum_sequence_size;
  int missing_pair_ub;
  int valid_ub;
  int sub_valid_ub;
  int matrix_ub;
};

struct MakeNJTreesP2Header {
  char magic[8];
  u32 version;
  int resolve_root;
  int nseqs;
  int next_no;
  int seed;
  int name_length;
  int sequence_length;
  int trace_sequences_ub;
  int first_matrix_ub;
  int second_matrix_ub;
  int first_adjusted_matrix_ub;
  int second_adjusted_matrix_ub;
};

struct MarkOutsidesHeader {
  char magic[8];
  u32 version;
  int done_sequence_ub;
  int next_no;
  int xover_rows_ub;
  int maximum_current_xover;
  int xover_struct_bytes;
};

struct MakeSDMP2Header {
  char magic[8];
  u32 version;
  int next_no;
  int sequence_length;
};

struct FillRmatHeader {
  char magic[8];
  u32 version;
  int y;
  int next_no;
  int result_matrix_ub1;
  int result_matrix_ub2;
  int distance_matrix_ub1;
  int distance_matrix_ub2;
  int distance_matrix_ub3;
};
#pragma pack(pop)

enum SectionId {
  STORE_LPV_IN = 1,
  ANALYSIS_LIST_IN = 2,
  REDO_LIST_IN = 3,
  DISTANCE_IN = 4,
  TREE_DISTANCE_IN = 5,
  FSS_RDP_IN = 6,
  COMPRESSED_SEQUENCE_IN = 7,
  SEQUENCE_DATA_IN = 8,
  PROBABILITY_ESTIMATE_IN = 9,
  FACT_THREE_IN = 10,
  FACT_IN = 11,
  REDO_LIST_OUT = 101,
  STORE_LPV_OUT = 102,
  RESULT_OUT = 103,
  END_MARKER = 0xffffffffU
};

enum FindSubSeqPB3SectionId {
  FSP_AH_IN = 1,
  FSP_COMPRESSED_SEQUENCE_IN = 2,
  FSP_XOVER_SEQUENCE_IN = 3,
  FSP_FSS_RDP_IN = 4,
  FSP_AH_OUT = 101,
  FSP_XOVER_SEQUENCE_OUT = 102,
  FSP_RESULT_OUT = 103
};

enum XOHomologyPSectionId {
  XOH_XOVER_SEQUENCE_IN = 1,
  XOH_HOMOLOGY_IN = 2,
  XOH_HOMOLOGY_OUT = 101,
  XOH_RESULT_OUT = 102
};

enum FindNextPSectionId {
  FNP_HOMOLOGY_IN = 1,
  FNP_RESULT_OUT = 101
};

enum DefineEventP2SectionId {
  DEP_SCALARS_IN = 1,
  DEP_XOVER_SEQUENCE_IN = 2,
  DEP_HOMOLOGY_IN = 3,
  DEP_SCALARS_OUT = 101,
  DEP_RESULT_OUT = 102
};

enum ProbCalcP2SectionId {
  PCP2_FACT_THREE_IN = 1,
  PCP2_RESULT_OUT = 101
};

enum FindSubSeqPB4SectionId {
  FSP4_AH_IN = 1,
  FSP4_COMPRESSED_SEQUENCE_IN = 2,
  FSP4_XOVER_SEQUENCE_IN = 3,
  FSP4_XDIFFPOS_IN = 4,
  FSP4_XPOSDIFF_IN = 5,
  FSP4_FSS_RDP_IN = 6,
  FSP4_AH_OUT = 101,
  FSP4_XOVER_SEQUENCE_OUT = 102,
  FSP4_XDIFFPOS_OUT = 103,
  FSP4_XPOSDIFF_OUT = 104,
  FSP4_RESULT_OUT = 105
};

enum FindFirstCOPSectionId {
  FFCO_HOMOLOGY_IN = 1,
  FFCO_RESULT_OUT = 101
};

enum ProbCalcPSectionId {
  PCP_FACT_IN = 1,
  PCP_RESULT_OUT = 101
};

enum CleanXOSNWSectionId {
  CXO_XOVER_SEQUENCE_IN = 1,
  CXO_XOVER_SEQUENCE_OUT = 101,
  CXO_RESULT_OUT = 102
};

enum MakeTestPVsSectionId {
  MTP_DONE_SEQUENCE_IN = 1,
  MTP_CURRENT_XOVER_IN = 2,
  MTP_XOVER_LIST_IN = 3,
  MTP_TEST_PVS_IN = 4,
  MTP_DONE_SEQUENCE_OUT = 101,
  MTP_CURRENT_XOVER_OUT = 102,
  MTP_XOVER_LIST_OUT = 103,
  MTP_TEST_PVS_OUT = 104,
  MTP_RESULT_OUT = 105
};

enum FindBestRecSignalP2SectionId {
  FBR_LOW_P_IN = 1,
  FBR_DONE_SEQUENCE_IN = 2,
  FBR_TRACE_IN = 3,
  FBR_CURRENT_XOVER_IN = 4,
  FBR_TEST_PVS_IN = 5,
  FBR_LOW_P_OUT = 101,
  FBR_DONE_SEQUENCE_OUT = 102,
  FBR_TRACE_OUT = 103,
  FBR_CURRENT_XOVER_OUT = 104,
  FBR_TEST_PVS_OUT = 105,
  FBR_RESULT_OUT = 106
};

enum UFDistSectionId {
  UFD_VALID_IN = 1,
  UFD_DIFFERENCES_IN = 2,
  UFD_BREAKPOINT_DISTANCE_IN = 3,
  UFD_REMAINDER_DISTANCE_IN = 4,
  UFD_SEQUENCES_IN = 5,
  UFD_SEQUENCE_DATA_IN = 6,
  UFD_BREAKPOINT_DISTANCE_OUT = 101,
  UFD_REMAINDER_DISTANCE_OUT = 102,
  UFD_RESULT_OUT = 103
};

enum SuperDistP2SectionId {
  SDP2_AVERAGE_IN = 1,
  SDP2_PAIR_DIFF_IN = 2,
  SDP2_PAIR_VALID_IN = 3,
  SDP2_DISTANCE_IN = 4,
  SDP2_REDO_IN = 5,
  SDP2_CATEGORY_COUNT_IN = 6,
  SDP2_ISEQ14_IN = 7,
  SDP2_ISEQ04_IN = 8,
  SDP2_ISEQ13_IN = 9,
  SDP2_ISEQ03_IN = 10,
  SDP2_ISEQ12_IN = 11,
  SDP2_ISEQ02_IN = 12,
  SDP2_ISEQ11_IN = 13,
  SDP2_VALID14_IN = 14,
  SDP2_DIFF14_IN = 15,
  SDP2_VALID13_IN = 16,
  SDP2_DIFF13_IN = 17,
  SDP2_VALID12_IN = 18,
  SDP2_DIFF12_IN = 19,
  SDP2_VALID11_IN = 20,
  SDP2_DIFF11_IN = 21,
  SDP2_DIFF04_IN = 22,
  SDP2_DIFF03_IN = 23,
  SDP2_DIFF02_IN = 24,
  SDP2_AVERAGE_OUT = 101,
  SDP2_PAIR_DIFF_OUT = 102,
  SDP2_PAIR_VALID_OUT = 103,
  SDP2_DISTANCE_OUT = 104,
  SDP2_RESULT_OUT = 105
};

enum CheckMatrixPSectionId {
  CMP_MINIMUMS_IN = 1,
  CMP_SEQUENCES_IN = 2,
  CMP_MISSING_PAIR_IN = 3,
  CMP_VALID_IN = 4,
  CMP_SUB_VALID_IN = 5,
  CMP_FIRST_MATRIX_IN = 6,
  CMP_SECOND_MATRIX_IN = 7,
  CMP_FIRST_TOTAL_IN = 8,
  CMP_SECOND_TOTAL_IN = 9,
  CMP_MINIMUMS_OUT = 101,
  CMP_MISSING_PAIR_OUT = 102,
  CMP_FIRST_MATRIX_OUT = 103,
  CMP_SECOND_MATRIX_OUT = 104,
  CMP_FIRST_TOTAL_OUT = 105,
  CMP_SECOND_TOTAL_OUT = 106,
  CMP_RESULT_OUT = 107
};

enum MakeNJTreesP2SectionId {
  MNJ_SEQUENCES_IN = 1,
  MNJ_MIN_PAIR_IN = 2,
  MNJ_SEQUENCE_PAIR_IN = 3,
  MNJ_OUTLIER_IN = 4,
  MNJ_TRACE_SEQUENCES_IN = 5,
  MNJ_FIRST_MATRIX_IN = 6,
  MNJ_SECOND_MATRIX_IN = 7,
  MNJ_FIRST_ADJUSTED_MATRIX_IN = 8,
  MNJ_SECOND_ADJUSTED_MATRIX_IN = 9,
  MNJ_REDO_LIST_IN = 10,
  MNJ_FIRST_HOLDER_IN = 11,
  MNJ_SECOND_HOLDER_IN = 12,
  MNJ_TEMP_FIRST_MATRIX_IN = 13,
  MNJ_TEMP_SECOND_MATRIX_IN = 14,
  MNJ_SEQUENCES_OUT = 101,
  MNJ_MIN_PAIR_OUT = 102,
  MNJ_SEQUENCE_PAIR_OUT = 103,
  MNJ_OUTLIER_OUT = 104,
  MNJ_TRACE_SEQUENCES_OUT = 105,
  MNJ_FIRST_MATRIX_OUT = 106,
  MNJ_SECOND_MATRIX_OUT = 107,
  MNJ_FIRST_ADJUSTED_MATRIX_OUT = 108,
  MNJ_SECOND_ADJUSTED_MATRIX_OUT = 109,
  MNJ_REDO_LIST_OUT = 110,
  MNJ_FIRST_HOLDER_OUT = 111,
  MNJ_SECOND_HOLDER_OUT = 112,
  MNJ_TEMP_FIRST_MATRIX_OUT = 113,
  MNJ_TEMP_SECOND_MATRIX_OUT = 114,
  MNJ_RESULT_OUT = 115
};

enum MarkOutsidesSectionId {
  MO_DONE_SEQUENCE_IN = 1,
  MO_CURRENT_XOVER_IN = 2,
  MO_XOVER_LIST_IN = 3,
  MO_DONE_SEQUENCE_OUT = 101,
  MO_CURRENT_XOVER_OUT = 102,
  MO_XOVER_LIST_OUT = 103,
  MO_RESULT_OUT = 104
};

enum MakeSDMP2SectionId {
  MSD_START_POSITIONS_IN = 1,
  MSD_END_POSITIONS_IN = 2,
  MSD_SEQUENCES_IN = 3,
  MSD_COMPARISON_MATRIX_IN = 4,
  MSD_MISSING_DATA_IN = 5,
  MSD_SEQUENCE_DATA_IN = 6,
  MSD_SUMMARY_MATRIX_IN = 7,
  MSD_DISTANCE_MATRIX_IN = 8,
  MSD_SUMMARY_MATRIX_OUT = 101,
  MSD_DISTANCE_MATRIX_OUT = 102,
  MSD_RESULT_OUT = 103
};

enum FillRmatSectionId {
  FRM_RESULT_MATRIX_IN = 1,
  FRM_DISTANCE_MATRIX_IN = 2,
  FRM_POSITIONS_IN = 3,
  FRM_RESULT_MATRIX_OUT = 101,
  FRM_DISTANCE_MATRIX_OUT = 102,
  FRM_POSITIONS_OUT = 103,
  FRM_RESULT_OUT = 104
};

static void write_bytes(HANDLE file, const void *buffer, DWORD bytes) {
  DWORD written = 0;
  WriteFile(file, buffer, bytes, &written, (void *)0);
}

static void write_section(HANDLE file, u32 id, const void *data, u32 bytes) {
  struct SectionHeader section;
  section.id = id;
  section.bytes = bytes;
  write_bytes(file, &section, (DWORD)sizeof(section));
  if (bytes != 0) write_bytes(file, data, bytes);
}

int STDCALL XOHomologyPCapture(
    short inlyer, int sequence_length, int xover_length, short xover_window,
    char *xover_sequence, int *homology) {
  static XOHomologyPFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const int stride = sequence_length + xover_window * 2;

  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (XOHomologyPFn)GetProcAddress(module, "XOHomologyP");
  }
  if (!original) return 0;
  ++invocation;

  if (invocation == 1) {
    struct XOHomologyPHeader header;
    const char magic[8] = {'X', 'O', 'H', 'O', 'M', 'P', '\0', '\0'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.inlyer = inlyer;
    header.sequence_length = sequence_length;
    header.xover_length = xover_length;
    header.xover_window = xover_window;
    file = CreateFileA(
        "xohomology-p-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, XOH_XOVER_SEQUENCE_IN, xover_sequence, (u32)(3 * stride));
      write_section(
          file, XOH_HOMOLOGY_IN, homology,
          (u32)(3 * stride * (int)sizeof(int)));
      CloseHandle(file);
    }
  }

  {
    const int result = original(
        inlyer, sequence_length, xover_length, xover_window, xover_sequence,
        homology);
    if (invocation == 1) {
      file = CreateFileA(
          "xohomology-p-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(
            file, XOH_HOMOLOGY_OUT, homology,
            (u32)(3 * stride * (int)sizeof(int)));
        write_section(file, XOH_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end, (DWORD)sizeof(end));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL FindNextPCapture(
    int homology_ub, int start, int high, int med, int low, int xover_length,
    int xover_window, int *homology) {
  static FindNextPFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (FindNextPFn)GetProcAddress(module, "FindNextP");
  }
  if (!original) return -1;
  ++invocation;
  if (invocation == 1) {
    struct FindNextPHeader header;
    const char magic[8] = {'F', 'I', 'N', 'D', 'N', 'X', 'T', '\0'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.homology_ub = homology_ub;
    header.start = start;
    header.high = high;
    header.med = med;
    header.low = low;
    header.xover_length = xover_length;
    header.xover_window = xover_window;
    file = CreateFileA(
        "find-next-p-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(
          file, FNP_HOMOLOGY_IN, homology,
          (u32)(3 * (homology_ub + 1) * (int)sizeof(int)));
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        homology_ub, start, high, med, low, xover_length, xover_window,
        homology);
    if (invocation == 1) {
      file = CreateFileA(
          "find-next-p-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, FNP_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end, (DWORD)sizeof(end));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL DefineEventP2Capture(
    int homology_ub, int short_output, int long_winded, int med, int high,
    int low, int target_x, int circular, int xx, int xover_window,
    int sequence_length, int xover_length, int sequence_daughter,
    int sequence_minor, int *end_flag, int *begin, int *end, int *ncommon,
    int *event_length, char *xover_sequence, int *homology) {
  static DefineEventP2Fn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  int scalar_state[5];

  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (DefineEventP2Fn)GetProcAddress(module, "DefineEventP2");
  }
  if (!original) return -1;
  ++invocation;
  scalar_state[0] = *end_flag;
  scalar_state[1] = *begin;
  scalar_state[2] = *end;
  scalar_state[3] = *ncommon;
  scalar_state[4] = *event_length;

  if (invocation == 1) {
    struct DefineEventP2Header header;
    const char magic[8] = {'D', 'E', 'F', 'E', 'V', 'P', '2', '\0'};
    const int xover_stride = sequence_length + 1 + xover_window * 2;
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.homology_ub = homology_ub;
    header.short_output = short_output;
    header.long_winded = long_winded;
    header.med = med;
    header.high = high;
    header.low = low;
    header.target_x = target_x;
    header.circular = circular;
    header.xx = xx;
    header.xover_window = xover_window;
    header.sequence_length = sequence_length;
    header.xover_length = xover_length;
    header.sequence_daughter = sequence_daughter;
    header.sequence_minor = sequence_minor;
    file = CreateFileA(
        "define-event-p2-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, DEP_SCALARS_IN, scalar_state, (u32)sizeof(scalar_state));
      write_section(
          file, DEP_XOVER_SEQUENCE_IN, xover_sequence,
          (u32)(3 * xover_stride));
      write_section(
          file, DEP_HOMOLOGY_IN, homology,
          (u32)(3 * (homology_ub + 1) * (int)sizeof(int)));
      CloseHandle(file);
    }
  }

  {
    const int result = original(
        homology_ub, short_output, long_winded, med, high, low, target_x,
        circular, xx, xover_window, sequence_length, xover_length,
        sequence_daughter, sequence_minor, end_flag, begin, end, ncommon,
        event_length, xover_sequence, homology);
    if (invocation == 1) {
      scalar_state[0] = *end_flag;
      scalar_state[1] = *begin;
      scalar_state[2] = *end;
      scalar_state[3] = *ncommon;
      scalar_state[4] = *event_length;
      file = CreateFileA(
          "define-event-p2-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(
            file, DEP_SCALARS_OUT, scalar_state, (u32)sizeof(scalar_state));
        write_section(file, DEP_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

double STDCALL ProbCalcP2Capture(
    double *fact_three, int fact_three_ub, int xover_length,
    int number_in_common, double individual_probability,
    int informative_length) {
  static ProbCalcP2Fn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (ProbCalcP2Fn)GetProcAddress(module, "ProbCalcP2");
  }
  if (!original) return 0.0;
  ++invocation;
  if (invocation == 1) {
    struct ProbCalcP2Header header;
    const char magic[8] = {'P', 'R', 'O', 'B', 'C', 'P', '2', '\0'};
    const int width = fact_three_ub + 1;
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.fact_three_ub = fact_three_ub;
    header.xover_length = xover_length;
    header.number_in_common = number_in_common;
    header.individual_probability = individual_probability;
    header.informative_length = informative_length;
    file = CreateFileA(
        "prob-calc-p2-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(
          file, PCP2_FACT_THREE_IN, fact_three,
          (u32)(width * width * width * (int)sizeof(double)));
      CloseHandle(file);
    }
  }
  {
    const double result = original(
        fact_three, fact_three_ub, xover_length, number_in_common,
        individual_probability, informative_length);
    if (invocation == 1) {
      file = CreateFileA(
          "prob-calc-p2-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, PCP2_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL FindSubSeqPB4Capture(
    int *ah, int fss_ub, int xover_window, int compressed_sequence_ub,
    int sequence_length, int next_no, int seq1, int seq2, int seq3,
    unsigned char *compressed_sequence, int xover_sequence_ub,
    char *xover_sequence, int *xdiffpos, int *xposdiff,
    unsigned char *fss_rdp) {
  static FindSubSeqPB4Fn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 position_bytes =
      (u32)((sequence_length + 201) * (int)sizeof(int));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (FindSubSeqPB4Fn)GetProcAddress(module, "FindSubSeqPB4");
  }
  if (!original) return 0;
  ++invocation;
  if (invocation == 1) {
    struct FindSubSeqPB4Header header;
    const char magic[8] = {'F', 'S', 'P', 'B', '4', '\0', '\0', '\0'};
    const int fss_width = fss_ub + 1;
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.fss_ub = fss_ub;
    header.xover_window = xover_window;
    header.compressed_sequence_ub = compressed_sequence_ub;
    header.sequence_length = sequence_length;
    header.next_no = next_no;
    header.seq1 = seq1;
    header.seq2 = seq2;
    header.seq3 = seq3;
    header.xover_sequence_ub = xover_sequence_ub;
    file = CreateFileA(
        "find-subseq-pb4-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, FSP4_AH_IN, ah, (u32)(4 * (int)sizeof(int)));
      write_section(
          file, FSP4_COMPRESSED_SEQUENCE_IN, compressed_sequence,
          (u32)((compressed_sequence_ub + 1) * (next_no + 1)));
      write_section(
          file, FSP4_XOVER_SEQUENCE_IN, xover_sequence,
          (u32)((xover_sequence_ub + 1) * 3));
      write_section(file, FSP4_XDIFFPOS_IN, xdiffpos, position_bytes);
      write_section(file, FSP4_XPOSDIFF_IN, xposdiff, position_bytes);
      write_section(
          file, FSP4_FSS_RDP_IN, fss_rdp,
          (u32)(4 * fss_width * fss_width * fss_width));
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        ah, fss_ub, xover_window, compressed_sequence_ub, sequence_length,
        next_no, seq1, seq2, seq3, compressed_sequence, xover_sequence_ub,
        xover_sequence, xdiffpos, xposdiff, fss_rdp);
    if (invocation == 1) {
      file = CreateFileA(
          "find-subseq-pb4-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, FSP4_AH_OUT, ah, (u32)(4 * (int)sizeof(int)));
        write_section(
            file, FSP4_XOVER_SEQUENCE_OUT, xover_sequence,
            (u32)((xover_sequence_ub + 1) * 3));
        write_section(file, FSP4_XDIFFPOS_OUT, xdiffpos, position_bytes);
        write_section(file, FSP4_XPOSDIFF_OUT, xposdiff, position_bytes);
        write_section(file, FSP4_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL FindFirstCOPCapture(
    int x, int med, int high, int xover_length, int homology_ub,
    int *homology) {
  static FindFirstCOPFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (FindFirstCOPFn)GetProcAddress(module, "FindFirstCOP");
  }
  if (!original) return -1;
  ++invocation;
  if (invocation == 1) {
    struct FindFirstCOPHeader header;
    const char magic[8] = {'F', 'F', 'I', 'R', 'S', 'T', '\0', '\0'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.x = x;
    header.med = med;
    header.high = high;
    header.xover_length = xover_length;
    header.homology_ub = homology_ub;
    file = CreateFileA(
        "find-first-co-p-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(
          file, FFCO_HOMOLOGY_IN, homology,
          (u32)(3 * (homology_ub + 1) * (int)sizeof(int)));
      CloseHandle(file);
    }
  }
  {
    const int result = original(x, med, high, xover_length, homology_ub, homology);
    if (invocation == 1) {
      file = CreateFileA(
          "find-first-co-p-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, FFCO_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

double STDCALL ProbCalcPCapture(
    double *fact, int xover_length, int number_in_common,
    double individual_probability, int informative_length) {
  static ProbCalcPFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (ProbCalcPFn)GetProcAddress(module, "ProbCalcP");
  }
  if (!original) return 0.0;
  ++invocation;
  if (invocation == 1) {
    struct ProbCalcPHeader header;
    const char magic[8] = {'P', 'R', 'O', 'B', 'C', 'P', '\0', '\0'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.xover_length = xover_length;
    header.number_in_common = number_in_common;
    header.individual_probability = individual_probability;
    header.informative_length = informative_length;
    file = CreateFileA(
        "prob-calc-p-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, PCP_FACT_IN, fact, (u32)(172 * (int)sizeof(double)));
      CloseHandle(file);
    }
  }
  {
    const double result = original(
        fact, xover_length, number_in_common, individual_probability,
        informative_length);
    if (invocation == 1) {
      file = CreateFileA(
          "prob-calc-p-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, PCP_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL CleanXOSNWCapture(
    int xover_length, int xover_window, int xover_sequence_ub,
    char *xover_sequence) {
  static CleanXOSNWFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 sequence_bytes = (u32)(3 * (xover_sequence_ub + 1));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (CleanXOSNWFn)GetProcAddress(module, "CleanXOSNW");
  }
  if (!original) return 0;
  ++invocation;
  if (invocation == 1) {
    struct CleanXOSNWHeader header;
    const char magic[8] = {'C', 'L', 'N', 'X', 'O', 'S', 'N', 'W'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.xover_length = xover_length;
    header.xover_window = xover_window;
    header.xover_sequence_ub = xover_sequence_ub;
    file = CreateFileA(
        "clean-xosnw-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, CXO_XOVER_SEQUENCE_IN, xover_sequence, sequence_bytes);
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        xover_length, xover_window, xover_sequence_ub, xover_sequence);
    if (invocation == 1) {
      file = CreateFileA(
          "clean-xosnw-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, CXO_XOVER_SEQUENCE_OUT, xover_sequence, sequence_bytes);
        write_section(file, CXO_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL MakeTestPVsCapture(
    int done_sequence_ub, unsigned char *done_sequence, int next_no,
    int xover_rows_ub, int xover_slots_ub, short *current_xover,
    XOVERDEFINE *xover_list, double *test_pvs) {
  static MakeTestPVsFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 done_bytes =
      (u32)((done_sequence_ub + 1) * (xover_slots_ub + 1));
  const u32 current_bytes = (u32)((next_no + 1) * (int)sizeof(short));
  const u32 list_bytes = (u32)((xover_rows_ub + 1) * (xover_slots_ub + 1) *
                               (int)sizeof(XOVERDEFINE));
  const u32 probability_bytes =
      (u32)((next_no + 1) * (xover_slots_ub + 1) * (int)sizeof(double));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (MakeTestPVsFn)GetProcAddress(module, "MakeTestPVs");
  }
  if (!original) return 0;
  ++invocation;
  if (invocation == 1) {
    struct MakeTestPVsHeader header;
    const char magic[8] = {'M', 'K', 'T', 'E', 'S', 'T', 'P', 'V'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.done_sequence_ub = done_sequence_ub;
    header.next_no = next_no;
    header.xover_rows_ub = xover_rows_ub;
    header.xover_slots_ub = xover_slots_ub;
    header.xover_struct_bytes = (int)sizeof(XOVERDEFINE);
    file = CreateFileA(
        "make-test-pvs-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, MTP_DONE_SEQUENCE_IN, done_sequence, done_bytes);
      write_section(file, MTP_CURRENT_XOVER_IN, current_xover, current_bytes);
      write_section(file, MTP_XOVER_LIST_IN, xover_list, list_bytes);
      write_section(file, MTP_TEST_PVS_IN, test_pvs, probability_bytes);
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        done_sequence_ub, done_sequence, next_no, xover_rows_ub,
        xover_slots_ub, current_xover, xover_list, test_pvs);
    if (invocation == 1) {
      file = CreateFileA(
          "make-test-pvs-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, MTP_DONE_SEQUENCE_OUT, done_sequence, done_bytes);
        write_section(file, MTP_CURRENT_XOVER_OUT, current_xover, current_bytes);
        write_section(file, MTP_XOVER_LIST_OUT, xover_list, list_bytes);
        write_section(file, MTP_TEST_PVS_OUT, test_pvs, probability_bytes);
        write_section(file, MTP_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL FindBestRecSignalP2Capture(
    char done_target, int next_no, int probability_rows_ub,
    int probability_columns_ub, double *low_p, char *done_sequence,
    int *trace, short *current_xover, double *test_pvs) {
  static FindBestRecSignalP2Fn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 done_bytes =
      (u32)((next_no + 1) * (probability_rows_ub + 1));
  const u32 current_bytes = (u32)((next_no + 1) * (int)sizeof(short));
  const u32 probability_bytes =
      (u32)((probability_rows_ub + 1) * (probability_columns_ub + 1) *
            (int)sizeof(double));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) {
      original = (FindBestRecSignalP2Fn)GetProcAddress(module, "FindBestRecSignalP2");
    }
  }
  if (!original) return 0;
  ++invocation;
  if (invocation == 1) {
    struct FindBestRecSignalP2Header header;
    const char magic[8] = {'F', 'B', 'R', 'S', 'I', 'G', '2', '\0'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.done_target = done_target;
    header.next_no = next_no;
    header.probability_rows_ub = probability_rows_ub;
    header.probability_columns_ub = probability_columns_ub;
    file = CreateFileA(
        "find-best-rec-signal-p2-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, FBR_LOW_P_IN, low_p, (u32)sizeof(double));
      write_section(file, FBR_DONE_SEQUENCE_IN, done_sequence, done_bytes);
      write_section(file, FBR_TRACE_IN, trace, (u32)(2 * (int)sizeof(int)));
      write_section(file, FBR_CURRENT_XOVER_IN, current_xover, current_bytes);
      write_section(file, FBR_TEST_PVS_IN, test_pvs, probability_bytes);
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        done_target, next_no, probability_rows_ub, probability_columns_ub,
        low_p, done_sequence, trace, current_xover, test_pvs);
    if (invocation == 1) {
      file = CreateFileA(
          "find-best-rec-signal-p2-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, FBR_LOW_P_OUT, low_p, (u32)sizeof(double));
        write_section(file, FBR_DONE_SEQUENCE_OUT, done_sequence, done_bytes);
        write_section(file, FBR_TRACE_OUT, trace, (u32)(2 * (int)sizeof(int)));
        write_section(file, FBR_CURRENT_XOVER_OUT, current_xover, current_bytes);
        write_section(file, FBR_TEST_PVS_OUT, test_pvs, probability_bytes);
        write_section(file, FBR_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL UFDistCapture(
    int sequence_length, int begin, int end, int pair_matrix_ub,
    float *valid, float *differences, float *breakpoint_distance,
    float *remainder_distance, int *sequences, int sequence_data_ub,
    short *sequence_data) {
  static UFDistFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 pair_bytes = (u32)((pair_matrix_ub + 1) * (pair_matrix_ub + 1) *
                               (int)sizeof(float));
  const u32 sequence_bytes =
      (u32)((sequence_data_ub + 1) * (pair_matrix_ub + 1) *
            (int)sizeof(short));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (UFDistFn)GetProcAddress(module, "UFDist");
  }
  if (!original) return 0;
  ++invocation;
  if (invocation == 1) {
    struct UFDistHeader header;
    const char magic[8] = {'U', 'F', 'D', 'I', 'S', 'T', '\0', '\0'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.sequence_length = sequence_length;
    header.begin = begin;
    header.end = end;
    header.pair_matrix_ub = pair_matrix_ub;
    header.sequence_data_ub = sequence_data_ub;
    file = CreateFileA(
        "ufdist-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, UFD_VALID_IN, valid, pair_bytes);
      write_section(file, UFD_DIFFERENCES_IN, differences, pair_bytes);
      write_section(
          file, UFD_BREAKPOINT_DISTANCE_IN, breakpoint_distance,
          (u32)(3 * (int)sizeof(float)));
      write_section(
          file, UFD_REMAINDER_DISTANCE_IN, remainder_distance,
          (u32)(3 * (int)sizeof(float)));
      write_section(file, UFD_SEQUENCES_IN, sequences, (u32)(3 * (int)sizeof(int)));
      write_section(file, UFD_SEQUENCE_DATA_IN, sequence_data, sequence_bytes);
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        sequence_length, begin, end, pair_matrix_ub, valid, differences,
        breakpoint_distance, remainder_distance, sequences, sequence_data_ub,
        sequence_data);
    if (invocation == 1) {
      file = CreateFileA(
          "ufdist-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(
            file, UFD_BREAKPOINT_DISTANCE_OUT, breakpoint_distance,
            (u32)(3 * (int)sizeof(float)));
        write_section(
            file, UFD_REMAINDER_DISTANCE_OUT, remainder_distance,
            (u32)(3 * (int)sizeof(float)));
        write_section(file, UFD_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int super_dist_p2_write_inputs(
    HANDLE file, int next_no, int ub14, int ub04, int ub13, int ub03,
    int ub12, int ub02, int ub11, double *average, float *pair_diff,
    float *pair_valid, float *distance, short *redo, int *category_count,
    short *iseq14, short *iseq04, short *iseq13, short *iseq03,
    short *iseq12, short *iseq02, short *iseq11, char *valid14,
    char *diff14, char *valid13, char *diff13, char *valid12, char *diff12,
    char *valid11, char *diff11, char *diff04, char *diff03, char *diff02) {
  const u32 matrix_bytes =
      (u32)((next_no + 1) * (next_no + 1) * (int)sizeof(float));
  write_section(file, SDP2_AVERAGE_IN, average, (u32)sizeof(double));
  write_section(file, SDP2_PAIR_DIFF_IN, pair_diff, matrix_bytes);
  write_section(file, SDP2_PAIR_VALID_IN, pair_valid, matrix_bytes);
  write_section(file, SDP2_DISTANCE_IN, distance, matrix_bytes);
  write_section(
      file, SDP2_REDO_IN, redo,
      (u32)((next_no + 1) * (int)sizeof(short)));
  write_section(
      file, SDP2_CATEGORY_COUNT_IN, category_count,
      (u32)(9 * (int)sizeof(int)));
#define WRITE_ISEQ(id, pointer, ub) \
  write_section(file, id, pointer, \
                (u32)((ub + 1) * (next_no + 1) * (int)sizeof(short)))
  WRITE_ISEQ(SDP2_ISEQ14_IN, iseq14, ub14);
  WRITE_ISEQ(SDP2_ISEQ04_IN, iseq04, ub04);
  WRITE_ISEQ(SDP2_ISEQ13_IN, iseq13, ub13);
  WRITE_ISEQ(SDP2_ISEQ03_IN, iseq03, ub03);
  WRITE_ISEQ(SDP2_ISEQ12_IN, iseq12, ub12);
  WRITE_ISEQ(SDP2_ISEQ02_IN, iseq02, ub02);
  WRITE_ISEQ(SDP2_ISEQ11_IN, iseq11, ub11);
#undef WRITE_ISEQ
  write_section(file, SDP2_VALID14_IN, valid14, (u32)(626 * 626));
  write_section(file, SDP2_DIFF14_IN, diff14, (u32)(626 * 626));
  write_section(file, SDP2_VALID13_IN, valid13, (u32)(1025 * 1025));
  write_section(file, SDP2_DIFF13_IN, diff13, (u32)(1025 * 1025));
  write_section(file, SDP2_VALID12_IN, valid12, (u32)(730 * 730));
  write_section(file, SDP2_DIFF12_IN, diff12, (u32)(730 * 730));
  write_section(file, SDP2_VALID11_IN, valid11, (u32)(1025 * 1025));
  write_section(file, SDP2_DIFF11_IN, diff11, (u32)(1025 * 1025));
  write_section(file, SDP2_DIFF04_IN, diff04, (u32)(1025 * 1025));
  write_section(file, SDP2_DIFF03_IN, diff03, (u32)(730 * 730));
  write_section(file, SDP2_DIFF02_IN, diff02, (u32)(1025 * 1025));
  return 1;
}

double STDCALL SuperDistPCapture(
    int x, int next_no, int ub14, int ub04, int ub13, int ub03, int ub12,
    int ub02, int ub11, double *average, float *pair_diff, float *pair_valid,
    float *distance, short *redo, int *category_count, short *iseq14,
    short *iseq04, short *iseq13, short *iseq03, short *iseq12,
    short *iseq02, short *iseq11, char *valid14, char *diff14,
    char *valid13, char *diff13, char *valid12, char *diff12, char *valid11,
    char *diff11, char *diff04, char *diff03, char *diff02) {
  static SuperDistPFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 matrix_bytes =
      (u32)((next_no + 1) * (next_no + 1) * (int)sizeof(float));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (SuperDistPFn)GetProcAddress(module, "SuperDistP");
  }
  if (!original) return 0.0;
  ++invocation;
  if (invocation == 1) {
    struct SuperDistP2Header header;
    const char magic[8] = {'S', 'U', 'P', 'D', 'I', 'S', 'T', '1'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.x = x;
    header.next_no = next_no;
    header.ub14 = ub14;
    header.ub04 = ub04;
    header.ub13 = ub13;
    header.ub03 = ub03;
    header.ub12 = ub12;
    header.ub02 = ub02;
    header.ub11 = ub11;
    file = CreateFileA(
        "super-dist-p-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      super_dist_p2_write_inputs(
          file, next_no, ub14, ub04, ub13, ub03, ub12, ub02, ub11, average,
          pair_diff, pair_valid, distance, redo, category_count, iseq14,
          iseq04, iseq13, iseq03, iseq12, iseq02, iseq11, valid14, diff14,
          valid13, diff13, valid12, diff12, valid11, diff11, diff04, diff03,
          diff02);
      CloseHandle(file);
    }
  }
  {
    const double result = original(
        x, next_no, ub14, ub04, ub13, ub03, ub12, ub02, ub11, average,
        pair_diff, pair_valid, distance, redo, category_count, iseq14, iseq04,
        iseq13, iseq03, iseq12, iseq02, iseq11, valid14, diff14, valid13,
        diff13, valid12, diff12, valid11, diff11, diff04, diff03, diff02);
    if (invocation == 1) {
      file = CreateFileA(
          "super-dist-p-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, SDP2_AVERAGE_OUT, average, (u32)sizeof(double));
        write_section(file, SDP2_PAIR_DIFF_OUT, pair_diff, matrix_bytes);
        write_section(file, SDP2_PAIR_VALID_OUT, pair_valid, matrix_bytes);
        write_section(file, SDP2_DISTANCE_OUT, distance, matrix_bytes);
        write_section(file, SDP2_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

double STDCALL SuperDistP2Capture(
    int x, int next_no, int ub14, int ub04, int ub13, int ub03, int ub12,
    int ub02, int ub11, double *average, float *pair_diff, float *pair_valid,
    float *distance, short *redo, int *category_count, short *iseq14,
    short *iseq04, short *iseq13, short *iseq03, short *iseq12,
    short *iseq02, short *iseq11, char *valid14, char *diff14,
    char *valid13, char *diff13, char *valid12, char *diff12, char *valid11,
    char *diff11, char *diff04, char *diff03, char *diff02) {
  static SuperDistP2Fn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 matrix_bytes =
      (u32)((next_no + 1) * (next_no + 1) * (int)sizeof(float));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (SuperDistP2Fn)GetProcAddress(module, "SuperDistP2");
  }
  if (!original) return 0.0;
  ++invocation;
  if (invocation == 1) {
    struct SuperDistP2Header header;
    const char magic[8] = {'S', 'U', 'P', 'D', 'I', 'S', 'T', '2'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.x = x;
    header.next_no = next_no;
    header.ub14 = ub14;
    header.ub04 = ub04;
    header.ub13 = ub13;
    header.ub03 = ub03;
    header.ub12 = ub12;
    header.ub02 = ub02;
    header.ub11 = ub11;
    file = CreateFileA(
        "super-dist-p2-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      super_dist_p2_write_inputs(
          file, next_no, ub14, ub04, ub13, ub03, ub12, ub02, ub11, average,
          pair_diff, pair_valid, distance, redo, category_count, iseq14,
          iseq04, iseq13, iseq03, iseq12, iseq02, iseq11, valid14, diff14,
          valid13, diff13, valid12, diff12, valid11, diff11, diff04, diff03,
          diff02);
      CloseHandle(file);
    }
  }
  {
    const double result = original(
        x, next_no, ub14, ub04, ub13, ub03, ub12, ub02, ub11, average,
        pair_diff, pair_valid, distance, redo, category_count, iseq14, iseq04,
        iseq13, iseq03, iseq12, iseq02, iseq11, valid14, diff14, valid13,
        diff13, valid12, diff12, valid11, diff11, diff04, diff03, diff02);
    if (invocation == 1) {
      file = CreateFileA(
          "super-dist-p2-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, SDP2_AVERAGE_OUT, average, (u32)sizeof(double));
        write_section(file, SDP2_PAIR_DIFF_OUT, pair_diff, matrix_bytes);
        write_section(file, SDP2_PAIR_VALID_OUT, pair_valid, matrix_bytes);
        write_section(file, SDP2_DISTANCE_OUT, distance, matrix_bytes);
        write_section(file, SDP2_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL CheckMatrixPCapture(
    int *minimums, int *sequences, int next_no, int sco,
    int minimum_sequence_size, int missing_pair_ub,
    unsigned char *missing_pair, int valid_ub, float *valid,
    int sub_valid_ub, float *sub_valid, int matrix_ub, float *first_matrix,
    float *second_matrix, int *first_total, int *second_total) {
  static CheckMatrixPFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 minimum_bytes = (u32)((next_no + 1) * (int)sizeof(int));
  const u32 missing_bytes =
      (u32)((missing_pair_ub + 1) * (missing_pair_ub + 1));
  const u32 valid_bytes =
      (u32)((valid_ub + 1) * (valid_ub + 1) * (int)sizeof(float));
  const u32 sub_valid_bytes =
      (u32)((sub_valid_ub + 1) * (sub_valid_ub + 1) * (int)sizeof(float));
  const u32 matrix_bytes =
      (u32)((matrix_ub + 1) * (matrix_ub + 1) * (int)sizeof(float));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (CheckMatrixPFn)GetProcAddress(module, "CheckMatrixP");
  }
  if (!original) return 0;
  ++invocation;
  if (invocation == 1) {
    struct CheckMatrixPHeader header;
    const char magic[8] = {'C', 'H', 'K', 'M', 'A', 'T', 'P', '\0'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.next_no = next_no;
    header.sco = sco;
    header.minimum_sequence_size = minimum_sequence_size;
    header.missing_pair_ub = missing_pair_ub;
    header.valid_ub = valid_ub;
    header.sub_valid_ub = sub_valid_ub;
    header.matrix_ub = matrix_ub;
    file = CreateFileA(
        "check-matrix-p-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, CMP_MINIMUMS_IN, minimums, minimum_bytes);
      write_section(file, CMP_SEQUENCES_IN, sequences, (u32)(3 * (int)sizeof(int)));
      write_section(file, CMP_MISSING_PAIR_IN, missing_pair, missing_bytes);
      write_section(file, CMP_VALID_IN, valid, valid_bytes);
      write_section(file, CMP_SUB_VALID_IN, sub_valid, sub_valid_bytes);
      write_section(file, CMP_FIRST_MATRIX_IN, first_matrix, matrix_bytes);
      write_section(file, CMP_SECOND_MATRIX_IN, second_matrix, matrix_bytes);
      write_section(file, CMP_FIRST_TOTAL_IN, first_total, minimum_bytes);
      write_section(file, CMP_SECOND_TOTAL_IN, second_total, minimum_bytes);
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        minimums, sequences, next_no, sco, minimum_sequence_size,
        missing_pair_ub, missing_pair, valid_ub, valid, sub_valid_ub,
        sub_valid, matrix_ub, first_matrix, second_matrix, first_total,
        second_total);
    if (invocation == 1) {
      file = CreateFileA(
          "check-matrix-p-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, CMP_MINIMUMS_OUT, minimums, minimum_bytes);
        write_section(file, CMP_MISSING_PAIR_OUT, missing_pair, missing_bytes);
        write_section(file, CMP_FIRST_MATRIX_OUT, first_matrix, matrix_bytes);
        write_section(file, CMP_SECOND_MATRIX_OUT, second_matrix, matrix_bytes);
        write_section(file, CMP_FIRST_TOTAL_OUT, first_total, minimum_bytes);
        write_section(file, CMP_SECOND_TOTAL_OUT, second_total, minimum_bytes);
        write_section(file, CMP_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL MakeNJTreesP2Capture(
    int resolve_root, int nseqs, int next_no, int *sequences,
    unsigned char *min_pair, unsigned char *sequence_pair, int seed,
    int name_length, int sequence_length, int trace_sequences_ub,
    int *outlier, int *trace_sequences, int first_matrix_ub,
    float *first_matrix, int second_matrix_ub, float *second_matrix,
    int first_adjusted_matrix_ub, float *first_adjusted_matrix,
    int second_adjusted_matrix_ub, float *second_adjusted_matrix,
    int *redo_list, char *first_holder, char *second_holder,
    float *temp_first_matrix, float *temp_second_matrix) {
  static MakeNJTreesP2Fn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 trace_bytes =
      (u32)((trace_sequences_ub + 1) * (next_no + 1) * (int)sizeof(int));
  const u32 first_bytes =
      (u32)((first_matrix_ub + 1) * (first_matrix_ub + 1) * (int)sizeof(float));
  const u32 second_bytes =
      (u32)((second_matrix_ub + 1) * (second_matrix_ub + 1) * (int)sizeof(float));
  const u32 first_adjusted_bytes =
      (u32)((first_adjusted_matrix_ub + 1) * (first_adjusted_matrix_ub + 1) *
            (int)sizeof(float));
  const u32 second_adjusted_bytes =
      (u32)((second_adjusted_matrix_ub + 1) * (second_adjusted_matrix_ub + 1) *
            (int)sizeof(float));
  const u32 temp_bytes =
      (u32)((nseqs + 1) * (nseqs + 1) * (int)sizeof(float));
  /* VB6 ReDim's argument is an inclusive upper bound. */
  const u32 holder_bytes = (u32)((nseqs + 1) * 80 + 1);
  const u32 redo_bytes = (u32)((next_no + 1) * (int)sizeof(int));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (MakeNJTreesP2Fn)GetProcAddress(module, "MakeNJTreesP2");
  }
  if (!original) return 0;
  ++invocation;
#define WRITE_MNJ_STATE(suffix) \
  write_section(file, MNJ_SEQUENCES_##suffix, sequences, (u32)(3 * (int)sizeof(int))); \
  write_section(file, MNJ_MIN_PAIR_##suffix, min_pair, 3U); \
  write_section(file, MNJ_SEQUENCE_PAIR_##suffix, sequence_pair, 3U); \
  write_section(file, MNJ_OUTLIER_##suffix, outlier, (u32)(3 * (int)sizeof(int))); \
  write_section(file, MNJ_TRACE_SEQUENCES_##suffix, trace_sequences, trace_bytes); \
  write_section(file, MNJ_FIRST_MATRIX_##suffix, first_matrix, first_bytes); \
  write_section(file, MNJ_SECOND_MATRIX_##suffix, second_matrix, second_bytes); \
  write_section(file, MNJ_FIRST_ADJUSTED_MATRIX_##suffix, first_adjusted_matrix, first_adjusted_bytes); \
  write_section(file, MNJ_SECOND_ADJUSTED_MATRIX_##suffix, second_adjusted_matrix, second_adjusted_bytes); \
  write_section(file, MNJ_REDO_LIST_##suffix, redo_list, redo_bytes); \
  write_section(file, MNJ_FIRST_HOLDER_##suffix, first_holder, holder_bytes); \
  write_section(file, MNJ_SECOND_HOLDER_##suffix, second_holder, holder_bytes); \
  write_section(file, MNJ_TEMP_FIRST_MATRIX_##suffix, temp_first_matrix, temp_bytes); \
  write_section(file, MNJ_TEMP_SECOND_MATRIX_##suffix, temp_second_matrix, temp_bytes)
  if (invocation == 1) {
    struct MakeNJTreesP2Header header;
    const char magic[8] = {'M', 'A', 'K', 'E', 'N', 'J', 'P', '2'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.resolve_root = resolve_root;
    header.nseqs = nseqs;
    header.next_no = next_no;
    header.seed = seed;
    header.name_length = name_length;
    header.sequence_length = sequence_length;
    header.trace_sequences_ub = trace_sequences_ub;
    header.first_matrix_ub = first_matrix_ub;
    header.second_matrix_ub = second_matrix_ub;
    header.first_adjusted_matrix_ub = first_adjusted_matrix_ub;
    header.second_adjusted_matrix_ub = second_adjusted_matrix_ub;
    file = CreateFileA(
        "make-nj-trees-p2-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      WRITE_MNJ_STATE(IN);
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        resolve_root, nseqs, next_no, sequences, min_pair, sequence_pair, seed,
        name_length, sequence_length, trace_sequences_ub, outlier,
        trace_sequences, first_matrix_ub, first_matrix, second_matrix_ub,
        second_matrix, first_adjusted_matrix_ub, first_adjusted_matrix,
        second_adjusted_matrix_ub, second_adjusted_matrix, redo_list,
        first_holder, second_holder, temp_first_matrix, temp_second_matrix);
    if (invocation == 1) {
      file = CreateFileA(
          "make-nj-trees-p2-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        WRITE_MNJ_STATE(OUT);
        write_section(file, MNJ_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
#undef WRITE_MNJ_STATE
    return result;
  }
}

int STDCALL MarkOutsidesCapture(
    int done_sequence_ub, unsigned char *done_sequence, int next_no,
    int xover_rows_ub, short *current_xover, XOVERDEFINE *xover_list) {
  static MarkOutsidesFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  int maximum_current_xover = 0;
  int x;
  u32 done_bytes;
  u32 current_bytes;
  u32 xover_bytes;
  for (x = 0; x <= next_no; ++x) {
    if (current_xover[x] > maximum_current_xover)
      maximum_current_xover = current_xover[x];
  }
  done_bytes = (u32)((done_sequence_ub + 1) *
                     (maximum_current_xover + 1));
  current_bytes = (u32)((next_no + 1) * (int)sizeof(short));
  xover_bytes = (u32)((xover_rows_ub + 1) *
                      (maximum_current_xover + 1) *
                      (int)sizeof(XOVERDEFINE));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module)
      original = (MarkOutsidesFn)GetProcAddress(module, "MarkOutsides");
  }
  if (!original) return 0;
  ++invocation;
  if (invocation == 1) {
    struct MarkOutsidesHeader header;
    const char magic[8] = {'M', 'A', 'R', 'K', 'O', 'U', 'T', 'S'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.done_sequence_ub = done_sequence_ub;
    header.next_no = next_no;
    header.xover_rows_ub = xover_rows_ub;
    header.maximum_current_xover = maximum_current_xover;
    header.xover_struct_bytes = (int)sizeof(XOVERDEFINE);
    file = CreateFileA(
        "mark-outsides-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, MO_DONE_SEQUENCE_IN, done_sequence, done_bytes);
      write_section(file, MO_CURRENT_XOVER_IN, current_xover, current_bytes);
      write_section(file, MO_XOVER_LIST_IN, xover_list, xover_bytes);
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        done_sequence_ub, done_sequence, next_no, xover_rows_ub,
        current_xover, xover_list);
    if (invocation == 1) {
      file = CreateFileA(
          "mark-outsides-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, MO_DONE_SEQUENCE_OUT, done_sequence, done_bytes);
        write_section(file, MO_CURRENT_XOVER_OUT, current_xover, current_bytes);
        write_section(file, MO_XOVER_LIST_OUT, xover_list, xover_bytes);
        write_section(file, MO_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL MakeSDMP2Capture(
    int next_no, int sequence_length, int *start_positions,
    int *end_positions, int *sequences, int *comparison_matrix,
    unsigned char *missing_data, short *sequence_data, double *summary_matrix,
    double *distance_matrix) {
  static MakeSDMP2Fn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const u32 alignment_cells =
      (u32)((next_no + 1) * (sequence_length + 1));
  const u32 summary_bytes =
      (u32)(9 * (next_no + 1) * (int)sizeof(double));
  const u32 distance_bytes =
      (u32)(45 * (next_no + 1) * (int)sizeof(double));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (MakeSDMP2Fn)GetProcAddress(module, "MakeSDMP2");
  }
  if (!original) return 0;
  ++invocation;
  if (invocation == 1) {
    struct MakeSDMP2Header header;
    const char magic[8] = {'M', 'A', 'K', 'E', 'S', 'D', 'M', 'P'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.next_no = next_no;
    header.sequence_length = sequence_length;
    file = CreateFileA(
        "make-sdmp2-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, MSD_START_POSITIONS_IN, start_positions,
                    (u32)(5 * (int)sizeof(int)));
      write_section(file, MSD_END_POSITIONS_IN, end_positions,
                    (u32)(5 * (int)sizeof(int)));
      write_section(file, MSD_SEQUENCES_IN, sequences,
                    (u32)(3 * (int)sizeof(int)));
      write_section(file, MSD_COMPARISON_MATRIX_IN, comparison_matrix,
                    (u32)(6 * (int)sizeof(int)));
      write_section(file, MSD_MISSING_DATA_IN, missing_data, alignment_cells);
      write_section(file, MSD_SEQUENCE_DATA_IN, sequence_data,
                    (u32)(alignment_cells * (int)sizeof(short)));
      write_section(file, MSD_SUMMARY_MATRIX_IN, summary_matrix, summary_bytes);
      write_section(file, MSD_DISTANCE_MATRIX_IN, distance_matrix,
                    distance_bytes);
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        next_no, sequence_length, start_positions, end_positions, sequences,
        comparison_matrix, missing_data, sequence_data, summary_matrix,
        distance_matrix);
    if (invocation == 1) {
      file = CreateFileA(
          "make-sdmp2-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, MSD_SUMMARY_MATRIX_OUT, summary_matrix,
                      summary_bytes);
        write_section(file, MSD_DISTANCE_MATRIX_OUT, distance_matrix,
                      distance_bytes);
        write_section(file, MSD_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL FillRmatCapture(
    int y, int next_no, int result_matrix_ub1, int result_matrix_ub2,
    int distance_matrix_ub1, int distance_matrix_ub2,
    int distance_matrix_ub3, double *result_matrix, double *distance_matrix,
    unsigned char *positions) {
  static FillRmatFn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  const char *file_name;
  const u32 result_bytes =
      (u32)((result_matrix_ub1 + 1) * (result_matrix_ub2 + 1) *
            (next_no + 1) * (int)sizeof(double));
  const u32 distance_bytes =
      (u32)(3 * (distance_matrix_ub1 + 1) * (distance_matrix_ub2 + 1) *
            (distance_matrix_ub3 + 1) * (int)sizeof(double));
  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (FillRmatFn)GetProcAddress(module, "FillRmat");
  }
  if (!original) return 0;
  ++invocation;
  if (invocation == 1) file_name = "fill-rmat-y0-v1.bin";
  else if (invocation == 2) file_name = "fill-rmat-y1-v1.bin";
  else if (invocation == 3) file_name = "fill-rmat-y2-v1.bin";
  else file_name = (const char *)0;
  if (file_name) {
    struct FillRmatHeader header;
    const char magic[8] = {'F', 'I', 'L', 'L', 'R', 'M', 'A', 'T'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.y = y;
    header.next_no = next_no;
    header.result_matrix_ub1 = result_matrix_ub1;
    header.result_matrix_ub2 = result_matrix_ub2;
    header.distance_matrix_ub1 = distance_matrix_ub1;
    header.distance_matrix_ub2 = distance_matrix_ub2;
    header.distance_matrix_ub3 = distance_matrix_ub3;
    file = CreateFileA(
        file_name, GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
        (void *)0, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, FRM_RESULT_MATRIX_IN, result_matrix, result_bytes);
      write_section(file, FRM_DISTANCE_MATRIX_IN, distance_matrix,
                    distance_bytes);
      write_section(file, FRM_POSITIONS_IN, positions, 2U);
      CloseHandle(file);
    }
  }
  {
    const int result = original(
        y, next_no, result_matrix_ub1, result_matrix_ub2,
        distance_matrix_ub1, distance_matrix_ub2, distance_matrix_ub3,
        result_matrix, distance_matrix, positions);
    if (file_name) {
      file = CreateFileA(
          file_name, GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
          (void *)0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end_marker = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, FRM_RESULT_MATRIX_OUT, result_matrix, result_bytes);
        write_section(file, FRM_DISTANCE_MATRIX_OUT, distance_matrix,
                      distance_bytes);
        write_section(file, FRM_POSITIONS_OUT, positions, 2U);
        write_section(file, FRM_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end_marker, (DWORD)sizeof(end_marker));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL FindSubSeqPB3Capture(
    int *ah, int fss_ub, int xover_window, int compressed_sequence_ub,
    int sequence_length, int next_no, int seq1, int seq2, int seq3,
    unsigned char *compressed_sequence, int xover_sequence_ub,
    char *xover_sequence, unsigned char *fss_rdp) {
  static FindSubSeqPB3Fn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;

  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) {
      original = (FindSubSeqPB3Fn)GetProcAddress(module, "FindSubSeqPB3");
    }
  }
  if (!original) return 0;
  ++invocation;

  if (invocation == 1) {
    struct FindSubSeqPB3Header header;
    const char magic[8] = {'F', 'S', 'P', 'B', '3', '\0', '\0', '\0'};
    const int fss_width = fss_ub + 1;
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.fss_ub = fss_ub;
    header.xover_window = xover_window;
    header.compressed_sequence_ub = compressed_sequence_ub;
    header.sequence_length = sequence_length;
    header.next_no = next_no;
    header.seq1 = seq1;
    header.seq2 = seq2;
    header.seq3 = seq3;
    header.xover_sequence_ub = xover_sequence_ub;
    file = CreateFileA(
        "find-subseq-pb3-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, FSP_AH_IN, ah, (u32)(4 * (int)sizeof(int)));
      write_section(
          file, FSP_COMPRESSED_SEQUENCE_IN, compressed_sequence,
          (u32)((compressed_sequence_ub + 1) * (next_no + 1)));
      write_section(
          file, FSP_XOVER_SEQUENCE_IN, xover_sequence,
          (u32)((xover_sequence_ub + 1) * 3));
      write_section(
          file, FSP_FSS_RDP_IN, fss_rdp,
          (u32)(4 * fss_width * fss_width * fss_width));
      CloseHandle(file);
    }
  }

  {
    const int result = original(
        ah, fss_ub, xover_window, compressed_sequence_ub, sequence_length,
        next_no, seq1, seq2, seq3, compressed_sequence, xover_sequence_ub,
        xover_sequence, fss_rdp);
    if (invocation == 1) {
      file = CreateFileA(
          "find-subseq-pb3-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, FSP_AH_OUT, ah, (u32)(4 * (int)sizeof(int)));
        write_section(
            file, FSP_XOVER_SEQUENCE_OUT, xover_sequence,
            (u32)((xover_sequence_ub + 1) * 3));
        write_section(file, FSP_RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end, (DWORD)sizeof(end));
        CloseHandle(file);
      }
    }
    return result;
  }
}

int STDCALL AlistRDP4Capture(
    int store_lpv_ub, double *store_lpv, short *analysis_list,
    int list_length, int start, int end, int next_no, double sub_threshold,
    unsigned char *redo_list, int circular, int mc_correction, int mc_flag,
    double lowest_probability, int target_x, int sequence_length,
    int short_output, int distance_ub, float *distance,
    int tree_distance_ub, float *tree_distance, int fss_rdp_ub,
    int compressed_sequence_ub, unsigned char *compressed_sequence,
    short *sequence_data, int xover_window, short xover_window_x,
    unsigned char *fss_rdp, int probability_file_flag,
    int probability_one_ub, int probability_two_ub,
    double *probability_estimate, int fact_three_ub, double *fact_three,
    double *fact) {
  static AlistRDP4Fn original;
  static int invocation;
  HANDLE file = INVALID_HANDLE_VALUE;
  struct CaptureHeader header;
  u32 store_bytes;
  u32 redo_bytes;

  if (!original) {
    HMODULE module = LoadLibraryA("DNA5_original.dll");
    if (module) original = (AlistRDP4Fn)GetProcAddress(module, "AlistRDP4");
  }
  if (!original) return 0;

  ++invocation;
  store_bytes = (u32)((store_lpv_ub + 1) * (next_no + 1) * (int)sizeof(double));
  redo_bytes = (u32)(list_length + 1);

  if (invocation == 1) {
    const int fss_width = fss_rdp_ub + 1;
    const int fact_width = fact_three_ub + 1;
    const int probability_width_one = probability_one_ub + 1;
    const int probability_width_two = probability_two_ub + 1;
    const char magic[8] = {'A', 'L', 'R', 'D', 'P', '4', '\0', '\0'};
    int i;
    for (i = 0; i < 8; ++i) header.magic[i] = magic[i];
    header.version = 1;
    header.store_lpv_ub = store_lpv_ub;
    header.list_length = list_length;
    header.start = start;
    header.end = end;
    header.next_no = next_no;
    header.sub_threshold = sub_threshold;
    header.circular = circular;
    header.mc_correction = mc_correction;
    header.mc_flag = mc_flag;
    header.lowest_probability = lowest_probability;
    header.target_x = target_x;
    header.sequence_length = sequence_length;
    header.short_output = short_output;
    header.distance_ub = distance_ub;
    header.tree_distance_ub = tree_distance_ub;
    header.fss_rdp_ub = fss_rdp_ub;
    header.compressed_sequence_ub = compressed_sequence_ub;
    header.xover_window = xover_window;
    header.xover_window_x = xover_window_x;
    header.probability_file_flag = probability_file_flag;
    header.probability_one_ub = probability_one_ub;
    header.probability_two_ub = probability_two_ub;
    header.fact_three_ub = fact_three_ub;

    file = CreateFileA(
        "alist-rdp4-v1.bin", GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
    if (file != INVALID_HANDLE_VALUE) {
      write_bytes(file, &header, (DWORD)sizeof(header));
      write_section(file, STORE_LPV_IN, store_lpv, store_bytes);
      write_section(
          file, ANALYSIS_LIST_IN, analysis_list,
          (u32)(3 * (list_length + 1) * (int)sizeof(short)));
      write_section(file, REDO_LIST_IN, redo_list, redo_bytes);
      write_section(
          file, DISTANCE_IN, distance,
          (u32)((distance_ub + 1) * (distance_ub + 1) * (int)sizeof(float)));
      write_section(
          file, TREE_DISTANCE_IN, tree_distance,
          (u32)((tree_distance_ub + 1) * (tree_distance_ub + 1) *
                (int)sizeof(float)));
      write_section(
          file, FSS_RDP_IN, fss_rdp,
          (u32)(4 * fss_width * fss_width * fss_width));
      write_section(
          file, COMPRESSED_SEQUENCE_IN, compressed_sequence,
          (u32)((compressed_sequence_ub + 1) * (next_no + 1)));
      /* Seqnum is normally ReDim(length, NextNo + 1) at this call site. */
      write_section(
          file, SEQUENCE_DATA_IN, sequence_data,
          (u32)((sequence_length + 1) * (next_no + 2) * (int)sizeof(short)));
      /* VB declares ProbEstimate(UBPE1, UBPE2, 50); the third UB is fixed. */
      write_section(
          file, PROBABILITY_ESTIMATE_IN, probability_estimate,
          (u32)(probability_width_one * probability_width_two * 51 *
                (int)sizeof(double)));
      write_section(
          file, FACT_THREE_IN, fact_three,
          (u32)(fact_width * fact_width * fact_width * (int)sizeof(double)));
      /* RDP initializes this globally as ReDim Fact(171). */
      write_section(file, FACT_IN, fact, (u32)(172 * (int)sizeof(double)));
      CloseHandle(file);
    }
  }

  {
    const int result = original(
        store_lpv_ub, store_lpv, analysis_list, list_length, start, end,
        next_no, sub_threshold, redo_list, circular, mc_correction, mc_flag,
        lowest_probability, target_x, sequence_length, short_output,
        distance_ub, distance, tree_distance_ub, tree_distance, fss_rdp_ub,
        compressed_sequence_ub, compressed_sequence, sequence_data,
        xover_window, xover_window_x, fss_rdp, probability_file_flag,
        probability_one_ub, probability_two_ub, probability_estimate,
        fact_three_ub, fact_three, fact);

    if (invocation == 1) {
      file = CreateFileA(
          "alist-rdp4-v1.bin", GENERIC_WRITE,
          FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0, OPEN_EXISTING,
          FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
      if (file != INVALID_HANDLE_VALUE) {
        const struct SectionHeader end = {END_MARKER, 0};
        SetFilePointer(file, 0, (long *)0, FILE_END);
        write_section(file, REDO_LIST_OUT, redo_list, redo_bytes);
        write_section(file, STORE_LPV_OUT, store_lpv, store_bytes);
        write_section(file, RESULT_OUT, &result, (u32)sizeof(result));
        write_bytes(file, &end, (DWORD)sizeof(end));
        CloseHandle(file);
      }
    }
    return result;
  }
}
