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
