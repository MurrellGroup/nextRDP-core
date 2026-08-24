typedef void *HANDLE;
typedef void *HMODULE;
typedef void *FARPROC;
typedef unsigned long DWORD;
typedef int BOOL;

int _fltused;

#define STDCALL __attribute__((stdcall))
#define DLLIMPORT __declspec(dllimport)
#define GENERIC_WRITE 0x40000000UL
#define FILE_SHARE_READ 0x00000001UL
#define FILE_SHARE_WRITE 0x00000002UL
#define OPEN_ALWAYS 4UL
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

typedef int(STDCALL *MakeVarMap2Fn)(
    int, int, int, short *, short *, int *, int *, int *, int *);
typedef int(STDCALL *MakeCntHit2Fn)(
    int, int, int, int, int, int, int, double *, short *, float *, int *);

static void write_bytes(HANDLE file, const void *data, DWORD bytes) {
  DWORD written = 0;
  WriteFile(file, data, bytes, &written, (void *)0);
}

static HANDLE open_trace(const char *name) {
  HANDLE file = CreateFileA(
      name, GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0,
      OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
  if (file != INVALID_HANDLE_VALUE) {
    SetFilePointer(file, 0, (long *)0, FILE_END);
  }
  return file;
}

int STDCALL MakeVarMap2Capture(
    int next_no, int sequence_length, int map_upper_bound,
    short *sequence_data, short *variable_site_map,
    int *variable_region_positions, int *variable_positions,
    int *sequences, int *comparison_matrix) {
  static MakeVarMap2Fn original;
  static unsigned int invocation;
  if (!original) {
    HMODULE module = LoadLibraryA("dna_original.dll");
    if (module) {
      original = (MakeVarMap2Fn)GetProcAddress(module, "MakeVarMap2");
    }
  }
  if (!original) return 0;

  ++invocation;
  const int count = next_no + 1;
  HANDLE file = open_trace("make-var-map2-trace.bin");
  if (file != INVALID_HANDLE_VALUE) {
    const unsigned int header[] = {
        0x32564d43U, invocation, (unsigned int)next_no,
        (unsigned int)sequence_length, (unsigned int)map_upper_bound};
    write_bytes(file, header, sizeof(header));
    write_bytes(file, sequences, 3U * sizeof(int));
    write_bytes(file, comparison_matrix, 6U * sizeof(int));
    write_bytes(
        file, sequence_data,
        (DWORD)((sequence_length + 1) * count * (int)sizeof(short)));
  }

  const int result = original(
      next_no, sequence_length, map_upper_bound, sequence_data,
      variable_site_map, variable_region_positions, variable_positions,
      sequences, comparison_matrix);
  if (file != INVALID_HANDLE_VALUE) {
    write_bytes(file, &result, sizeof(result));
    write_bytes(
        file, variable_site_map,
        (DWORD)(3 * (map_upper_bound + 1) * count * (int)sizeof(short)));
    write_bytes(
        file, variable_region_positions,
        (DWORD)((sequence_length + 1) * (int)sizeof(int)));
    write_bytes(
        file, variable_positions,
        (DWORD)((sequence_length + 1) * (int)sizeof(int)));
    CloseHandle(file);
  }
  return result;
}

int STDCALL MakeCntHit2Capture(
    int beginning, int ending, int smoothing_window, int next_no,
    int variable_sites, int sequence_length, int map_upper_bound,
    double *count_hits, short *variable_site_map,
    float *variable_site_smooth, int *variable_region_positions) {
  static MakeCntHit2Fn original;
  static unsigned int invocation;
  if (!original) {
    HMODULE module = LoadLibraryA("dna_original.dll");
    if (module) {
      original = (MakeCntHit2Fn)GetProcAddress(module, "MakeCntHit2");
    }
  }
  if (!original) return 0;

  ++invocation;
  const int count = next_no + 1;
  HANDLE file = open_trace("make-cnt-hit2-trace.bin");
  if (file != INVALID_HANDLE_VALUE) {
    const unsigned int header[] = {
        0x32484343U, invocation, (unsigned int)beginning,
        (unsigned int)ending, (unsigned int)smoothing_window,
        (unsigned int)next_no, (unsigned int)variable_sites,
        (unsigned int)sequence_length, (unsigned int)map_upper_bound};
    write_bytes(file, header, sizeof(header));
    write_bytes(file, count_hits, (DWORD)(6 * count * (int)sizeof(double)));
    write_bytes(
        file, variable_site_map,
        (DWORD)(3 * (map_upper_bound + 1) * count * (int)sizeof(short)));
    write_bytes(
        file, variable_region_positions,
        (DWORD)((sequence_length + 1) * (int)sizeof(int)));
  }

  const int result = original(
      beginning, ending, smoothing_window, next_no, variable_sites,
      sequence_length, map_upper_bound, count_hits, variable_site_map,
      variable_site_smooth, variable_region_positions);
  if (file != INVALID_HANDLE_VALUE) {
    write_bytes(file, &result, sizeof(result));
    write_bytes(file, count_hits, (DWORD)(6 * count * (int)sizeof(double)));
    write_bytes(
        file, variable_site_smooth,
        (DWORD)(3 * (variable_sites + 1) * count * (int)sizeof(float)));
    CloseHandle(file);
  }
  return result;
}
