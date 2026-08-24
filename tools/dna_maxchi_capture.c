/* Scalar boundary tracer for the old DNA.dll helpers called by MCXoverF. */

typedef void *HANDLE;
typedef void *HMODULE;
typedef void *FARPROC;
typedef unsigned long DWORD;

#define STDCALL __attribute__((stdcall))
#define DLLIMPORT __attribute__((dllimport))
#define INVALID_HANDLE_VALUE ((HANDLE)(long)-1)
#define GENERIC_WRITE 0x40000000UL
#define FILE_SHARE_READ 1UL
#define FILE_SHARE_WRITE 2UL
#define CREATE_ALWAYS 2UL
#define OPEN_EXISTING 3UL
#define FILE_ATTRIBUTE_NORMAL 0x80UL
#define FILE_END 2UL

DLLIMPORT HMODULE STDCALL LoadLibraryA(const char *name);
DLLIMPORT FARPROC STDCALL GetProcAddress(HMODULE module, const char *name);
DLLIMPORT HANDLE STDCALL CreateFileA(
    const char *, DWORD, DWORD, void *, DWORD, DWORD, HANDLE);
DLLIMPORT int STDCALL WriteFile(HANDLE, const void *, DWORD, DWORD *, void *);
DLLIMPORT DWORD STDCALL SetFilePointer(HANDLE, long, long *, DWORD);
DLLIMPORT int STDCALL CloseHandle(HANDLE);

typedef int(STDCALL *FindSideFn)(
    int, int, int, int, int, int, int, int, char *, double *, double *);
typedef int(STDCALL *OptFn)(
    int, double, int, int, int, int, int, int, char *, char *);
typedef int(STDCALL *DestroyPeaksFn)(
    int, int, int, int, int, double *, double *);

#pragma pack(push, 1)
struct TraceRecord {
  int function;
  int values[12];
  double doubles[3];
};
#pragma pack(pop)

static int trace_invocation;

static void append_record(struct TraceRecord *record) {
  HANDLE file = CreateFileA(
      "dna-maxchi-helpers.bin", GENERIC_WRITE,
      FILE_SHARE_READ | FILE_SHARE_WRITE, (void *)0,
      trace_invocation == 0 ? CREATE_ALWAYS : OPEN_EXISTING,
      FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
  if (file != INVALID_HANDLE_VALUE) {
    DWORD written = 0;
    if (trace_invocation != 0)
      SetFilePointer(file, 0, (long *)0, FILE_END);
    WriteFile(file, record, (DWORD)sizeof(*record), &written, (void *)0);
    CloseHandle(file);
  }
  ++trace_invocation;
}

int STDCALL FindSideCapture(
    int top_left, int top_right, int sequence_length, int left, int right,
    int window, int informative, int comparison, char *scores,
    double *high_left, double *high_right) {
  static FindSideFn original;
  struct TraceRecord record = {0};
  int result;
  if (!original) {
    HMODULE module = LoadLibraryA("dna_original.dll");
    if (module) original = (FindSideFn)GetProcAddress(module, "FindSide");
  }
  if (!original) return 0;
  result = original(
      top_left, top_right, sequence_length, left, right, window, informative,
      comparison, scores, high_left, high_right);
  record.function = 1;
  record.values[0] = top_left;
  record.values[1] = top_right;
  record.values[2] = sequence_length;
  record.values[3] = left;
  record.values[4] = right;
  record.values[5] = window;
  record.values[6] = informative;
  record.values[7] = comparison;
  record.values[8] = result;
  record.doubles[0] = *high_left;
  record.doubles[1] = *high_right;
  append_record(&record);
  return result;
}

static int capture_opt(
    int function, const char *name, OptFn *original, int edge,
    double high, int top, int maximum, int comparison, int window,
    int informative, int sequence_length, char *scores, char *missing_map) {
  struct TraceRecord record = {0};
  int result;
  if (!*original) {
    HMODULE module = LoadLibraryA("dna_original.dll");
    if (module) *original = (OptFn)GetProcAddress(module, name);
  }
  if (!*original) return 0;
  result = (*original)(
      edge, high, top, maximum, comparison, window, informative,
      sequence_length, scores, missing_map);
  record.function = function;
  record.values[0] = edge;
  record.values[1] = top;
  record.values[2] = maximum;
  record.values[3] = comparison;
  record.values[4] = window;
  record.values[5] = informative;
  record.values[6] = sequence_length;
  record.values[7] = result;
  record.doubles[0] = high;
  append_record(&record);
  return result;
}

int STDCALL OptLeftBPMCCapture(
    int left, double high, int top, int maximum, int comparison, int window,
    int informative, int sequence_length, char *scores, char *missing_map) {
  static OptFn original;
  return capture_opt(
      2, "OptLeftBPMC", &original, left, high, top, maximum, comparison,
      window, informative, sequence_length, scores, missing_map);
}

int STDCALL OptRightBPMCCapture(
    int right, double high, int top, int maximum, int comparison, int window,
    int informative, int sequence_length, char *scores, char *missing_map) {
  static OptFn original;
  return capture_opt(
      3, "OptRightBPMC", &original, right, high, top, maximum, comparison,
      window, informative, sequence_length, scores, missing_map);
}

int STDCALL DestroyPeaksCapture(
    int comparison, int informative, int sequence_length, int left, int right,
    double *smooth, double *chi_values) {
  static DestroyPeaksFn original;
  struct TraceRecord record = {0};
  int result;
  if (!original) {
    HMODULE module = LoadLibraryA("dna_original.dll");
    if (module)
      original = (DestroyPeaksFn)GetProcAddress(module, "DestroyPeaks");
  }
  if (!original) return 0;
  result = original(
      comparison, informative, sequence_length, left, right, smooth,
      chi_values);
  record.function = 4;
  record.values[0] = comparison;
  record.values[1] = informative;
  record.values[2] = sequence_length;
  record.values[3] = left;
  record.values[4] = right;
  record.values[5] = result;
  append_record(&record);
  return result;
}
