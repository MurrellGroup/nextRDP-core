#include <windows.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#pragma pack(push, 1)
struct MakeSDMP2Header {
  char magic[8];
  uint32_t version;
  int32_t next_no;
  int32_t sequence_length;
};

struct FillRmatHeader {
  char magic[8];
  uint32_t version;
  int32_t y;
  int32_t next_no;
  int32_t result_matrix_ub1;
  int32_t result_matrix_ub2;
  int32_t distance_matrix_ub1;
  int32_t distance_matrix_ub2;
  int32_t distance_matrix_ub3;
};

struct CalCRChainHeader {
  char magic[8];
  uint32_t version;
  int32_t next_no;
};

struct SectionHeader {
  uint32_t id;
  uint32_t bytes;
};
#pragma pack(pop)

typedef double(__stdcall *CalCRFn)(
    double, int, int, int *, int *, float *, float *, double *, double *,
    float *);

static void *read_section(const char *path, uint32_t wanted, uint32_t *bytes) {
  FILE *file = fopen(path, "rb");
  struct SectionHeader section;
  void *data;
  if (!file) return NULL;
  if (fseek(file, 0, SEEK_SET) != 0) {
    fclose(file);
    return NULL;
  }
  /* Every source fixture used here has a 20- or 40-byte packed header. */
  {
    char magic[8];
    fread(magic, 1, 8, file);
    if (memcmp(magic, "MAKESDMP", 8) == 0) {
      fseek(file, (long)sizeof(struct MakeSDMP2Header), SEEK_SET);
    } else if (memcmp(magic, "FILLRMAT", 8) == 0) {
      fseek(file, (long)sizeof(struct FillRmatHeader), SEEK_SET);
    } else {
      fclose(file);
      return NULL;
    }
  }
  while (fread(&section, sizeof(section), 1, file) == 1) {
    if (section.id == 0xffffffffU) break;
    if (section.id == wanted) {
      data = calloc(section.bytes ? section.bytes : 1, 1);
      if (!data || fread(data, 1, section.bytes, file) != section.bytes) {
        free(data);
        fclose(file);
        return NULL;
      }
      fclose(file);
      *bytes = section.bytes;
      return data;
    }
    fseek(file, (long)section.bytes, SEEK_CUR);
  }
  fclose(file);
  return NULL;
}

static int write_section(
    FILE *file, uint32_t id, const void *data, uint32_t bytes) {
  const struct SectionHeader section = {id, bytes};
  return fwrite(&section, sizeof(section), 1, file) == 1 &&
      fwrite(data, 1, bytes, file) == bytes;
}

int main(int argc, char **argv) {
  struct MakeSDMP2Header source_header;
  struct CalCRChainHeader output_header;
  struct SectionHeader end = {0xffffffffU, 0};
  uint32_t bytes;
  int *sequences;
  int *comparison;
  double *rmat[3];
  float *rcorr;
  float *rinv;
  float *trcorr;
  double int_value[2] = {0.0, 0.0};
  double results[3];
  FILE *file;
  HMODULE library;
  CalCRFn function;
  int y;
  size_t sequence_count;
  if (argc != 6) {
    fprintf(stderr, "usage: capture-calcr-chain make-sdmp fill0 fill1 fill2 output\n");
    return 2;
  }
  file = fopen(argv[1], "rb");
  if (!file || fread(&source_header, sizeof(source_header), 1, file) != 1) {
    if (file) fclose(file);
    fprintf(stderr, "cannot read MakeSDMP2 fixture\n");
    return 1;
  }
  fclose(file);
  sequences = (int *)read_section(argv[1], 3, &bytes);
  comparison = (int *)read_section(argv[1], 4, &bytes);
  for (y = 0; y < 3; ++y) {
    rmat[y] = (double *)read_section(argv[2 + y], 101, &bytes);
  }
  if (!sequences || !comparison || !rmat[0] || !rmat[1] || !rmat[2]) {
    fprintf(stderr, "cannot read source fixture sections\n");
    return 1;
  }
  sequence_count = (size_t)(source_header.next_no + 1);
  rcorr = (float *)calloc(9 * sequence_count, sizeof(float));
  rinv = (float *)calloc(9 * sequence_count, sizeof(float));
  trcorr = (float *)calloc(45 * sequence_count, sizeof(float));
  if (!rcorr || !rinv || !trcorr) return 1;
  library = LoadLibraryA("DNA.dll");
  function = library ? (CalCRFn)GetProcAddress(library, "CalCR") : NULL;
  if (!function) {
    fprintf(stderr, "cannot load DNA.dll!CalCR\n");
    return 1;
  }
  for (y = 0; y < 3; ++y) {
    int_value[0] = 0.0;
    int_value[1] = 0.0;
    results[y] = function(
        1.0e-14, source_header.next_no, y, sequences, comparison, rcorr,
        rinv, int_value, rmat[y], trcorr);
  }
  memcpy(output_header.magic, "CALCR3\0\0", 8);
  output_header.version = 1;
  output_header.next_no = source_header.next_no;
  file = fopen(argv[5], "wb");
  if (!file) return 1;
  fwrite(&output_header, sizeof(output_header), 1, file);
  write_section(file, 1, sequences, 3U * (uint32_t)sizeof(int));
  write_section(file, 2, comparison, 6U * (uint32_t)sizeof(int));
  for (y = 0; y < 3; ++y) {
    write_section(
        file, (uint32_t)(3 + y), rmat[y],
        (uint32_t)(18 * sequence_count * sizeof(double)));
  }
  write_section(
      file, 101, rcorr, (uint32_t)(9 * sequence_count * sizeof(float)));
  write_section(
      file, 102, rinv, (uint32_t)(9 * sequence_count * sizeof(float)));
  write_section(
      file, 103, trcorr, (uint32_t)(45 * sequence_count * sizeof(float)));
  write_section(file, 104, int_value, (uint32_t)sizeof(int_value));
  write_section(file, 105, results, (uint32_t)sizeof(results));
  fwrite(&end, sizeof(end), 1, file);
  fclose(file);
  FreeLibrary(library);
  return 0;
}

