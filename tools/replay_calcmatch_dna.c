#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

typedef int(__stdcall *MakeVarMap2Fn)(
    int, int, int, short *, short *, int *, int *, int *, int *);
typedef int(__stdcall *MakeCntHit2Fn)(
    int, int, int, int, int, int, int, double *, short *, float *, int *);

static void read_exact(FILE *file, void *data, size_t bytes) {
  if (fread(data, 1, bytes, file) != bytes) {
    fprintf(stderr, "truncated fixture\n");
    exit(2);
  }
}

static size_t first_difference(
    const void *actual, const void *expected, size_t bytes) {
  const unsigned char *a = (const unsigned char *)actual;
  const unsigned char *e = (const unsigned char *)expected;
  size_t index;
  for (index = 0; index < bytes; ++index) {
    if (a[index] != e[index]) return index;
  }
  return bytes;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: replay_calcmatch_dna fixture dna.dll\n");
    return 2;
  }
  FILE *file = fopen(argv[1], "rb");
  if (!file) {
    fprintf(stderr, "cannot open fixture\n");
    return 2;
  }
  unsigned int header[10];
  int sequences[3];
  int comparison[6];
  read_exact(file, header, sizeof(header));
  read_exact(file, sequences, sizeof(sequences));
  read_exact(file, comparison, sizeof(comparison));
  if (header[0] != 0x594d4343U) {
    fprintf(stderr, "bad fixture magic\n");
    return 2;
  }
  const int next_no = (int)header[3];
  const int count = next_no + 1;
  const int sequence_length = (int)header[6];
  const int map_upper_bound = 160;
  const size_t sequence_cells = (size_t)(sequence_length + 1) * count;
  const size_t map_cells = (size_t)3 * (map_upper_bound + 1) * count;
  const size_t position_cells = (size_t)sequence_length + 2;
  short *sequence_data = (short *)calloc(sequence_cells, sizeof(short));
  short *expected_map = (short *)calloc(map_cells, sizeof(short));
  int *expected_regions = (int *)calloc(position_cells, sizeof(int));
  int *expected_positions = (int *)calloc(position_cells, sizeof(int));
  read_exact(file, sequence_data, sequence_cells * sizeof(short));
  int expected_variable_sites;
  read_exact(file, &expected_variable_sites, sizeof(expected_variable_sites));
  read_exact(file, expected_map, map_cells * sizeof(short));
  read_exact(file, expected_regions, position_cells * sizeof(int));
  read_exact(
      file, expected_positions,
      (size_t)(sequence_length + 1) * sizeof(int));
  double *expected_hits = (double *)calloc((size_t)6 * count, sizeof(double));
  read_exact(file, expected_hits, (size_t)6 * count * sizeof(double));
  const size_t smooth_cells =
      (size_t)3 * (expected_variable_sites + 1) * count;
  float *expected_smooth = (float *)calloc(smooth_cells, sizeof(float));
  read_exact(file, expected_smooth, smooth_cells * sizeof(float));
  fclose(file);

  HMODULE module = LoadLibraryA(argv[2]);
  if (!module) {
    fprintf(stderr, "cannot load dna.dll: %lu\n", GetLastError());
    return 2;
  }
  MakeVarMap2Fn make_var_map =
      (MakeVarMap2Fn)GetProcAddress(module, "MakeVarMap2");
  MakeCntHit2Fn make_count_hits =
      (MakeCntHit2Fn)GetProcAddress(module, "MakeCntHit2");
  if (!make_var_map || !make_count_hits) {
    fprintf(stderr, "missing exports\n");
    return 2;
  }

  short *actual_map = (short *)calloc(map_cells, sizeof(short));
  int *actual_regions = (int *)calloc(position_cells, sizeof(int));
  int *actual_positions = (int *)calloc(position_cells, sizeof(int));
  actual_regions[sequence_length + 1] = sequence_length + 1;
  const int actual_variable_sites = make_var_map(
      next_no, sequence_length, map_upper_bound, sequence_data, actual_map,
      actual_regions, actual_positions, sequences, comparison);
  const size_t map_difference = first_difference(
      actual_map, expected_map, map_cells * sizeof(short));
  const size_t region_difference = first_difference(
      actual_regions, expected_regions,
      (size_t)(sequence_length + 1) * sizeof(int));
  const size_t position_difference = first_difference(
      actual_positions, expected_positions,
      (size_t)(sequence_length + 1) * sizeof(int));
  printf(
      "MakeVarMap2 result=%s map=%s regions=%s positions=%s\n",
      actual_variable_sites == expected_variable_sites ? "PASS" : "FAIL",
      map_difference == map_cells * sizeof(short) ? "PASS" : "FAIL",
      region_difference == (size_t)(sequence_length + 1) * sizeof(int)
          ? "PASS" : "FAIL",
      position_difference == (size_t)(sequence_length + 1) * sizeof(int)
          ? "PASS" : "FAIL");
  if (map_difference != map_cells * sizeof(short)) {
    printf("map-first-byte-difference=%u\n", (unsigned int)map_difference);
  }

  double *actual_hits = (double *)calloc((size_t)6 * count, sizeof(double));
  float *actual_smooth = (float *)calloc(smooth_cells, sizeof(float));
  const int count_result = make_count_hits(
      (int)header[7], (int)header[8], (int)header[9], next_no,
      actual_variable_sites, sequence_length, map_upper_bound, actual_hits,
      actual_map, actual_smooth, actual_regions);
  const size_t hit_difference = first_difference(
      actual_hits, expected_hits, (size_t)6 * count * sizeof(double));
  const size_t smooth_difference = first_difference(
      actual_smooth, expected_smooth, smooth_cells * sizeof(float));
  printf(
      "MakeCntHit2 result=%d hits=%s smooth=%s first-hit-byte=%u "
      "first-smooth-byte=%u\n",
      count_result,
      hit_difference == (size_t)6 * count * sizeof(double) ? "PASS" : "FAIL",
      smooth_difference == smooth_cells * sizeof(float) ? "PASS" : "FAIL",
      (unsigned int)hit_difference, (unsigned int)smooth_difference);
  printf(
      "role1-seq4 outside=%.17g/%.17g inside=%.17g/%.17g "
      "product=%.17g/%.17g\n",
      actual_hits[1 + 4 * 6], expected_hits[1 + 4 * 6],
      actual_hits[4 + 4 * 6], expected_hits[4 + 4 * 6],
      actual_hits[1 + 4 * 6] * actual_hits[4 + 4 * 6],
      expected_hits[1 + 4 * 6] * expected_hits[4 + 4 * 6]);
  return 0;
}
