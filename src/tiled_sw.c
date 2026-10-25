#if !defined(_WIN32) && !defined(_WIN64)
#define _POSIX_C_SOURCE 200112L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#define WINDOWS 1
#include <windows.h>
#include <malloc.h>
#else
#define WINDOWS 0
#include <time.h>
#include <sys/time.h>
#endif

static void *aligned_allocate(size_t size, size_t alignment) {
#if WINDOWS
    return _aligned_malloc(size, alignment);
#else
    void *ptr = NULL;
    if (posix_memalign(&ptr, alignment, size) != 0) return NULL;
    return ptr;
#endif
}

static void aligned_release(void *ptr) {
#if WINDOWS
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

static int parse_positive_int(const char *text, int *value) {
    char *end = NULL;
    long parsed = strtol(text, &end, 10);
    if (text == end || *end != '\0' || parsed < 1 || parsed > INT_MAX)
        return 0;
    *value = (int)parsed;
    return 1;
}

// --------------------------------------------------------
// High-resolution timer
// --------------------------------------------------------
static double walltime_seconds(void) {
#if WINDOWS
    static double freq = 0.0;
    if (freq == 0.0) {
        LARGE_INTEGER f;
        QueryPerformanceFrequency(&f);
        freq = (double)f.QuadPart;
    }
    LARGE_INTEGER t;
    QueryPerformanceCounter(&t);
    return (double)t.QuadPart / freq;
#else
#if defined(CLOCK_MONOTONIC)
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
        return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
#endif
}

// --------------------------------------------------------
// FASTA reader
// --------------------------------------------------------
typedef struct { char *a; char *b; } Pair;

static char *read_file_all(const char *path, size_t *out_len) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;

    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    rewind(f);

    char *buf = (char*)malloc((size_t)sz + 1);
    if (!buf) { fclose(f); return NULL; }

    size_t r = fread(buf, 1, (size_t)sz, f);
    buf[r] = '\0';

    if (out_len) *out_len = r;
    fclose(f);
    return buf;
}

static Pair *read_pairs(const char *fname, int *out_n) {
    *out_n = 0;

    size_t len;
    char *buf = read_file_all(fname, &len);
    if (!buf) return NULL;

    // Normalize newlines (remove \r)
    char *wp = buf;
    for (char *p = buf; *p; ++p)
        if (*p != '\r') *wp++ = *p;
    *wp = '\0';

    char *p = buf;
    char *seqs[65536];
    int seqc = 0;

    while (*p) {
        while (*p && *p != '>') ++p;
        if (!*p) break;
        ++p;

        while (*p && *p != '\n') ++p;
        if (*p == '\n') ++p;

        size_t cap = 256;
        char *acc = (char*)malloc(cap);
        size_t acc_len = 0;

        while (*p && *p != '>') {
            char *nl = strchr(p, '\n');
            size_t chunk = nl ? (size_t)(nl - p) : strlen(p);

            for (size_t i = 0; i < chunk; ++i) {
                char c = p[i];
                if (c == ' ' || c == '\t' || c == '\n' || c == '\r') continue;
                if (c >= 'a' && c <= 'z') c -= 32; // upper-case

                if (acc_len + 2 >= cap) {
                    cap *= 2;
                    acc = (char*)realloc(acc, cap);
                    if (!acc) { free(buf); return NULL; }
                }
                acc[acc_len++] = c;
            }

            if (!nl) { p += chunk; break; }
            p = nl + 1;
        }

        acc[acc_len] = '\0';
        seqs[seqc++] = acc;
    }

    if (seqc < 2) {
        for (int i = 0; i < seqc; ++i) free(seqs[i]);
        free(buf);
        return NULL;
    }

    int pairs = seqc / 2;
    Pair *out = (Pair*)malloc(sizeof(Pair) * pairs);
    if (!out) {
        for (int i = 0; i < seqc; ++i) free(seqs[i]);
        free(buf);
        return NULL;
    }

    for (int i = 0; i < pairs; ++i) {
        out[i].a = seqs[2 * i];
        out[i].b = seqs[2 * i + 1];
    }

    *out_n = pairs;
    free(buf);
    return out;
}

// --------------------------------------------------------
// Smith–Waterman core (int16, vectorization-friendly)
// --------------------------------------------------------
static int smith_waterman_score_fast(
    const char * restrict A,
    const char * restrict B,
    int16_t match,
    int16_t mismatch,
    int16_t gap,
    int16_t * restrict prev,
    int16_t * restrict cur,
    int m
) {
    int n = (int)strlen(A);
    if (n == 0 || m == 0) return 0;

#pragma omp simd
    for (int j = 0; j <= m; ++j)
        prev[j] = 0;

    int16_t best = 0;

    for (int i = 1; i <= n; ++i) {
        cur[0] = 0;
        char ai = A[i - 1];

#pragma omp simd
        for (int j = 1; j <= m; ++j) {
            int16_t diag = prev[j - 1] + (ai == B[j - 1] ? match : mismatch);
            int16_t up   = prev[j] + gap;
            int16_t left = cur[j - 1] + gap;

            int16_t v = diag;
            if (up > v) v = up;
            if (left > v) v = left;
            if (v < 0) v = 0;

            cur[j] = v;
            if (v > best) best = v;
        }

        int16_t *tmp = prev;
        prev = cur;
        cur  = tmp;
    }

    return (int)best;
}

// --------------------------------------------------------
// MAIN: one harness-controlled benchmark configuration
// --------------------------------------------------------
int main(int argc, char **argv) {
    int threads = 1;
    int chunk_size = 1;
    if (argc < 2 || argc > 4 ||
        (argc >= 3 && !parse_positive_int(argv[2], &threads)) ||
        (argc >= 4 && !parse_positive_int(argv[3], &chunk_size))) {
        fprintf(stderr, "Usage: %s pairs.fa [threads [schedule_chunk]]\n", argv[0]);
        return 1;
    }

#ifndef _OPENMP
    if (threads > 1) {
        fprintf(stderr, "This binary was built without OpenMP; rebuild with -fopenmp\n");
        return 1;
    }
#else
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#endif

    const char *pairs_file = argv[1];

    int npairs = 0;
    Pair *pairs = read_pairs(pairs_file, &npairs);
    if (!pairs || npairs == 0) {
        fprintf(stderr, "Error: cannot open %s or no pairs found\n", pairs_file);
        return 1;
    }

    int max_m = 0;
    for (int i = 0; i < npairs; ++i) {
        int lb = (int)strlen(pairs[i].b);
        if (lb > max_m) max_m = lb;
    }
    if (max_m <= 0) {
        fprintf(stderr, "No valid sequences in %s\n", pairs_file);
        return 1;
    }

    size_t row_capacity = ((size_t)max_m + 1 + 15) / 16 * 16;
    long long score_sum = 0;
    int allocation_failed = 0;
    double t0 = walltime_seconds();

#ifdef _OPENMP
#pragma omp parallel reduction(+:score_sum) reduction(|:allocation_failed)
        {
            int16_t *prev = (int16_t*)aligned_allocate(row_capacity * sizeof(int16_t), 64);
            int16_t *cur  = (int16_t*)aligned_allocate(row_capacity * sizeof(int16_t), 64);
            if (!prev || !cur) allocation_failed = 1;

#pragma omp for schedule(dynamic,chunk_size)
            for (int i = 0; i < npairs; ++i) {
                if (prev && cur) {
                    int m = (int)strlen(pairs[i].b);
                    int sc = smith_waterman_score_fast(
                        pairs[i].a, pairs[i].b,
                        2, -2, -1,
                        prev, cur, m
                    );
                    score_sum += sc;
                }
            }

            aligned_release(prev);
            aligned_release(cur);
        }
#else
    int16_t *prev = (int16_t*)aligned_allocate(row_capacity * sizeof(int16_t), 64);
    int16_t *cur = (int16_t*)aligned_allocate(row_capacity * sizeof(int16_t), 64);
    if (!prev || !cur) {
        fprintf(stderr, "Error: unable to allocate DP rows\n");
        aligned_release(prev);
        aligned_release(cur);
        return 1;
    }
    for (int i = 0; i < npairs; ++i) {
        int m = (int)strlen(pairs[i].b);
        int sc = smith_waterman_score_fast(
            pairs[i].a, pairs[i].b,
            2, -2, -1,
            prev, cur, m
        );
        score_sum += sc;
    }

    aligned_release(prev);
    aligned_release(cur);
#endif

    if (allocation_failed) {
        fprintf(stderr, "Error: unable to allocate DP rows for all OpenMP threads\n");
        return 1;
    }

    double elapsed = walltime_seconds() - t0;
    printf("threads=%d  chunk=%d  pairs=%d  time=%.6f  score_sum=%lld\n",
           threads, chunk_size, npairs, elapsed, score_sum);

    // cleanup
    for (int i = 0; i < npairs; ++i) {
        free(pairs[i].a);
        free(pairs[i].b);
    }
    free(pairs);

    return 0;
}
