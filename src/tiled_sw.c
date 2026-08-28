#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <malloc.h>   // _aligned_malloc, _aligned_free

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#define WINDOWS 1
#include <windows.h>
#else
#define WINDOWS 0
#include <time.h>
#include <sys/time.h>
#endif

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
// MAIN: automatic thread benchmarking
// --------------------------------------------------------
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s pairs.fa\n", argv[0]);
        return 1;
    }

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

    int thread_list[] = {1, 2, 4, 8, 16};
    int nconfigs = sizeof(thread_list) / sizeof(thread_list[0]);

    printf("=== Tiled Smith-Waterman Auto Benchmark ===\n");
    printf("Dataset: %s  (pairs=%d, max_len_B=%d)\n\n", pairs_file, npairs, max_m);

    double baseline_time = 0.0;

    for (int c = 0; c < nconfigs; ++c) {
        int threads = thread_list[c];

#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif

        long long score_sum = 0;
        double t0 = walltime_seconds();

#ifdef _OPENMP
#pragma omp parallel reduction(+:score_sum)
        {
            int16_t *prev = (int16_t*)_aligned_malloc((size_t)(max_m + 1) * sizeof(int16_t), 64);
            int16_t *cur  = (int16_t*)_aligned_malloc((size_t)(max_m + 1) * sizeof(int16_t), 64);

#pragma omp for schedule(dynamic,1)
            for (int i = 0; i < npairs; ++i) {
                int m = (int)strlen(pairs[i].b);
                int sc = smith_waterman_score_fast(
                    pairs[i].a, pairs[i].b,
                    2, -2, -1,
                    prev, cur, m
                );
                score_sum += sc;
            }

            _aligned_free(prev);
            _aligned_free(cur);
        }
#else
        int16_t *prev = (int16_t*)_aligned_malloc((size_t)(max_m + 1) * sizeof(int16_t), 64);
        int16_t *cur  = (int16_t*)_aligned_malloc((size_t)(max_m + 1) * sizeof(int16_t), 64);

        for (int i = 0; i < npairs; ++i) {
            int m = (int)strlen(pairs[i].b);
            int sc = smith_waterman_score_fast(
                pairs[i].a, pairs[i].b,
                2, -2, -1,
                prev, cur, m
            );
            score_sum += sc;
        }

        _aligned_free(prev);
        _aligned_free(cur);
#endif

        double t1 = walltime_seconds();
        double elapsed = t1 - t0;

        if (c == 0) baseline_time = elapsed;
        double speedup = (baseline_time > 0.0) ? (baseline_time / elapsed) : 1.0;

        printf("threads=%2d  pairs=%d  time=%.6f  score_sum=%lld  speedup=%.3f\n",
               threads, npairs, elapsed, score_sum, speedup);
    }

    // cleanup
    for (int i = 0; i < npairs; ++i) {
        free(pairs[i].a);
        free(pairs[i].b);
    }
    free(pairs);

    return 0;
}
