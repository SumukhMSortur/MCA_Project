#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <malloc.h>

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

#define NEG_INF_16 ((int16_t)-30000)

// --------------------------------------------------------
// Timer
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
        return ts.tv_sec + ts.tv_nsec*1e-9;
#endif
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec + tv.tv_usec*1e-6;
#endif
}

// --------------------------------------------------------
// FASTA reader
// --------------------------------------------------------
typedef struct { char *a; char *b; } Pair;

static char* read_file_all(const char *path,size_t *out_len){
    FILE *f=fopen(path,"rb");
    if(!f) return NULL;
    fseek(f,0,SEEK_END);
    long sz=ftell(f);
    rewind(f);
    char *buf=(char*)malloc(sz+1);
    if(!buf){ fclose(f); return NULL; }
    size_t r=fread(buf,1,sz,f);
    buf[r]='\0';
    fclose(f);
    *out_len=r;
    return buf;
}

static Pair* read_pairs(const char *fname,int *out_n){
    *out_n=0;
    size_t l;
    char *buf=read_file_all(fname,&l);
    if(!buf) return NULL;

    char *wp=buf;
    for(char *p=buf;*p;p++){
        if(*p!='\r') *wp++=*p;
    }
    *wp='\0';

    char *p=buf;
    char *seqs[65536];
    int seqc=0;

    while(*p){
        while(*p && *p!='>') p++;
        if(!*p) break;
        p++;
        while(*p && *p!='\n') p++;
        if(*p=='\n') p++;

        size_t cap=256;
        char *acc=(char*)malloc(cap);
        size_t acc_len=0;

        while(*p && *p!='>'){
            char *nl=strchr(p,'\n');
            size_t chunk=nl?(size_t)(nl-p):strlen(p);

            for(size_t i=0;i<chunk;i++){
                char c=p[i];
                if(c==' '||c=='\n'||c=='\t'||c=='\r') continue;
                if(c>='a' && c<='z') c-=32;
                if(acc_len+2>=cap){
                    cap*=2;
                    acc=(char*)realloc(acc,cap);
                    if(!acc){ free(buf); return NULL; }
                }
                acc[acc_len++]=c;
            }

            if(!nl){ p+=chunk; break; }
            p=nl+1;
        }
        acc[acc_len]='\0';
        seqs[seqc++]=acc;
    }

    if(seqc<2){
        for(int i=0;i<seqc;i++) free(seqs[i]);
        free(buf);
        return NULL;
    }

    int pairs=seqc/2;
    Pair *out=(Pair*)malloc(sizeof(Pair)*pairs);
    if(!out){
        for(int i=0;i<seqc;i++) free(seqs[i]);
        free(buf);
        return NULL;
    }

    for(int i=0;i<pairs;i++){
        out[i].a=seqs[2*i];
        out[i].b=seqs[2*i+1];
    }

    *out_n=pairs;
    free(buf);
    return out;
}

// --------------------------------------------------------
// Banded DP core
// --------------------------------------------------------
static int banded_score_fast(
    const char * restrict A,
    const char * restrict B,
    int band,
    int16_t match,
    int16_t mismatch,
    int16_t gap,
    int16_t * restrict prev,
    int16_t * restrict cur,
    int W
){
    int n = (int)strlen(A);
    int m = (int)strlen(B);

#pragma omp simd
    for(int k=0;k<W;k++)
        prev[k]=NEG_INF_16;

    int16_t best = 0;

    for(int i=0;i<=n;i++){
        int jmin = i - band;
        if(jmin<0) jmin=0;
        int jmax = i + band;
        if(jmax>m) jmax=m;

#pragma omp simd
        for(int k=0;k<W;k++)
            cur[k]=NEG_INF_16;

        for(int jj=jmin;jj<=jmax;jj++){
            int idx = jj - (i - band);
            if(idx<0 || idx>=W) continue;

            int16_t score = NEG_INF_16;

            if(i==0 && jj==0){
                score = 0;
            } else {
                if(i>0 && jj>0){
                    int pd = (jj-1) - ((i-1)-band);
                    if(pd>=0 && pd<W){
                        int16_t d = prev[pd] +
                                    (A[i-1]==B[jj-1] ? match : mismatch);
                        if(d>score) score=d;
                    }
                }
                if(i>0){
                    int pu = jj - ((i-1)-band);
                    if(pu>=0 && pu<W){
                        int16_t u = prev[pu] + gap;
                        if(u>score) score=u;
                    }
                }
                if(jj>0){
                    int pl = (jj-1) - (i-band);
                    if(pl>=0 && pl<W){
                        int16_t l = cur[pl] + gap;
                        if(l>score) score=l;
                    }
                }
            }

            cur[idx] = score;
            if(score>best) best=score;
        }

        int16_t *tmp=prev; prev=cur; cur=tmp;
    }

    if(best<0) best=0;
    return (int)best;
}

// --------------------------------------------------------
// MAIN: automatic thread benchmark
// --------------------------------------------------------
int main(int argc,char **argv){
    if(argc<3){
        fprintf(stderr,"Usage: %s pairs.fa <band>\n",argv[0]);
        return 1;
    }

    const char *fname = argv[1];
    int band = atoi(argv[2]);
    if(band<0) band=0;

    int npairs = 0;
    Pair *pairs = read_pairs(fname,&npairs);
    if(!pairs || npairs==0){
        fprintf(stderr,"Error: cannot open %s or no pairs found\n",fname);
        return 1;
    }

    int W = 2*band + 1;
    if(W <= 0){
        fprintf(stderr,"Band too small\n");
        return 1;
    }

    int thread_list[] = {1,2,4,8,16};
    int nconfigs = sizeof(thread_list)/sizeof(thread_list[0]);

    printf("=== WFA/Banded DP Auto Benchmark ===\n");
    printf("Dataset: %s  (pairs=%d, band=%d, W=%d)\n\n", fname, npairs, band, W);

    double baseline_time = 0.0;

    for(int c=0;c<nconfigs;c++){
        int threads = thread_list[c];

#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif

        long long score_sum = 0;
        double t0 = walltime_seconds();

#ifdef _OPENMP
#pragma omp parallel reduction(+:score_sum)
        {
            int16_t *prev = (int16_t*)_aligned_malloc((size_t)W*sizeof(int16_t),64);
            int16_t *cur  = (int16_t*)_aligned_malloc((size_t)W*sizeof(int16_t),64);

#pragma omp for schedule(dynamic,1)
            for(int i=0;i<npairs;i++){
                int sc = banded_score_fast(
                    pairs[i].a, pairs[i].b,
                    band,
                    2, -2, -1,
                    prev, cur, W
                );
                score_sum += sc;
            }

            _aligned_free(prev);
            _aligned_free(cur);
        }
#else
        int16_t *prev = (int16_t*)_aligned_malloc((size_t)W*sizeof(int16_t),64);
        int16_t *cur  = (int16_t*)_aligned_malloc((size_t)W*sizeof(int16_t),64);

        for(int i=0;i<npairs;i++){
            int sc = banded_score_fast(
                pairs[i].a, pairs[i].b,
                band,
                2, -2, -1,
                prev, cur, W
            );
            score_sum += sc;
        }

        _aligned_free(prev);
        _aligned_free(cur);
#endif

        double t1 = walltime_seconds();
        double elapsed = t1 - t0;

        if(c==0) baseline_time = elapsed;
        double speedup = (baseline_time>0.0) ? (baseline_time/elapsed) : 1.0;

        printf("threads=%2d  pairs=%d  time=%.6f  score_sum=%lld  speedup=%.3f\n",
               threads, npairs, elapsed, score_sum, speedup);
    }

    for(int i=0;i<npairs;i++){
        free(pairs[i].a);
        free(pairs[i].b);
    }
    free(pairs);

    return 0;
}
