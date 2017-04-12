/* -*- Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <unistd.h>
#include <sys/time.h>

#include "bitmap.h"
#include "shared.h"

double densityL = 0.1;
double densityR = 0.1;

uint64_t numOnesL = 0;
uint64_t numOnesR = 0;

const int numIters = 10;
const int numQueries = 1000000;
uint64_t queries[numQueries];

uint32_t seed = 1;

inline uint32_t xRand()
{
    return seed = (279470273ULL * seed) % 4294967291ULL;
}

inline uint64_t xRand64()
{
    return (uint64_t) xRand() << 32 | xRand();
}

uint64_t* createRandomBits(uint64_t nbits, uint32_t thresholdL, uint32_t thresholdR)
{
    fprintf(stderr, "nbits to create: %" PRIu64 "\n", nbits);
    fprintf(stderr, "allocated bits: %" PRIu64 " bytes\n", nbits/8);

    uint64_t* bits = NULL;
    assert(posix_memalign((void **) &bits, 4096, nbits / 8) == 0);

    for (uint64_t i = 0; i < nbits / 2; i++) {
        if (xRand() < thresholdL) { 
            bits[i / 64] |= 1UL << (i % 64);
            ++numOnesL;
        } else {
            bits[i / 64] &= ~(1ULL << (i % 64));
        }
    }
    for (uint64_t i = nbits / 2; i < nbits; i++) {
        if (xRand() < thresholdR) {
            bits[i / 64] |= 1ULL << (i % 64);
            ++numOnesR;
        } else {
            bits[i / 64] &= ~(1ULL << (i % 64));
        }
    }

    return bits;
}

enum benchmode {
    BENCH_RANK,
    BENCH_SELECT,
};

int main(int argc, char **argv)
{
    extern int optind;
    int ch;

    uint64_t nbits;
    benchmode mode = BENCH_RANK;

    while ((ch = getopt(argc, argv, "sn:d:")) != -1) {
        switch (ch) {
        case 's':
            mode = BENCH_SELECT;
            break;
        case 'n':
            nbits = atoi(optarg);
            nbits = 1ULL << nbits;
            break;
        case 'd':
            densityL = densityR = atof(optarg);
            break;
        }
    }

    printf("benchmode: %s\n", mode == BENCH_RANK ? "rank" : "select");

    uint32_t thresholdL = (uint32_t) (UINT32_MAX * densityL);
    uint32_t thresholdR = (uint32_t) (UINT32_MAX * densityR);

    uint64_t* bits = createRandomBits(nbits, thresholdL, thresholdR);
    BitmapPoppy* bitmap = new BitmapPoppy(bits, nbits);
    uint64_t dummy = 0x1234567890ABCDEF;

    if (mode == BENCH_RANK) {
        for (int i = 0; i < numQueries; i++) {
            queries[i] = xRand64() % nbits + 1;
        }
    } else {
        assert(mode == BENCH_SELECT);

        for (int i = 0; i < numQueries / 2; i++) {
            queries[i] = xRand64() % numOnesL + 1;
        }
        for (int i = numQueries / 2; i < numQueries; i++) {
            queries[i] = xRand64() % numOnesR + 1 + numOnesL;
        }
    }

    struct timeval tv_start, tv_end;
    gettimeofday(&tv_start, NULL);

    if (mode == BENCH_RANK) {
        for (int iter = 0; iter < numIters; iter++)
            for (int i = 0; i < numQueries; i++)
                dummy ^= bitmap->rank(queries[i]);
    } else {
        assert(mode == BENCH_SELECT);

        for (int iter = 0; iter < numIters; iter++)
            for (int i = 0; i < numQueries; i++)
                dummy ^= bitmap->select(queries[i]);
    }
    gettimeofday(&tv_end, NULL);

    double elapsed_seconds = timeval_diff(&tv_start, &tv_end);
    printf("%" PRIu64 " ops, %.2f seconds, ns/op: %.2f\n", 
           (uint64_t) numIters * numQueries, 
           elapsed_seconds,
           elapsed_seconds * 1000000000 / ((uint64_t) numIters * numQueries));

    if (dummy == 42) printf("42\n");

    return 0;
}
