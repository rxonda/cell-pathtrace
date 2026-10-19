#include <stdio.h>
#include "spu-rng.h"

int main(void) {
    static const unsigned int expected[3][4] = {
        {1793u, 2050u, 2307u, 2564u},
        {2098313u, 2098442u, 2098571u, 6295068u},
        {5521925u, 5784069u, 6046213u, 6344223u}
    };
    static const unsigned int seeded[4] = {
        1713311875u, 2366773723u, 3812440945u, 1511383879u
    };
    static const unsigned int next_pixel[4] = {
        399893900u, 207984398u, 2137489156u, 3986754366u
    };
    SpuRng rng;
    int step, lane;
    rng.s0 = (vector unsigned int){1u, 2u, 3u, 4u};
    rng.s1 = (vector unsigned int){5u, 6u, 7u, 8u};
    rng.s2 = (vector unsigned int){9u, 10u, 11u, 12u};
    rng.s3 = (vector unsigned int){13u, 14u, 15u, 16u};
    for (step = 0; step < 3; ++step) {
        vector unsigned int actual = spu_rng_next(&rng);
        for (lane = 0; lane < 4; ++lane)
            if (spu_extract(actual, lane) != expected[step][lane]) {
                fprintf(stderr, "xoshiro lane %d step %d failed\n", lane, step);
                return 1;
            }
    }
    spu_rng_seed(&rng, 12345u, 42u);
    {
        vector unsigned int actual = spu_rng_next(&rng);
        for (lane = 0; lane < 4; ++lane)
            if (spu_extract(actual, lane) != seeded[lane]) {
                fprintf(stderr, "SplitMix64 seed lane %d failed\n", lane);
                return 1;
            }
    }
    spu_rng_seed(&rng, 12345u, 43u);
    {
        vector unsigned int actual = spu_rng_next(&rng);
        for (lane = 0; lane < 4; ++lane)
            if (spu_extract(actual, lane) != next_pixel[lane]) {
                fprintf(stderr, "Next pixel seed lane %d failed\n", lane);
                return 1;
            }
    }
    spu_rng_seed(&rng, 12345u, 42u);
    for (step = 0; step < 10000; ++step) {
        vector float actual = spu_rng_uniform(&rng);
        for (lane = 0; lane < 4; ++lane) {
            float value = spu_extract(actual, lane);
            if (!(value >= 0.0f && value < 1.0f)) {
                fprintf(stderr, "Uniform range lane %d failed\n", lane);
                return 1;
            }
        }
    }
    puts("SPU xoshiro128++ test passed");
    return 0;
}
