#ifndef SPU_RNG_H
#define SPU_RNG_H

#include <stdint.h>
#include <spu_intrinsics.h>

typedef struct {
    vector unsigned int s0, s1, s2, s3;
} SpuRng;

/* SplitMix64 is used only to initialize the four independent lane states. */
static inline uint64_t splitmix64_next(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

static inline void spu_rng_seed(SpuRng *rng, uint64_t render_seed,
                                uint64_t pixel_id) {
    unsigned int words[4][4];
    int lane;
    for (lane = 0; lane < 4; ++lane) {
        uint64_t stream = pixel_id * 4 + lane;
        uint64_t seed = render_seed ^
            (stream * UINT64_C(0x9e3779b97f4a7c15));
        uint64_t a = splitmix64_next(&seed);
        uint64_t b = splitmix64_next(&seed);
        words[0][lane] = (unsigned int)a;
        words[1][lane] = (unsigned int)(a >> 32);
        words[2][lane] = (unsigned int)b;
        words[3][lane] = (unsigned int)(b >> 32);
        if (!(words[0][lane] | words[1][lane] |
              words[2][lane] | words[3][lane]))
            words[0][lane] = 1;
    }
    rng->s0 = (vector unsigned int){words[0][0], words[0][1], words[0][2], words[0][3]};
    rng->s1 = (vector unsigned int){words[1][0], words[1][1], words[1][2], words[1][3]};
    rng->s2 = (vector unsigned int){words[2][0], words[2][1], words[2][2], words[2][3]};
    rng->s3 = (vector unsigned int){words[3][0], words[3][1], words[3][2], words[3][3]};
}

static inline vector unsigned int spu_rng_next(SpuRng *rng) {
    /* xoshiro128++ 1.0, one 32-bit result per SIMD lane. */
    vector unsigned int sum = rng->s0 + rng->s3;
    vector unsigned int result =
        spu_or(spu_sl(sum, spu_splats(7u)),
               spu_rlmask(sum, spu_splats(-25))) + rng->s0;
    vector unsigned int t = spu_sl(rng->s1, spu_splats(9u));
    rng->s2 = spu_xor(rng->s2, rng->s0);
    rng->s3 = spu_xor(rng->s3, rng->s1);
    rng->s1 = spu_xor(rng->s1, rng->s2);
    rng->s0 = spu_xor(rng->s0, rng->s3);
    rng->s2 = spu_xor(rng->s2, t);
    rng->s3 = spu_or(spu_sl(rng->s3, spu_splats(11u)),
                     spu_rlmask(rng->s3, spu_splats(-21)));
    return result;
}

static inline vector float spu_rng_uniform(SpuRng *rng) {
    /* Mantissa construction yields [0, 1) without integer-to-float conversion. */
    vector unsigned int exponent = {0x3f800000u, 0x3f800000u,
                                     0x3f800000u, 0x3f800000u};
    vector unsigned int bits = spu_rng_next(rng);
    bits = spu_rlmask(bits, spu_splats(-9));
    bits = spu_or(bits, exponent);
    return (vector float)bits - spu_splats(1.0f);
}

#endif
