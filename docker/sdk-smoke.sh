#!/usr/bin/env bash
set -euo pipefail

cat > /tmp/sdk-spu-smoke.c <<'C'
#include <spu_intrinsics.h>
#include <simdmath.h>
#include <simdmath/sqrtf4.h>
int main(void) {
    vector float x = spu_splats(4.0f);
    volatile vector float y = _sqrtf4(x);
    (void)y;
    return 0;
}
C

cat > /tmp/sdk-ppu-smoke.c <<'C'
#include <stdio.h>
#include <libspe2.h>
int main(void) {
    puts("Cell PPU compiler is ready");
    return spe_context_create(0, NULL) == NULL;
}
C

spu-gcc -O2 /tmp/sdk-spu-smoke.c -o /tmp/sdk-spu-smoke -lm
ppu-gcc -m64 -pthread /tmp/sdk-ppu-smoke.c -o /tmp/sdk-ppu-smoke -lspe2
file /tmp/sdk-spu-smoke /tmp/sdk-ppu-smoke
