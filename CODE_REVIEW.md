# Code review and implementation status

This review began with the historical PPU/SPU source. The changes below now
cross-compile with the Cell SDK 3.1 Docker image and have been exercised in
IBM's Cell full-system simulator. Simulator elapsed time is not a Cell hardware
performance measurement.

## Correctness changes applied

| Original finding | Change |
| --- | --- |
| Uninitialized SPU row pointer table and unchecked allocations | Allocate two row buffers from an initialized pointer array, check both allocations, and free them after DMA completes. |
| Invalid sphere ID in an all-miss SIMD lane | Initialize intersection IDs and assign a safe sphere ID to every miss lane before gathering scalar sphere data. |
| SIMD sphere table disagreed with the scalar scene | Rebuild the SIMD sphere table from the scalar spheres at SPU startup. |
| Camera tent filter lacked the lower branch's `-1` | Restore the signed tent distribution in both camera coordinates. |
| Glass Schlick term lost `R0` | Use `R0 + (1-R0)c^5` for the reflection probability. |
| Worker thread signatures, overwritten handles, unsynchronized completion polling | Use a correct pthread function signature, join each SPE worker, and check SPE exit status before writing output. |
| Dead mailbox stack controller with unchecked pop and allocation | Remove it; the SPU code never sends the corresponding mailbox messages. |
| Unchecked arguments and PPM I/O | Validate samples, SPE count, and seed; check `fopen`, writes, and `fclose`. |
| Unnecessary random draws and material math | Draw roulette values only after depth five. Execute glass or diffuse calculations only when an active SIMD lane needs that material. Mixed-material packs still compute the selected material for all four lanes and mask the result. |
| Very sparse low-sample images from the small light source | Sample the light sphere by solid angle at each diffuse hit, test shadow visibility, and add the Lambertian contribution with its sampling weight. Suppress the corresponding next-bounce emitter hit to avoid counting it twice. |

The roulette termination path now returns the masked accumulated radiance,
including emission already collected by live lanes. The simulator render is a
functional test; image quality still needs comparison with a scalar reference
at useful sample counts.

On a 16 x 12, four-sample preview, direct light sampling raised the nonblack
pixel count from 18 to 191 of 192. Mean absolute RGB channel difference from
the checked-in 5,000-sample image, downsampled to 16 x 12, fell from about 64
to 17. A 32 x 24, eight-sample preview had an error of about 15 by the same
measure. These comparisons show a closer preview, not estimator convergence.
The 64 x 48, eight-sample preview had 3,067 nonblack pixels of 3,072 and an
average channel difference of about 15 from the target downsampled to that
size. The checked-in target is 640 x 480 at 5,000 samples per pixel, so the
preview still has much less spatial and sampling detail.
At 16 samples per pixel at the same 64 x 48 size, the error was about 14.7
instead of 15.1. More samples reduced some noise, but this small change also
indicates a systematic difference from the checked-in target. The target may
have been produced by a different source revision; visual comparison alone
cannot establish which estimator is correct.

A 64 x 48 render at 5,000 samples per pixel, with two SPE workers and seed
`1`, completed in the Cell full-system simulator. The container exited with
status zero after about 14 hours 55 minutes on the tested x86_64 host. The
output PPM loads as a complete 64 x 48 RGB image. Against the checked-in
640 x 480 image downsampled with a box filter, mean absolute RGB channel
difference is 14.45. The new image's mean RGB is (88.3, 75.3, 88.3), versus
(73.1, 62.6, 74.2) for the downsampled example. This persistent brightness
difference does not appear to be sampling noise.

At native size, the checked-in example visibly has grain and repeating
ripples on the walls and floor; the new 64 x 48 render looks smoother. Direct
light sampling should reduce variance at diffuse surfaces. The changed RNG
may affect structured artifacts, but these runs do not isolate its effect.
Different output resolutions also prevent a controlled noise comparison.

## SIMD random generator

`spu-rng.h` implements the authors' [xoshiro128++ 1.0 transition](https://prng.di.unimi.it/xoshiro128plusplus.c)
with four `vector unsigned int` state registers. Each SIMD lane has its own
128-bit state. Two SplitMix64 outputs initialize the four words for each
`(pixel, lane)` stream from a 64-bit render seed. State continues across sample
packs for that pixel. Pixel seeds do not depend on which SPE claims a row.

The generator turns the upper 23 output bits into a floating-point mantissa,
then subtracts one to produce values in `[0,1)`. It avoids the old serial
`erand48()` calls and double-to-float insertion. Two streams could theoretically
overlap; independently mixed seeds do not provide the mathematical separation
guarantee of xoshiro's jump operation or a counter-based generator.

The PPU accepts an optional third argument for the decimal or `0x` render
seed. `docker/run-render.sh` accepts it as a third decimal argument; both
default to `1`. For example:

```sh
CELL_WIDTH=16 CELL_HEIGHT=12 bash docker/build-project.sh
bash docker/run-render.sh 4 1 12345
```

The dedicated `tests/spu-rng-test.c` compares three xoshiro steps in all four
SIMD lanes against scalar reference values, checks SplitMix64-based stream
seeding for two pixels, and verifies 40,000 generated uniforms stay in `[0,1)`.
Run `bash docker/test-rng.sh` to build and run it in the Cell guest. The test
passed there. Two 16 x 12, four-sample renders with the same seed produced
byte-identical PPM files using one and two SPE workers.

## Remaining validation

- Compare the rendered scene against a scalar reference on identical camera
  samples. This will test intersections, throughput, and material branches
  beyond the focused fixes above, and help explain the remaining brightness
  difference from the checked-in target.
- Compare noise at matched resolution and scene settings. Stratified camera
  sampling may reduce variance further, but it changes the estimator and
  requires its own correctness checks.
- Measure RNG and complete-render time on Cell hardware before claiming a
  speedup. The simulator confirms execution, not real SPE throughput.
