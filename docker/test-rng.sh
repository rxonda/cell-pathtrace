#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p "$repo_dir/build"
test_dir=$(mktemp -d "$repo_dir/build/rng-test.XXXXXX")
image=${CELL_DOCKER_IMAGE:-cell-pathtrace:dev-3.1}

docker run --rm --user "$(id -u):$(id -g)" \
    -v "$repo_dir:/project:ro" -v "$test_dir:/build" "$image" \
    bash -c 'spu-gcc -O3 -I/project /project/tests/spu-rng-test.c -o /build/rng-test-spu && ppu-gcc -Wall -m64 /project/tests/rng-ppu-test.c -o /build/rng-test-ppu -lspe2'

docker run --rm --network none --user "$(id -u):$(id -g)" \
    -v "$repo_dir:/project:ro" -v "$test_dir:/build" "$image" \
    systemsim -cell -n -f /project/docker/rng-test.tcl > "$test_dir/test.log" 2>&1

if ! rg -q 'SPU xoshiro128\+\+ test passed' "$test_dir/test.log" ||
   ! rg -q 'RNG launcher passed' "$test_dir/test.log"; then
    echo "RNG test failed; see $test_dir/test.log" >&2
    exit 1
fi
echo "SPU RNG test passed; log: $test_dir/test.log"
