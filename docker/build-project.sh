#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=${CELL_BUILD_DIR:-"$repo_dir/build"}
width=${CELL_WIDTH:-640}
height=${CELL_HEIGHT:-480}
image=${CELL_DOCKER_IMAGE:-cell-pathtrace:dev-3.1}

for value in "$width" "$height"; do
    if [[ ! "$value" =~ ^[0-9]+$ ]] || (( value < 2 )); then
        echo 'CELL_WIDTH and CELL_HEIGHT must be integers of at least 2' >&2
        exit 2
    fi
done
mkdir -p "$build_dir"
build_dir=$(cd "$build_dir" && pwd)

docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$repo_dir:/project:ro" \
    -v "$build_dir:/build" \
    -e "CELL_WIDTH=$width" \
    -e "CELL_HEIGHT=$height" \
    "$image" \
    bash -c 'spu-gcc -O3 /project/forwardFloat-spu.c -o /build/forwardFloatSPU -lm && ppu-gcc -Wall -m64 -maltivec -pthread -DWIDTH="$CELL_WIDTH" -DHEIGHT="$CELL_HEIGHT" /project/smallpt-ppu.c -o /build/smallptPPU -lspe2 && file /build/forwardFloatSPU /build/smallptPPU'

printf 'Built %s x %s binaries in %s\n' "$width" "$height" "$build_dir"
