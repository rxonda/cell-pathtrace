#!/usr/bin/env bash
set -euo pipefail

samples=${1:-4}
spus=${2:-1}
seed=${3:-1}
if [[ ! "$samples" =~ ^[0-9]+$ ]] || (( samples < 4 || samples % 4 != 0 )); then
    echo 'samples must be a positive multiple of four' >&2
    exit 2
fi
if [[ ! "$spus" =~ ^[1-6]$ ]]; then
    echo 'SPE workers must be between one and six' >&2
    exit 2
fi
if [[ ! "$seed" =~ ^[0-9]+$ ]]; then
    echo 'seed must be an unsigned decimal integer' >&2
    exit 2
fi

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=${CELL_BUILD_DIR:-"$repo_dir/build"}
mkdir -p "$build_dir"
build_dir=$(cd "$build_dir" && pwd)
for binary in forwardFloatSPU smallptPPU; do
    if [[ ! -s "$build_dir/$binary" ]]; then
        echo "Missing $build_dir/$binary; run bash docker/build-project.sh first." >&2
        exit 2
    fi
done

image=${CELL_DOCKER_IMAGE:-cell-pathtrace:dev-3.1}
if [[ -e "$build_dir/image.ppm" ]]; then
    echo "Refusing to overwrite $build_dir/image.ppm" >&2
    exit 2
fi

docker run --rm \
    --network none \
    --user "$(id -u):$(id -g)" \
    -v "$repo_dir:/project:ro" \
    -v "$build_dir:/build" \
    -e "CELL_SAMPLES=$samples" \
    -e "CELL_SPUS=$spus" \
    -e "CELL_SEED=$seed" \
    "$image" \
    systemsim -cell -n -f /project/docker/render.tcl

if [[ ! -s "$build_dir/image.ppm" ]]; then
    echo 'The guest did not produce image.ppm; inspect the simulator output.' >&2
    exit 1
fi
echo "Rendered image: $build_dir/image.ppm"
