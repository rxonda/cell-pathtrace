#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cache_dir=${CELL_SDK_RPM_CACHE:-/tmp/cell-sdk-rpms}
simulator_archive=${1:-}

if ! docker info >/dev/null 2>&1; then
    echo 'Docker daemon is unavailable; start Docker and check your account can access it.' >&2
    exit 2
fi

if [[ -n "$simulator_archive" ]]; then
    bash "$repo_dir/docker/build.sh" "$simulator_archive"
fi

if ! docker image inspect cell-pathtrace:systemsim-3.1 >/dev/null 2>&1; then
    echo 'Provide the simulator archive: bash docker/build-sdk.sh /path/to/IBM-Full-System-Simulator.zip' >&2
    exit 2
fi

bash "$repo_dir/docker/download-sdk-rpms.sh"
(cd "$cache_dir" && sha256sum -c "$repo_dir/docker/sdk-rpms.sha256")

docker build \
    --build-context "cell_sdk_rpms=$cache_dir" \
    -f "$repo_dir/docker/Dockerfile.sdk" \
    -t cell-pathtrace:dev-3.1 \
    -t cell-pathtrace:sdk-3.1 \
    "$repo_dir/docker"
