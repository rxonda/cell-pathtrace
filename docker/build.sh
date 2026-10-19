#!/usr/bin/env bash
set -euo pipefail

archive=${1:-}
if [[ ! -f "$archive" ]]; then
    echo "Usage: docker/build.sh /path/to/IBM-Full-System-Simulator.zip" >&2
    exit 2
fi

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
rpm_dir=${CELL_DOCKER_CACHE:-/tmp/cell-docker-rpms}
mkdir -p "$rpm_dir"
unzip -j -n "$archive" \
    '*/systemsim-cell-3.1-25.f9.x86_64.rpm' \
    '*/sysroot_image-3.1-1.noarch.rpm' -d "$rpm_dir"

printf '%s\n' \
    'd488ce454138c8740a817ef117dc65ad5a1f68a373a141da3a872f69f9426bbe  systemsim-cell-3.1-25.f9.x86_64.rpm' \
    '65e45773ff0d42bc1b64a732d44a0ef3af8e0ee35553947e2ee73f30d85ca872  sysroot_image-3.1-1.noarch.rpm' \
    | (cd "$rpm_dir" && sha256sum -c -)

docker build \
    --build-context "cell_rpms=$rpm_dir" \
    -t cell-pathtrace:systemsim-3.1 \
    "$repo_dir/docker"
