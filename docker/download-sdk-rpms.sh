#!/usr/bin/env bash
set -euo pipefail

cache_dir=${CELL_SDK_RPM_CACHE:-/tmp/cell-sdk-rpms}
mkdir -p "$cache_dir"

open_base=https://ps3linux.net/cellsdk-repos/cellsdk-open-rhel_3.1.0.0.0/x86_64
devel_base=https://ps3linux.net/cellsdk-repos/cellsdk-devel-rhel_3.1.0.0.0/x86_64
fedora_base=https://archives.fedoraproject.org/pub/archive/fedora/linux/releases/9/Everything/ppc/os/Packages

download() {
    local url=$1
    local name=${url##*/}
    if [[ ! -s "$cache_dir/$name" ]]; then
        curl -L --fail --silent --show-error --retry 2 --output "$cache_dir/$name.part" "$url"
        mv "$cache_dir/$name.part" "$cache_dir/$name"
    fi
}

for name in \
    ppu-binutils-2.18.50-21.i686.rpm \
    ppu-gcc-4.1.1-166.i686.rpm \
    spu-binutils-2.18.50-21.i686.rpm \
    spu-gcc-4.1.1-166.i686.rpm \
    spu-newlib-1.16.0-17.i686.rpm; do
    download "$open_base/$name"
done

for name in \
    cell-libs-cross-3.1-1.noarch.rpm \
    cell-libs-cross-devel-3.1-1.noarch.rpm \
    simdmath-cross-devel-3.1-1.noarch.rpm; do
    download "$devel_base/$name"
done

for name in \
    glibc-2.8-3.ppc64.rpm \
    glibc-devel-2.8-3.ppc64.rpm \
    glibc-headers-2.8-3.ppc.rpm \
    kernel-headers-2.6.25-14.fc9.ppc.rpm; do
    download "$fedora_base/$name"
done

download https://downloads.sourceforge.net/project/libspe/libspe/libspe2/libspe2-2.2.80-132.tar.gz

printf 'Downloaded SDK build inputs to %s\n' "$cache_dir"
