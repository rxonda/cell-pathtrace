# Cell Linux emulation for this project

This program uses `libspe2` on a Linux PPU to launch an SPU image. For an
end-to-end run on x86_64, use IBM's Cell **Full System Simulator** in Linux
mode, together with the matching Cell SDK and guest disk image. The simulator
boots a Cell Linux guest with the SPU runtime. A PS3 homebrew toolchain alone
only cross-compiles code; RPCS3 targets the PS3 environment rather than this
Linux `libspe2` program. QEMU's PowerPC machine list does not include a Cell
machine with SPEs.

IBM's SDK 3.1 documentation lists Fedora 9 x86_64 as a development and
Full-System-Simulator execution platform. The combined Docker image below was
built and tested on an x86_64 Linux host. It contains the historical simulator
and the Cell SDK 3.1 PPU and SPU cross compilers. The simulated guest is Fedora
9 for PowerPC.

## Combined Docker image (tested)

The [historical simulator archive](https://archive.org/details/ibm-full-system-simulator-for-the-cbe-processor)
contains both `systemsim-cell-3.1-25.f9.x86_64.rpm` and the matching
`sysroot_image-3.1-1.noarch.rpm`. Download the ZIP, then build the local image:

```sh
curl -L --fail -o /tmp/cell-systemsim.zip \
  'https://archive.org/download/ibm-full-system-simulator-for-the-cbe-processor/IBM%20Full-System%20Simulator%20for%20the%20CBE%20Processor.zip'
bash docker/build-sdk.sh /tmp/cell-systemsim.zip
```

The script builds a simulator base image, then adds the SDK to produce
`cell-pathtrace:dev-3.1`. Both tools are in that final image. If the simulator
base image is already built, `bash docker/build-sdk.sh` reuses it.
The SDK build uses the simulator's Fedora 9 guest disk for the matching
`libspe2` runtime. It downloads Cell SDK 3.1 compiler,
SIMD math, and library RPMs from the [archived Cell repositories](https://ps3linux.net/cellsdk-repos/),
Fedora 9 PowerPC headers and libraries from the [Fedora archive](https://archives.fedoraproject.org/pub/archive/fedora/linux/releases/9/Everything/ppc/os/Packages/),
and the matching [libspe2 source headers](https://sourceforge.net/projects/libspe/files/libspe/libspe2/).
The scripts check pinned SHA-256 hashes before building. Downloaded archives
stay in `/tmp/cell-docker-rpms` and `/tmp/cell-sdk-rpms` unless the corresponding
cache variables specify other directories.

```sh
docker run --rm -v "$PWD:/work:ro" cell-pathtrace:dev-3.1 \
  bash /work/docker/sdk-smoke.sh
docker run --rm -it -w /work -v "$PWD:/work" cell-pathtrace:dev-3.1
```

The smoke test compiled and linked a SIMD SPU program and a 64-bit PPU
`libspe2` program. `file` identified them as Cell SPU and PowerPC executables.
This image contains the GNU cross compilers, target sysroot, SIMD math headers,
project-relevant Cell libraries, simulator, and guest disk. It does not contain
the entire IBM Developer ISO or IDE. Build binaries in the container shell and
run them inside its simulated Cell Linux guest.

Both project sources compile in this image. Run these commands in the
container from `/work`:

```sh
spu-gcc -O3 forwardFloat-spu.c -o forwardFloatSPU -lm
ppu-gcc -Wall -m64 -maltivec -pthread smallpt-ppu.c -o smallptPPU -lspe2
```

Boot the guest and check the SPE filesystem:

```sh
docker run --rm --network none -v "$PWD:/project:ro" cell-pathtrace:dev-3.1 \
  systemsim -cell -n -f /project/docker/smoke.tcl
```

The test boots Fedora 9, prints `ppc64`, checks that `/spu` is mounted as
`spufs`, then exits. It completed on this x86_64 host. Run a headless
interactive simulator console with:

```sh
docker run --rm -it --network none -v "$PWD:/project:ro" cell-pathtrace:dev-3.1 \
  systemsim -cell -n
```

At the `systemsim %` prompt, enter `mysim go` to boot the guest. For scripted
guest commands, copy the pattern in `docker/smoke.tcl`: source the simulator's
`systemsim.tcl`, call `fast_wait_for_linux_prompt`, queue commands with
`mysim console create UART0 in string`, then call `mysim go`. The simulator
needs a writable working directory for its copy-on-write guest disk files;
the image uses `/sim`, while `/project` is mounted read-only.

To build and run this path tracer, use the scripts from the repository root:

```sh
CELL_WIDTH=16 CELL_HEIGHT=12 bash docker/build-project.sh
bash docker/run-render.sh 4 1 12345
bash docker/test-rng.sh
```

The build script creates the SPU and PPU executables in `build/`. Without
`CELL_WIDTH` and `CELL_HEIGHT`, it builds the default 640 x 480 image. The
render script transfers both executables into the guest with IBM's `callthru`
utility, runs four samples per pixel on one SPE with seed `12345`, and copies `image.ppm` to
`build/image.ppm` on the host. It refuses to overwrite an existing result;
move it before another render. Set `CELL_BUILD_DIR` to choose another host
directory for both binaries and output. The `build/` directory is ignored by
Git. A 64 x 48 render at 5,000 samples per pixel completed in about 14 hours
55 minutes on the tested x86_64 host with two SPE workers and seed `1`. To
reproduce it without overwriting the small test output:

```sh
CELL_BUILD_DIR=build/render-5000-64x48 CELL_WIDTH=64 CELL_HEIGHT=48 \
  bash docker/build-project.sh
CELL_BUILD_DIR=build/render-5000-64x48 bash docker/run-render.sh 5000 2 1
```

The output is `build/render-5000-64x48/image.ppm`. The RNG and review changes
are documented in [CODE_REVIEW.md](CODE_REVIEW.md); image accuracy still needs
comparison with a scalar reference.
Direct light sampling makes low-sample previews visibly closer to the
checked-in example. A 64 x 48 preview with 16 samples per pixel takes much
longer than the small functional test. Its remaining brightness difference
from the example is discussed in [CODE_REVIEW.md](CODE_REVIEW.md).
The RNG test script builds a separate SPU test binary in `build/` and runs it
in the guest without replacing the renderer binaries.
The Cell guest has no compiler; the build script uses the compilers in the
container outside the guest.

## Fedora VM alternative (not tested here)

The following VM instructions follow IBM's SDK 3.1 guide. Keep these old
packages isolated rather than installing them on a modern Fedora host.

### Obtain matching components

Use the Fedora 9 / SDK 3.1 versions of all of these:

1. A Fedora 9 x86_64 installation image for a build VM.
2. `cell-install-3.1.0-0.0.noarch.rpm` and
   `CellSDK-Devel-Fedora_3.1.0.0.0.iso` (the Extras ISO is optional).
3. Optional for simulator use inside the VM: `systemsim-cell-3.1-25.f9.x86_64.rpm`
   and `sysroot_image-3.1-1.noarch.rpm`. Both are in the simulator ZIP above.

The [SDK 3.1 installation guide](https://www.ps3linux.net/ps3-filez/cellsdk-docs/3.1/install/SDK-Installation_Guide_v3.1.pdf)
names the packages and supported platforms. IBM distributed the simulator
separately from SDK 3.1. A [historical archive of the simulator](https://archive.org/details/ibm-full-system-simulator-for-the-cbe-processor)
and a [community index of SDK files](https://www.psdevwiki.com/ps3/Cell_Programming_IBM)
exist, but their package contents, integrity, licensing, and compatibility
must be checked before installation. Do not mix SDK 3.0 and 3.1 components.
An [archived copy of the 3.1 installer RPM](https://ps3linux.net/cellsdk-repos/cellsdk-devel-rhel_3.1.0.0.0/cell-install-3.1.0-0.0.noarch.rpm)
is reachable; its payload includes configuration for both Fedora and RHEL.
I could not verify a working download for the Fedora 9 Developer ISO. The
simulator ZIP does not contain it. Obtain that ISO before starting this VM
installation. The combined Docker image above assembles the required cross sysroot
from Fedora 9 PowerPC RPMs.

### Install and boot

1. Create a Fedora 9 x86_64 VM using your host's VM manager. Give it enough
   disk space for the SDK, simulator, guest image, and builds. Keep a VM
   snapshot after a successful install. A hardware-virtualized VM is
   preferable; software-emulating x86 as well as Cell will be slow.
2. Put the installer RPM and Developer ISO together in `/tmp/cellsdkiso`
   inside the VM. Install IBM's prerequisites, then run the SDK installer:

   ```sh
   sudo yum install rsync sed tcl wget
   sudo rpm -ivh /tmp/cellsdkiso/cell-install-3.1.0-0.0.noarch.rpm
   sudo /opt/cell/cellsdk --iso /tmp/cellsdkiso install
   export PATH=/opt/cell/toolchain/bin:$PATH
   command -v spu-gcc ppu-gcc
   ```

   The original Fedora 9 YUM mirror URLs may no longer work; point them at
   the [Fedora archive](https://archives.fedoraproject.org/pub/archive/fedora/linux/releases/9/Fedora/x86_64/)
   if package installation fails. IBM's installer also needs its matching
   open-source package repository. If that repository is unavailable, the
   installation is incomplete even when the Developer ISO is present.

   Read and accept the licenses presented by the installer. The simulator is
   already installed in the Docker image, so the VM only needs the SDK to
   cross-compile. If you also want to run the simulator inside the VM, install
   the matching `systemsim-cell` and `sysroot_image` RPMs. Confirm that
   `/opt/ibm/systemsim-cell/bin/systemsim` and
   `/opt/ibm/systemsim-cell/lib/cell/systemsim.tcl` exist before booting it.
3. If you installed the simulator in the VM, start it from a private directory:

   ```sh
   mkdir -p ~/cell-sim
   cd ~/cell-sim
   export PATH=/opt/ibm/systemsim-cell/bin:$PATH
   systemsim -cell -g
   ```

   In the simulator GUI select **Fast Mode**, then **Go** to boot Cell Linux.
   The console window should reach a Linux shell. The command-line interface
   also accepts `mysim go` at the `systemsim %` prompt. These steps come from
   IBM's [Programming Tutorial, simulator chapter](https://www.cresco.enea.it/LA1/cresco_sp14_ylichron/CBE-docs/CBE_Programming_Tutorial_v3.1.pdf).

### Build and run a smoke test

Use one of the SDK's supplied PPU/SPU `libspe2` examples first. It verifies
the guest, compiler, SPU loader, and runtime separately from this project's
known code defects. Build with the SDK's `spu-gcc` and PPU cross compiler in
the Fedora VM. On x86_64, the PPU compiler is `ppu-gcc`; host `gcc` produces
an x86 binary that cannot run in the Cell guest.

The corresponding cross-build commands are:

```sh
spu-gcc -O3 forwardFloat-spu.c -o forwardFloatSPU -lm
ppu-gcc -Wall -m64 -maltivec -pthread smallpt-ppu.c -o smallptPPU -lspe2
file forwardFloatSPU smallptPPU
```

Both commands were verified with the combined Docker image, and a small render
completed. The review and remaining validation are in [CODE_REVIEW.md](CODE_REVIEW.md).
Copy both binaries from the VM into this repository's `build/` directory on the
host, then run `bash docker/run-render.sh 4 1`.

To transfer a built PPU executable from the build VM's `/tmp` into the
simulated Linux console, IBM documents:

```sh
callthru source /tmp/example-ppu > example-ppu
chmod +x example-ppu
./example-ppu
```

Use the same transfer for this project's `smallptPPU` and
`forwardFloatSPU`, keeping both in the same guest directory because the PPU
program opens `./forwardFloatSPU`. Verify the guest has `libspe2` and
`spufs` before a render (for example, `mount | grep spufs` and
`ldd ./smallptPPU`). Start with a tiny image and four samples per pixel;
`CELL_WIDTH` and `CELL_HEIGHT` configure the build script for a practical
first simulation run.

## Scope of the simulator

Use functional or Fast Mode for correctness debugging. IBM's pipeline mode
can collect SPE instruction statistics, but simulated elapsed time is not a
substitute for a render benchmark on real Cell hardware. Retain a host-native
reference test for RNG and sampling math so most changes can be checked
without a full guest boot.
