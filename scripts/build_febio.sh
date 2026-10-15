#!/usr/bin/env bash
# Reproducible FEBio 4.13 build (official source, tag v4.13) with Intel MKL PARDISO.
# Verified on Ubuntu 24.04 / GCC 13 / CMake 3.28.  Installs to $PREFIX (default /opt/febio-4.13)
# and creates a wrapper $PREFIX/bin/febio4 that sets LD_LIBRARY_PATH.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PREFIX=${PREFIX:-/opt/febio-4.13}
WORK=${WORK:-$(mktemp -d)}
TAG=v4.13
EXPECTED_COMMIT=32ae206ff4881dfb54f62296cd1558e58ed9fcc6
MKL=${MKL:-/opt/mkl}
# MKL 2024.2 still ships the legacy sparse BLAS (mkl_dcsrtrsv) that FEBio 4.13 uses; MKL 2025+ does not.
python3 -m pip install --break-system-packages --ignore-installed --prefix "$MKL" \
  "mkl==2024.2.2" "mkl-devel==2024.2.2" "mkl-include==2024.2.2" "intel-openmp==2024.2.1"
MKLLIB="$MKL/local/lib"; [ -d "$MKLLIB" ] || MKLLIB="$MKL/lib"
MKLINC="$(dirname "$MKLLIB")/include"
for n in mkl_intel_lp64 mkl_intel_thread mkl_core; do ln -sf "lib$n.so.2" "$MKLLIB/lib$n.so"; done
git clone --depth 1 --branch "$TAG" https://github.com/febiosoftware/FEBio.git "$WORK/febio"
cd "$WORK/febio"
test "$(git rev-parse HEAD)" = "$EXPECTED_COMMIT" || { echo "unexpected commit"; exit 1; }
# Build patch (documented): allow FindLib to accept shared MKL libraries on Linux.
sed -i 's/find_library(TEMP NAMES lib${libName}.a lib${ARGV3}.a/find_library(TEMP NAMES lib${libName}.a lib${libName}.so lib${ARGV3}.a/' CMakeLists.txt
# Source patch (documented in docs/solver_selection.md): two restart (-r) bugs in FEBio 4.13 that
# crash every resumed model containing rigid bodies or logfile face_data:
#   * FERigidSolidDomain::Serialize skipped the domain material/elements -> null material in plot output
#   * FaceDataRecord did not serialize its surface -> null surface in logfile output
# The patch changes only serialization; it does not touch any constitutive, contact or solver code.
git apply "$SCRIPT_DIR/patches/febio-4.13-restart-fixes.patch"
mkdir -p cbuild && cd cbuild
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_MKL=ON -DMKLROOT="$(dirname "$MKLLIB")" -DMKL_INC="$MKLINC" \
  -DMKL_LIB_DIR="$MKLLIB" -DMKL_OMP_LIB="$MKLLIB/libiomp5.so" -DUSE_ZLIB=OFF -DUSE_STATIC_STDLIBS=OFF
make -j"$(nproc)" febio4
mkdir -p "$PREFIX/bin" "$PREFIX/lib"
cp bin/febio4 "$PREFIX/bin/febio4.bin"
cp lib/*.so "$PREFIX/lib/"
cp "$MKLLIB"/libmkl_core.so.2 "$MKLLIB"/libmkl_intel_lp64.so.2 "$MKLLIB"/libmkl_intel_thread.so.2 "$MKLLIB"/libiomp5.so \
   "$MKLLIB"/libmkl_def.so.2 "$MKLLIB"/libmkl_avx2.so.2 "$MKLLIB"/libmkl_avx512.so.2 "$PREFIX/lib/" 2>/dev/null || true
cat > "$PREFIX/bin/febio4" <<'W'
#!/bin/sh
HERE="$(dirname "$(readlink -f "$0")")"
export LD_LIBRARY_PATH="$HERE/../lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec "$HERE/febio4.bin" "$@"
W
chmod +x "$PREFIX/bin/febio4"
echo "FEBio installed: $PREFIX/bin/febio4  (export RLFEA_FEBIO=$PREFIX/bin/febio4 or put it on PATH)"
