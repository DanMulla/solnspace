#!/bin/bash
# Build the PolytopeSamplerMatlab PackedChol MEX files for macOS on Apple silicon.
#
# Why this exists
# ---------------
# PolytopeSamplerMatlab compiles PackedChol on first use via its own compile_solver.m,
# which calls MATLAB's `mex`. On macOS, MATLAB's C++ mex configuration
# (bin/maca64/mexopts/clang++_maca64.xml) resolves XCODE_AGREED_VERSION by reading
# Xcode's own preferences, so it requires a full Xcode install with an accepted
# licence. With only the Xcode Command Line Tools present, MATLAB reports
#
#     Error using compile_solver (line 29)
#     No C++ mex compiler is available.
#
# even though a perfectly good clang++ is installed. This script compiles the MEX
# files directly with that clang++, using the same flags MATLAB's mexopts XML
# specifies, and writes them into PolytopeSamplerMatlab/bin/. Once they are there,
# compile_solver finds them (exist(name)==3), probes them, and never invokes `mex`.
#
# The alternative fix is to install full Xcode from the App Store, run
#   sudo xcodebuild -license accept
# then `mex -setup C++` in MATLAB. That is the supported route; this script just
# avoids the ~10 GB download.
#
# Usage:  bash build_polytopesampler_macos.sh [/path/to/MATLAB.app]

set -euo pipefail

MR="${1:-/Applications/MATLAB_R2024b.app}"
ROOT="$(cd "$(dirname "$0")" && pwd)"
S="$ROOT/PolytopeSamplerMatlab/code/solver"
BIN="$ROOT/PolytopeSamplerMatlab/bin"
OBJ="$(mktemp -d)"
trap 'rm -rf "$OBJ"' EXIT

[ -d "$MR" ]  || { echo "MATLAB not found at $MR"; exit 1; }
[ -d "$S" ]   || { echo "solver sources not found at $S"; exit 1; }
command -v clang++ >/dev/null || { echo "clang++ not found; run: xcode-select --install"; exit 1; }

SDK="$(xcrun -sdk macosx --show-sdk-path)"
ARCH="$(uname -m)"
if [ "$ARCH" != "arm64" ]; then
  echo "This script targets Apple silicon (arm64); detected $ARCH." >&2
  exit 1
fi

# Flags mirror $MATLABROOT/bin/maca64/mexopts/clang++_maca64.xml
CXXF=(-fno-common -arch arm64 -mmacosx-version-min=12.0 -fexceptions -isysroot "$SDK"
      -fwrapv -ffp-contract=off -std=c++14 -stdlib=libc++ -O2 -DNDEBUG -DMATLAB_MEX_FILE)
LDF=(-Wl,-twolevel_namespace -arch arm64 -mmacosx-version-min=12.0 -Wl,-syslibroot,"$SDK"
     -framework Cocoa -bundle -stdlib=libc++ -O)
EXP=(-Wl,-exported_symbols_list,"$MR/extern/lib/maca64/mexFunction.map"
     -Wl,-exported_symbols_list,"$MR/extern/lib/maca64/c_exportsmexfileversion.map")
LIBS=(-L"$MR/bin/maca64" -weak-lmx -weak-lmex -weak-lmat
      -L"$MR/extern/bin/maca64" -weak-lMatlabDataArray)
QD=(util bits dd_real dd_const qd_real qd_const)

echo "MATLAB : $MR"
echo "SDK    : $SDK"
echo "output : $BIN"

# qd library and MATLAB's mex api version stub (mex adds the latter automatically)
for f in "${QD[@]}"; do
  clang++ -c "${CXXF[@]}" -I"$MR/extern/include" -I"$S" "$S/qd/$f.cc" -o "$OBJ/$f.o"
done
clang++ -c "${CXXF[@]}" -I"$MR/extern/include" \
        "$MR/extern/version/cpp_mexapi_version.cpp" -o "$OBJ/mexapi_version.o"

OBJS=("$OBJ/mexapi_version.o")
for f in "${QD[@]}"; do OBJS+=("$OBJ/$f.o"); done

# SIMD_LEN 0 and 4 are what sample.m requests (0 and opts.simdLen, default 4);
# 1 is built too so opts.simdLen = 1 also works. Name must match
# MexSolver.solverName(), which appends 'arm' on Apple silicon.
for L in 0 1 4; do
  clang++ -c "${CXXF[@]}" -DMATLAB_DEFAULT_RELEASE=R2018a -DSIMD_LEN="$L" \
    -I"$MR/extern/include" -I"$S" "$S/PackedChol.cpp" -o "$OBJ/PackedChol$L.o"
  clang++ "${LDF[@]}" "${EXP[@]}" "$OBJ/PackedChol$L.o" "${OBJS[@]}" "${LIBS[@]}" \
    -o "$BIN/PackedChol${L}arm.mexmaca64"
  echo "  built PackedChol${L}arm.mexmaca64"
done

echo
echo "Done. Verify in MATLAB with:"
echo "  addpath(genpath('$ROOT')); compile_solver(0); compile_solver(4)"
