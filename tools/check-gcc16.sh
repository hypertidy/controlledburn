#!/usr/bin/env bash
# Reproduce the CRAN GCC 16 -Warray-bounds WARN and confirm a fix.
#
# CRAN's r-devel-linux-x86_64-{debian,fedora}-gcc flavors use GCC 16. The
# r-hub gcc16 container (g++ 16.2.1 on Fedora, R-devel) is the same toolchain,
# so we can reproduce the compile-time warning locally and A/B a fix against
# a pinned image. Scoped to R CMD INSTALL: the WARN is emitted at the
# "whether package can be installed" step, during compilation, so a full
# --as-cran check (which drags in terra/GDAL etc.) is not needed just to
# confirm the warning is gone.
#
# Usage:
#   tools/check-gcc16.sh [REF]        # default REF: HEAD
#   tools/check-gcc16.sh HEAD HEAD~1  # compare two refs (post then pre)
#
# Exit status is non-zero if any built ref still emits the -Warray-bounds
# warning, so this is usable as a CI gate.

set -euo pipefail

IMAGE="ghcr.io/r-hub/containers/gcc16:latest"
REPO_ROOT="$(git -C "$(dirname "$0")/.." rev-parse --show-toplevel)"
REFS=("${@:-HEAD}")
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

docker pull "$IMAGE"

# Build the runtime deps once into a shared library the installs reuse.
mkdir -p "$WORK/rlib"
docker run --rm --user root -v "$WORK:/work" "$IMAGE" bash -lc '
  Rscript -e "install.packages(c(\"cpp11\",\"vaster\",\"wk\",\"ximage\"), lib=\"/work/rlib\", repos=\"https://cloud.r-project.org\")"
'

status=0
for ref in "${REFS[@]}"; do
  name="$(echo "$ref" | tr -c 'A-Za-z0-9._-' '_')"
  git -C "$REPO_ROOT" archive "$ref" | tar -x -C "$WORK" --one-top-level="src_$name"
  echo "#################### $ref ####################"
  docker run --rm --user root -v "$WORK:/work" -e R_LIBS_USER=/work/rlib "$IMAGE" bash -lc "
    cd /work/src_$name && R CMD INSTALL --preclean --no-docs . > /work/$name.log 2>&1
  " || true
  if grep -q "array-bounds\|outside array bounds" "$WORK/$name.log"; then
    echo "WARN present:"
    grep -n "outside array bounds" "$WORK/$name.log"
    status=1
  else
    echo "clean: no -Warray-bounds warning"
  fi
  echo
done

exit $status
