# Resolve AMBERHOME and every library pmemd.cuda needs, wherever the archive
# was extracted. Sourced by the MD scripts; exits non-zero if Amber or the
# CUDA engine cannot be found. Runs INSIDE the nvhpc container.
#
# Measured on this cluster, 2026-11-19..24, one job per lesson:
#  - the container has NO python, so nothing here may call it (job 188654);
#  - `pmemd.cuda` is a DANGLING symlink; the real binary is pmemd.cuda_SPFP,
#    so the engine is resolved by absolute path and tested with -x (194210);
#  - amber_full_libs.tar.gz extracts a COMPLETE system library tree to ~/usr,
#    libc included. Putting /workspace/usr/lib/x86_64-linux-gnu on
#    LD_LIBRARY_PATH gave every binary the wrong glibc and ended in a
#    segmentation fault. It must NEVER go on the library path.
#  - that tree's .so entries are ABSOLUTE symlinks (libopenblas.so.0 ->
#    /etc/alternatives/...), so they dangle inside a container even when the
#    directory is on the path. Hence the shim directory below: our own
#    symlinks, only the sonames Amber needs, no libc.
#  - the nvhpc image keeps libcufft/libcurand/libcusolver under its HPC SDK
#    tree, not on the default path.
for cand in /workspace/amber/amber24 /workspace/amber24 "${HOME:-/root}/amber/amber24"; do
  if [ -f "$cand/amber.sh" ]; then export AMBERHOME="$cand"; break; fi
done
if [ -z "${AMBERHOME:-}" ] || [ ! -f "$AMBERHOME/amber.sh" ]; then
  echo "ERROR: cannot find amber.sh (looked under /workspace/amber/amber24, /workspace/amber24)"; exit 1
fi
# amber.sh references unset vars (PERL5LIB, ...). Our callers run `set -u`,
# which turns those into fatal errors, so disable nounset across the source
# and restore it (all callers here use set -u).
set +u
source "$AMBERHOME/amber.sh"
set -u
echo "[amber_env] AMBERHOME=$AMBERHOME"

# 1. Amber's own libraries: libemil.so, libkmmd.so, libnetcdff.so.6
export LD_LIBRARY_PATH="$AMBERHOME/lib:${LD_LIBRARY_PATH:-}"

# 2. CUDA maths libraries, off the default path in this image
for m in /opt/nvidia/hpc_sdk/Linux_x86_64/*/math_libs/lib64; do
  [ -d "$m" ] && export LD_LIBRARY_PATH="$m:$LD_LIBRARY_PATH"
done

# 3. OpenBLAS, via a shim of our own holding the soname only
SHIM="$(dirname "$AMBERHOME")/libshim"
if [ ! -e "$SHIM/libopenblas.so.0" ]; then
  for r in /workspace/usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r*.so \
           /workspace/usr/lib/x86_64-linux-gnu/openblas-*/libopenblas*.so; do
    if [ -f "$r" ]; then mkdir -p "$SHIM"; ln -sf "$r" "$SHIM/libopenblas.so.0"; break; fi
  done
fi
if [ -e "$SHIM/libopenblas.so.0" ]; then
  export LD_LIBRARY_PATH="$SHIM:$LD_LIBRARY_PATH"
else
  echo "[amber_env] WARNING: no libopenblas.so.0 shim built -- pmemd will likely fail to load"
fi

# The MD engine, by absolute path: `pmemd.cuda` is a dangling symlink here and
# -x follows symlinks, so a broken one falls through to the real binary. CPU
# pmemd is deliberately NOT a fallback: a silent switch would turn a ~12 h GPU
# run into weeks and read as a hung job.
for p in "$AMBERHOME/bin/pmemd.cuda" "$AMBERHOME/bin/pmemd.cuda_SPFP" \
         "$AMBERHOME/bin/pmemd.cuda_DPFP"; do
  if [ -x "$p" ]; then export PMEMD="$p"; break; fi
done
if [ -z "${PMEMD:-}" ]; then
  echo "ERROR: no executable pmemd.cuda/_SPFP/_DPFP under $AMBERHOME/bin"; exit 1
fi
echo "[amber_env] PMEMD=$PMEMD"
