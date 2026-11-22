# Resolve AMBERHOME and the support libs wherever the archive was extracted,
# so a different unpack location (~/amber vs ~) can't break every script.
# Sourced by the MD scripts. Exits non-zero if Amber cannot be found.
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
for u in /workspace/amber/usr /workspace/usr "$(dirname "$AMBERHOME")/usr"; do
  if [ -d "$u/lib" ]; then
    export LD_LIBRARY_PATH="$u/lib:$u/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"; break
  fi
done
echo "[amber_env] AMBERHOME=$AMBERHOME"

# Resolve the MD engine by ABSOLUTE PATH, not by name. `pmemd.cuda` is listed in
# the archive's bin yet a PATH lookup for it failed with "command not found"
# (job 194210), so the generic name is not dependable here -- a dangling symlink
# would look exactly like this. -x follows symlinks, so a broken one falls
# through to the real SPFP/DPFP binary.
for p in "$AMBERHOME/bin/pmemd.cuda" "$AMBERHOME/bin/pmemd.cuda_SPFP" \
         "$AMBERHOME/bin/pmemd.cuda_DPFP"; do
  if [ -x "$p" ]; then export PMEMD="$p"; break; fi
done
if [ -z "${PMEMD:-}" ]; then
  echo "ERROR: no executable pmemd.cuda/_SPFP/_DPFP under $AMBERHOME/bin"; exit 1
fi
echo "[amber_env] PMEMD=$PMEMD"
