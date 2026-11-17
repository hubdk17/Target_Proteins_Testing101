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
