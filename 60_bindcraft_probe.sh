#!/bin/bash
# BindCraft feasibility probe. READ-ONLY: downloads nothing, installs nothing,
# needs no GPU and no Slurm allocation. Run it on the login node.
#
#   cd ~/rfd && bash 60_bindcraft_probe.sh
#
# Why a probe at all: BindCraft ships `install_bindcraft.sh --pkg_manager conda`
# and this cluster has no conda, so the shipped installer cannot be the route.
# The plan is to layer it on an existing container that already carries JAX and
# PyRosetta. That plan has four load-bearing assumptions, and this script is
# the one round that either confirms or kills all four at once:
#
#   A. the login node can reach the AF2 parameter host (compute nodes are
#      host-filtered, so anything needed must be fetched here)
#   B. we do not already have AF2 *multimer* parameters on disk
#   C. one of our containers carries a JAX that ColabDesign can use
#   D. there is room for 5.3 GB of parameters plus trajectory output
#
# Everything is printed. Nothing is inferred from silence: each check says
# either a value or why it could not get one.

set -uo pipefail   # deliberately not -e: a failing probe is a RESULT, not an
                   # error, and the remaining checks still need to run.

RFD="$HOME/rfd"
say() { printf '\n== %s ==\n' "$1"; }

printf 'host: %s   user: %s   date: %s\n' "$(hostname)" "$USER" "$(date -Is)"

# ---------------------------------------------------------------- A. host reach
say "A. can this node reach what BindCraft needs"
# The AF2 parameter tarball, GitHub (repo + ColabDesign), and PyPI. -I so we
# fetch headers only: a few hundred bytes, not 5.3 GB.
for URL in \
  "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar" \
  "https://github.com/martinpacesa/BindCraft" \
  "https://pypi.org/simple/" ; do
  CODE=$(curl -sS -o /dev/null -w '%{http_code}' -I -m 20 -L "$URL" 2>/dev/null)
  if [ -z "$CODE" ] || [ "$CODE" = "000" ]; then
    printf '  UNREACHABLE  %s\n' "$URL"
  else
    printf '  HTTP %-4s    %s\n' "$CODE" "$URL"
  fi
done
# Size of the tarball, if the host answered. Worth knowing before committing the
# disk, and it is the same request, so it is free.
SZ=$(curl -sS -I -m 20 -L "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar" 2>/dev/null \
     | awk 'BEGIN{IGNORECASE=1} /^content-length:/ {print $2}' | tr -d '\r' | tail -1)
if [ -n "${SZ:-}" ]; then
  printf '  parameter tarball: %s bytes (%.1f GB)\n' "$SZ" "$(echo "$SZ" | awk '{print $1/1073741824}')"
else
  printf '  parameter tarball size: not reported\n'
fi

# ------------------------------------------------------- B. parameters on disk
say "B. AF2 parameters already on disk"
# BindCraft needs the *multimer* weights. Our scoring bundle may only carry
# monomer/ptm ones, so look for both and name which is which.
FOUND=0
while IFS= read -r f; do
  FOUND=1
  printf '  %s  %s\n' "$(du -h "$f" 2>/dev/null | cut -f1)" "$f"
done < <(find "$HOME" -maxdepth 6 \( -name 'params_model_*multimer*.npz' -o -name 'params_model_*.npz' \) 2>/dev/null | head -30)
[ "$FOUND" = 0 ] && printf '  none found under %s (searched 6 levels deep)\n' "$HOME"
printf '  -- any line above containing "multimer" means assumption B is already satisfied --\n'

# ------------------------------------------------------------- C. containers
say "C. containers available"
ls -lh /home/hpc_apps/*.sqsh 2>/dev/null | awk '{printf "  %s  %s\n", $5, $9}' \
  || printf '  cannot list /home/hpc_apps\n'

say "C2. what our existing venvs already carry"
# Read the dist-info directory names rather than running pip: no container
# start, no GPU, and it works even for a venv built for another python.
# ColabDesign needs JAX; PyRosetta is the dependency we least want to rebuild.
shopt -s nullglob
for V in "$RFD"/venv*/ "$RFD"/*env*/ ; do
  [ -d "$V" ] || continue
  printf '  %s\n' "$V"
  for SP in "$V"lib/python*/site-packages ; do
    [ -d "$SP" ] || continue
    for PKG in jax jaxlib dm_haiku dm-haiku pyrosetta colabdesign torch numpy ; do
      for D in "$SP"/"$PKG"-*.dist-info "$SP"/"$PKG"*.egg-info ; do
        [ -e "$D" ] || continue
        printf '      %s\n' "$(basename "$D")"
      done
    done
  done
done
printf '  -- a venv listing both a jax-* and a pyrosetta-* line is the install target --\n'

# ------------------------------------------------------------------- D. disk
say "D. room"
df -h "$HOME" 2>/dev/null | tail -1 | awk '{printf "  %s: %s used of %s, %s free\n", $6, $3, $2, $4}'
printf '  need: ~5.3 GB parameters + trajectory output\n'

say "E. is anything already here"
for P in "$RFD/BindCraft" "$RFD/colabdesign" "$RFD/bindcraft"; do
  if [ -e "$P" ]; then printf '  EXISTS  %s\n' "$P"; else printf '  absent  %s\n' "$P"; fi
done

say "verdict inputs"
cat <<'TXT'
  Paste all of the above back. The decision it drives:

    A unreachable                -> stop. The parameters cannot be fetched from
                                    this cluster at all and the whole comparison
                                    is off until that changes.
    B shows a "multimer" file    -> skip the 5.3 GB download entirely.
    C2 shows jax + pyrosetta in
      one venv                   -> layer ColabDesign on it, no conda needed.
    C2 shows neither             -> one venv to build from scratch, which is the
                                    expensive branch; cap it at two rounds.

  Nothing here commits a GPU hour, so a bad answer costs only this paste.
TXT
