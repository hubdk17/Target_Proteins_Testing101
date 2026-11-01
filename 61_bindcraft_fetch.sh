#!/bin/bash
# BindCraft step 2: fetch everything, install nothing. LOGIN NODE ONLY.
#
#   cd ~/rfd && bash scripts/61_bindcraft_fetch.sh
#
# Why this runs on the login node and not in a job: compute nodes are
# host-filtered, so the 5.2 GB parameter tarball and both git clones have to be
# pulled here. 60_bindcraft_probe.sh confirmed this node reaches all three hosts
# (HTTP 200) and that we hold only monomer params - params_model_1_ptm.npz -
# so the multimer weights are a genuine download, not a path problem.
#
# What it does, all of it additive - it creates new paths and edits nothing that
# exists:
#   1. clones BindCraft
#   2. clones ColabDesign, so the later pip install needs no network from a
#      compute node
#   3. downloads and unpacks the AF2 parameter tarball into params/
#   4. copies venv_af2 to venv_bc
#
# Step 4 deserves its reason stated: venv_af2 is our working AF2 scoring
# environment, and its contents were hard won - jax 0.4.20 against
# jaxlib 0.4.20+cuda12.cudnn89, dm-haiku pinned at 0.0.13 because 0.0.14 reads
# a jax attribute that 0.4.20 does not have, numpy held at 1.26.4. Installing
# BindCraft's dependencies into that venv risks resolving any one of those pins
# away and breaking stage-3 scoring, which we would discover later and blame on
# something else. A copy costs disk we have (5.3 TB free) and makes that class
# of accident impossible.

set -euo pipefail

RFD="$HOME/rfd"
PARAMS_URL="https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"
EXPECT_BYTES=5587968000     # what the probe measured, to the byte
cd "$RFD"

say() { printf '\n== %s ==\n' "$1"; }

# --------------------------------------------------------------- 1. BindCraft
say "1. BindCraft source"
if [ -d "$RFD/BindCraft/.git" ]; then
  echo "  already cloned; leaving it alone"
else
  git clone --depth 1 https://github.com/martinpacesa/BindCraft.git "$RFD/BindCraft"
fi
git -C "$RFD/BindCraft" log -1 --format='  at %h %ad  %s' --date=short

# ------------------------------------------------------------- 2. ColabDesign
say "2. ColabDesign source"
# BindCraft drives AF2 through ColabDesign. Cloned rather than pip-installed
# from the network, so the install step can run offline inside the container.
if [ -d "$RFD/colabdesign/.git" ]; then
  echo "  already cloned; leaving it alone"
else
  git clone --depth 1 https://github.com/sokrypton/ColabDesign.git "$RFD/colabdesign"
fi
git -C "$RFD/colabdesign" log -1 --format='  at %h %ad  %s' --date=short

# ----------------------------------------------------------- 3. AF2 weights
say "3. AF2 parameters (5.2 GB)"
PDIR="$RFD/BindCraft/params"
mkdir -p "$PDIR"
TAR="$PDIR/alphafold_params_2022-12-06.tar"

if ls "$PDIR"/params_model_*multimer*.npz >/dev/null 2>&1; then
  echo "  multimer params already unpacked; skipping the download"
else
  # -C resumes a part-finished file, so an interrupted download costs only the
  # remainder. --retry rides out a blip rather than failing the whole step.
  echo "  downloading to $TAR  (minutes, not seconds)"
  curl -L --retry 5 --retry-delay 5 -C - -o "$TAR" "$PARAMS_URL"

  # Verify what ARRIVED, not what was asked for. A truncated tarball unpacks
  # partially and the failure would surface much later as a missing model.
  GOT=$(stat -c %s "$TAR")
  printf '  on disk: %s bytes   expected: %s\n' "$GOT" "$EXPECT_BYTES"
  if [ "$GOT" != "$EXPECT_BYTES" ]; then
    echo "  *** SIZE MISMATCH - not unpacking a partial tarball. Re-run this script;"
    echo "      curl -C - will resume from $GOT bytes."
    exit 1
  fi
  echo "  unpacking"
  tar -xf "$TAR" -C "$PDIR"
  rm -f "$TAR"          # 5.2 GB of nothing once unpacked
fi

echo "  parameter files now present:"
ls -1 "$PDIR" | sed 's/^/    /'
NMULTI=$(ls -1 "$PDIR" | grep -c multimer || true)
printf '  multimer files: %s  (BindCraft needs these; 0 means the fetch failed)\n' "$NMULTI"

# -------------------------------------------------------- 4. the venv copy
say "4. venv_bc (a copy of venv_af2)"
if [ -d "$RFD/venv_bc" ]; then
  echo "  already exists; leaving it alone"
else
  test -d "$RFD/venv_af2" || { echo "  venv_af2 missing - wrong cluster state, stop here"; exit 1; }
  echo "  copying (a few GB, PyRosetta is most of it)"
  cp -a "$RFD/venv_af2" "$RFD/venv_bc"
  # A venv records its own absolute path in the activate scripts and in every
  # console-script shebang. Repoint them, or activating venv_bc would quietly
  # put venv_af2 on PATH and we would be installing into the env we set out to
  # protect.
  grep -rl --binary-files=without-match "$RFD/venv_af2" "$RFD/venv_bc/bin" 2>/dev/null \
    | while read -r f; do sed -i "s|$RFD/venv_af2|$RFD/venv_bc|g" "$f"; done
  echo "  leftover references to venv_af2 inside venv_bc/bin (want 0):"
  printf '    %s\n' "$(grep -rl --binary-files=without-match "$RFD/venv_af2" "$RFD/venv_bc/bin" 2>/dev/null | wc -l)"
fi
du -sh "$RFD/venv_bc" 2>/dev/null | sed 's/^/  /'

say "done - nothing installed yet"
cat <<'TXT'
  Fetched only. No pip has run, no GPU was used, and venv_af2 is untouched.
  Paste this output back; the next step is a no-GPU job inside the container
  that tries the imports BindCraft needs and prints the ones that are missing.
  That list, not a guess, decides what we install.
TXT
