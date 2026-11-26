# Minimal stand-in for pdb4amber: the nvhpc container has NO python at all, so
# pdb4amber (a #!/usr/bin/python3 script) cannot run and the whole build died
# at its first line. Our inputs are narrow and known -- single model, heavy
# atoms only, no HETATM, no altlocs, standard residue names, two chains -- so
# the cleaning tleap actually needs is: keep ATOM records, drop hydrogens and
# waters/ions, drop altloc duplicates, and INSERT A TER AT EVERY CHAIN CHANGE.
# That TER is the whole point: without it tleap bonds the last target residue
# to the first binder residue (see md-chain-fusion-trap).
# Usage: awk -f 71_clean_pdb.awk in.pdb > clean.pdb
/^ENDMDL/ { exit }                       # first model only
!/^ATOM/  { next }                       # drops HETATM, ANISOU, CONECT, ...
{
  rn = substr($0,18,3); gsub(/ /,"",rn)
  if (rn=="WAT" || rn=="HOH" || rn=="NA" || rn=="CL" || rn=="SO4" || rn=="GOL") next
  alt = substr($0,17,1)
  if (alt!=" " && alt!="A") next         # keep conformer A only
  an = substr($0,13,4); gsub(/ /,"",an)
  el = substr($0,77,2); gsub(/ /,"",el)
  if (el=="H" || el=="D") next           # no hydrogens: tleap adds them
  if (el=="" && an ~ /^[0-9]?H/) next    # fallback when element column is blank
  ch = substr($0,22,1)
  if (NR_out>0 && ch!=prev_ch) print "TER"
  # blank the altloc so tleap never sees it, keep everything else verbatim
  print substr($0,1,16) " " substr($0,18)
  NR_out++; prev_ch = ch
}
END { if (NR_out>0) print "TER"; print "END" }
