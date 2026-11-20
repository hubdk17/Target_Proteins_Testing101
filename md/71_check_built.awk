# Two build gates on the tleap-built system, in awk because the nvhpc
# container has no python. Exits 1 if either fails, so the job aborts BEFORE
# any GPU time is spent.
#
#   1. CHIRALITY. tleap adds hydrogens and rebuilds side chains, and Amber
#      protein force fields carry no improper term on the CA centre, so a
#      D-residue is invisible to energy, RMSD, Rg, DSSP and every contact
#      metric. Gate: the N-C-CA-CB improper is ~ -122 deg for L (calibrated on
#      1YCR, 94/94 at -122.3 +/- 3.0); positive means D.
#   2. CHAIN SEPARATION. At least one C(i)-N(i+1) distance must exceed 2.0 A,
#      i.e. target and binder are not one covalent chain. A fused complex can
#      never come apart and every downstream metric still reads as stable.
#
# Usage: awk -f 71_check_built.awk sys_built.pdb
/^ATOM/ {
  rn = substr($0,18,3); gsub(/ /,"",rn)
  if (rn=="WAT"||rn=="HOH"||rn=="NA"||rn=="CL"||rn=="Na+"||rn=="Cl-") next
  key = substr($0,22,5)
  if (!(key in seen)) { seen[key]=1; n_res++; ord[n_res]=key; name[key]=rn }
  an = substr($0,13,4); gsub(/ /,"",an)
  if (an=="N"||an=="CA"||an=="C"||an=="CB") {
    x[key,an]=substr($0,31,8)+0; y[key,an]=substr($0,39,8)+0; z[key,an]=substr($0,47,8)+0
    has[key,an]=1
  }
}
END {
  # ---- gate 1: chirality -------------------------------------------------
  nchk=0; nd=0; bad=""
  for (i=1;i<=n_res;i++) {
    k=ord[i]
    if (name[k]=="GLY") continue
    if (!(has[k,"N"] && has[k,"CA"] && has[k,"C"] && has[k,"CB"])) continue
    # b1 = C-N, b2 = CA-C, b3 = CB-CA   (improper N-C-CA-CB)
    b1x=x[k,"C"]-x[k,"N"];  b1y=y[k,"C"]-y[k,"N"];  b1z=z[k,"C"]-z[k,"N"]
    b2x=x[k,"CA"]-x[k,"C"]; b2y=y[k,"CA"]-y[k,"C"]; b2z=z[k,"CA"]-z[k,"C"]
    b3x=x[k,"CB"]-x[k,"CA"];b3y=y[k,"CB"]-y[k,"CA"];b3z=z[k,"CB"]-z[k,"CA"]
    n1x=b1y*b2z-b1z*b2y; n1y=b1z*b2x-b1x*b2z; n1z=b1x*b2y-b1y*b2x
    n2x=b2y*b3z-b2z*b3y; n2y=b2z*b3x-b2x*b3z; n2z=b2x*b3y-b2y*b3x
    nb=sqrt(b2x*b2x+b2y*b2y+b2z*b2z); if (nb==0) continue
    u1=b2x/nb; u2=b2y/nb; u3=b2z/nb
    mx=n1y*u3-n1z*u2; my=n1z*u1-n1x*u3; mz=n1x*u2-n1y*u1
    v=atan2(mx*n2x+my*n2y+mz*n2z, n1x*n2x+n1y*n2y+n1z*n2z)*57.29577951308232
    nchk++
    if (v>0) { nd++; if (nd<=10) bad=bad " " name[k] k }
  }
  printf "chiral residues checked: %d  D-residues: %d%s\n", nchk, nd, bad
  # ---- gate 2: chain separation ------------------------------------------
  nbrk=0; at=""
  for (i=1;i<n_res;i++) {
    a=ord[i]; b=ord[i+1]
    if (!(has[a,"C"] && has[b,"N"])) continue
    dx=x[a,"C"]-x[b,"N"]; dy=y[a,"C"]-y[b,"N"]; dz=z[a,"C"]-z[b,"N"]
    if (sqrt(dx*dx+dy*dy+dz*dz)>2.0) { nbrk++; if (nbrk<=6) at=at " " (i+1) }
  }
  printf "protein residues: %d  chain breaks found: %d  at built-residue index(es):%s\n", n_res, nbrk, at
  if (nd>0)   { print "D-RESIDUES FOUND in built system -- aborting before MD"; exit 1 }
  if (nbrk==0){ print "NO CHAIN BREAK -- target and binder were fused into one chain by tleap"; \
                print "Fix: make sure a TER record separates the chains in clean.pdb, then rebuild"; exit 1 }
  print "both build gates PASSED"
}
