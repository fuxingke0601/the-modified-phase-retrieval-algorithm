#!/bin/bash
s=$1
#scp -r recovered.mtz recovered.map  fuxingke@fudemac:/Users/fuxingke/linux/fxk/heavy_atom_search/geng_talk/charge-flipping/exp-test/coot
#scp -r recovered.mtz recovered.map  fuxingke@fudemac:/Users/fuxingke/Desktop/coot
#prasa_inpmerge_unique1.mtz LABI F=EA SIGF=SIGEA  prasa_5c82_2.86_unique.mtz F=EA SIGF=SIGEA FA=FA SIGFA=SIGFA F=ecalculate.E_sigE.E 5c82_4se_sfall_ecal_unique1.mtz
../../bin/CF HKLIN  prasa_5oq2_2.61_unique.mtz  MAPOUT recovered_$s.map HKLOUT recovered_$s.mtz\
<<eof-DM
TITL charge flipping
LABI F=EA SIGF=SIGEA FA=FA SIGFA=SIGFA
FRAC 0.13
THRE 0.0
SAY 1
WEAK 0.28
SCUT 1500
CCTR 0.6217
NITR 498
NTRY 400
#WILS
#NWLS 50
#RWLS 3.5
CONS SE 8
#ECAL
#NECA 100
RESO 3.16
END
eof-DM
