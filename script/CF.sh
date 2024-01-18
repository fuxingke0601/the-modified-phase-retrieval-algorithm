#!/bin/bash

../bin/CF HKLIN 1000006.mtz MAPOUT recovered.map HKLOUT recovered.mtz\
<<eof-DM
TITL charge flipping
LABI F=F
FRAC 0.12
NITR 200
NTRY 20
WILS
NWLS 50
RWLS 3.0
CONS C 22 H 25 Cl 1 N 2 O 8
#ECAL
#NECA 100
#RESO 1.0
END
eof-DM
