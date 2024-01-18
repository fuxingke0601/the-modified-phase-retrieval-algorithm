#!/bin/bash

../../bin/CF HKLIN 1542463.mtz MAPOUT recovered.map HKLOUT recovered.mtz\
<<eof-DM
TITL charge flipping
LABI F=Fobs
FRAC 0.12
NITR 200
NTRY 20
WILS
NWLS 50
RWLS 3.0
CONS C 48 H 40 Mo 4 O 16 P 4
#ECAL
#NECA 100
#RESO 1.0
END
eof-DM
