#!/bin/sh

bp3 HKLIN 5ndh-sf.mtz HKLOUT nextrun.mtz << eof > nextrun.log

SITE   1 -0.986 -0.529 -0.197  
SITE   2 -0.175 -0.804 -0.054  
SITE   3 -0.333 -0.670 -0.048  
SITE   4 -0.105 -0.363 -0.208  
SITE   5 0.018 -0.452 -0.196  
SITE   6 -0.035 -0.565 -0.196  


XTAL bp3_phase
  CELL   48.727   57.825  100.862  90.000  90.000  90.000
  ATOM BR SITE 1
    OCCU   4.129
    BISO   9.868
  ATOM BR SITE 2
    OCCU   8.962
    BISO  69.174
  ATOM BR SITE 3
    OCCU   0.998
    BISO  18.419
  ATOM BR SITE 4
    OCCU   1.455
    BISO  31.680
  ATOM BR SITE 5
    OCCU   1.001
    BISO  24.662
  ATOM BR SITE 6
    OCCU   1.354
    BISO   4.700
  DNAMe bp3_phase
    COLUmn  F+= F(+)  SF+= SIGF(+)  F-= F(-)  SF-= SIGF(-)
    RESOlution    1.81   50.17
    BINS 13
    SDLU  0.282740  0.371238  0.533691  0.594276  0.404313  0.240816  0.149508  0.087847  0.001481  0.000000  0.001470  0.025690  0.002887 
    FORM  SE FP=  -8.831  FPP=   3.823
    FORM  BR FP=  -0.767  FPP=   1.283


TARGet SAD

# Refine all parameters
REFAll


eof
