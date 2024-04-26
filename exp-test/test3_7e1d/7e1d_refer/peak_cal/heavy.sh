#!/bin/sh

bp3 HKLIN ctruncate.mtz HKLOUT nextrun.mtz << eof > nextrun.log

SITE   1 0.339 0.577 0.061  
SITE   2 0.101 0.697 0.133  
SITE   3 0.137 0.863 0.186  
SITE   4 0.154 0.105 0.122  
SITE   5 0.198 0.084 0.121  
SITE   6 0.111 0.703 0.146  


XTAL bp3_phase
  CELL   35.814   52.762   97.397  90.000  90.000  90.000
  ATOM SE SITE 1
    OCCU   0.838
    BISO   2.500
  ATOM SE SITE 2
    OCCU   0.652
    BISO  13.936
  ATOM SE SITE 3
    OCCU   0.957
    BISO  20.508
  ATOM SE SITE 4
    OCCU   0.542
    BISO  18.231
  ATOM SE SITE 5
    OCCU   0.443
    BISO  21.432
  ATOM SE SITE 6
    OCCU   0.423
    BISO  21.845
  DNAMe bp3_phase
    COLUmn  F+= F(+)  SF+= SIGF(+)  F-= F(-)  SF-= SIGF(-)
    RESOlution    2.00   28.35
    BINS 10
    SDLU  0.916579  0.964786  0.907052  0.791205  0.752520  0.407454  0.536041  0.342736  0.066866  0.069867 
    FORM  SE FP=  -6.436  FPP=   3.836


TARGet SAD

# Refine all parameters
REFAll


eof
