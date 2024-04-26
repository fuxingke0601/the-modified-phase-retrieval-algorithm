#!/bin/sh

bp3 HKLIN 5c82-sf.mtz HKLOUT nextrun.mtz << eof > nextrun.log

SITE   1 0.382 0.360 0.061   NOREf  Y 
SITE   2 0.297 0.257 0.053  
SITE   3 0.324 0.171 0.010  
SITE   4 0.310 0.271 0.884  
SITE   5 0.316 0.304 0.006  
SITE   6 0.311 0.355 0.064  


XTAL bp3_phase
  CELL   96.536   44.641   42.244  90.000  97.665  90.000
  ATOM SE SITE 1
    OCCU   1.325
    BISO  10.827
  ATOM SE SITE 2
    OCCU   1.119
    BISO  11.275
  ATOM SE SITE 3
    OCCU   1.253
    BISO  16.339
  ATOM SE SITE 4
    OCCU   1.231
    BISO  21.190
  ATOM SE SITE 5
    OCCU   0.031
    BISO  19.305
  ATOM SE SITE 6
    OCCU   0.049
    BISO  18.723
  DNAMe bp3_phase
    COLUmn  F+= F(+)  SF+= SIGF(+)  F-= F(-)  SF-= SIGF(-)
    RESOlution    2.20   47.84
    BINS 8
    SDLU  0.881947  0.972128  0.956454  0.916702  0.939394  0.864545  0.680077  0.579174 
    FORM  SE FP=  -6.639  FPP=   3.837


TARGet SAD

# Refine all parameters
REFAll


eof
