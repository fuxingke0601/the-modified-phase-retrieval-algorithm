#!/bin/sh 
phenix.solve <<EOD > AutoSol_run_1_/TEMP0/solve.logit
access_file mirima.dat
WORKDIR AutoSol_run_1_/TEMP0/
LOGFILE solve.log
SOLVEFILE dataset_1_scale.log
SYMFILE p41212.sym
RESOLUTION 1000.0 1.55
CELL  57.6150016785 57.6150016785 150.807998657 90.0 90.0 90.0
readdenzo
premerged
nequiv_separate 1
require_nat

ratio_out 10.0

ratmin 0.0

ikeepflag 1

id_scale_ref 1

projectname project

crystalname crystal

datasetname dataset

MAD_ATOM SE
LAMBDA 1
RAWMADFILE ANISO_5t3g-sf_PHX.sca
ATOMNAME SE
WAVELENGTH 0.978
FPRIMV_MAD -8.0
FPRPRV_MAD 4.5
NRES 207
NANOMALOUS 15
SKIP_SOLVE
SAD  


EOD
