#!/bin/sh 
phenix.solve <<EOD > AutoSol_run_1_/TEMP0/solve.logit
access_file mirima.dat
WORKDIR AutoSol_run_1_/TEMP0/
LOGFILE solve.log
SOLVEFILE dataset_1_scale.log
SYMFILE p212121.sym
RESOLUTION 1000.0 2.0
CELL  35.813999176 52.7620010376 97.3970031738 90.0 90.0 90.0
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
RAWMADFILE ANISO_7e1d-sf_PHX.sca
ATOMNAME SE
WAVELENGTH 0.97852
FPRIMV_MAD -8.0
FPRPRV_MAD 4.5
NRES 192
NANOMALOUS 6
SKIP_SOLVE
SAD  


EOD
