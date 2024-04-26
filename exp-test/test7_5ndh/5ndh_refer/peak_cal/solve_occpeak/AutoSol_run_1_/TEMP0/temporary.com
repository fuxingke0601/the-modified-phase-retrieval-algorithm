#!/bin/sh 
phenix.solve <<EOD > AutoSol_run_1_/TEMP0/solve.logit
access_file mirima.dat
WORKDIR AutoSol_run_1_/TEMP0/
LOGFILE solve.log
SOLVEFILE dataset_1_scale.log
SYMFILE p212121.sym
RESOLUTION 1000.0 1.81
CELL  48.7270011902 57.8250007629 100.861999512 90.0 90.0 90.0
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

MAD_ATOM BR
LAMBDA 1
RAWMADFILE ANISO_5ndh-sf_PHX.sca
ATOMNAME BR
WAVELENGTH 0.9201
FPRIMV_MAD -9.94
FPRPRV_MAD 3.82
NRES 64
NANOMALOUS 6
SKIP_SOLVE
SAD  


EOD
