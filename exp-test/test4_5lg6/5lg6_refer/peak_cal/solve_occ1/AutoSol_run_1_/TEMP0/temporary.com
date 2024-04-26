#!/bin/sh 
phenix.solve <<EOD > AutoSol_run_1_/TEMP0/solve.logit
access_file mirima.dat
WORKDIR AutoSol_run_1_/TEMP0/
LOGFILE solve.log
SOLVEFILE dataset_1_scale.log
SYMFILE p21.sym
RESOLUTION 1000.0 2.5
CELL  66.4649963379 215.690002441 78.6299972534 90.0 91.8990020752 90.0
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
RAWMADFILE ANISO_5lg6-sf_PHX.sca
ATOMNAME SE
WAVELENGTH 0.9795
FPRIMV_MAD -8.0
FPRPRV_MAD 4.5
NRES 1936
NANOMALOUS 50
SKIP_SOLVE
SAD  


EOD
