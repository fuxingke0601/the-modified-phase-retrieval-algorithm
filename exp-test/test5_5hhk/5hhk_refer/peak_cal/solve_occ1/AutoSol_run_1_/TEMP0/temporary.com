#!/bin/sh 
phenix.solve <<EOD > AutoSol_run_1_/TEMP0/solve.logit
access_file mirima.dat
WORKDIR AutoSol_run_1_/TEMP0/
LOGFILE solve.log
SOLVEFILE dataset_1_scale.log
SYMFILE p21.sym
RESOLUTION 1000.0 1.4
CELL  42.0099983215 88.1200027466 79.6699981689 90.0 95.4489974976 90.0
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
RAWMADFILE ANISO_5hhk-sf_PHX.sca
ATOMNAME SE
WAVELENGTH 0.97918
FPRIMV_MAD -8.0
FPRPRV_MAD 4.5
NRES 610
NANOMALOUS 35
SKIP_SOLVE
SAD  


EOD
