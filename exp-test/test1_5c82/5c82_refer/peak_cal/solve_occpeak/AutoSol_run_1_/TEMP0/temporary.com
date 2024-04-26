#!/bin/sh 
phenix.solve <<EOD > AutoSol_run_1_/TEMP0/solve.logit
access_file mirima.dat
WORKDIR AutoSol_run_1_/TEMP0/
LOGFILE solve.log
SOLVEFILE dataset_1_scale.log
SYMFILE c2.sym
RESOLUTION 1000.0 2.2
CELL  96.5360031128 44.6409988403 42.2439994812 90.0 97.6650009155 90.0
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
RAWMADFILE ANISO_5c82-sf_PHX.sca
ATOMNAME SE
WAVELENGTH 0.97869
FPRIMV_MAD -8.0
FPRPRV_MAD 4.5
NRES 190
NANOMALOUS 6
SKIP_SOLVE
SAD  


EOD
