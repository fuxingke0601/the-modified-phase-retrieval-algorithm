#!/bin/sh 
phenix.solve <<EOD > AutoSol_run_1_/TEMP0/solve.logit
access_file mirima.dat
WORKDIR AutoSol_run_1_/TEMP0/
LOGFILE solve.log
SOLVEFILE dataset_1_scale.log
SYMFILE p21.sym
RESOLUTION 1000.0 1.85
CELL  87.8109970093 58.6749992371 95.7300033569 90.0 114.309997559 90.0
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
RAWMADFILE ANISO_4yf1-sf_PHX.sca
ATOMNAME SE
WAVELENGTH 0.97856
FPRIMV_MAD -8.0
FPRPRV_MAD 4.5
NRES 748
NANOMALOUS 10
SKIP_SOLVE
SAD  


EOD
