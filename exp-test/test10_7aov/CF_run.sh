#!/bin/bash
mtz=$1
weak=$2
scut=$3
cctr=$4
reso=$5
say=1
say=${6:-$say}
#echo "$say"
ntry=$7
nitr=$8
abso_path="."
abso_path=${9:-$abso_path}


path=$(readlink -f ${abso_path})
#echo "$path"
if [[ $# != 8 && $# != 9 ]];then
        echo "\$1 means mtz"
        echo "\$2 means the percentage of weak reflection in observing reflection"
	echo "\$3 means the number of strong reflection in direct method"
        echo "\$4 means the threshold of calc CC"
        echo "\$5 means the resolution"
	echo "\$6 means the resolution"
	echo "\$7 means the number of different initial random phases"
	echo "\$8  means the number of iteration"
	echo "\$9 path setting"
        echo "e.g.:$0 prasa_5oq2_2.61_unique.mtz 0.5 1300 0.6217 2.61 0 400 498"
        exit
fi
#scp -r recovered.mtz recovered.map  fuxingke@fudemac:/Users/fuxingke/linux/fxk/heavy_atom_search/geng_talk/charge-flipping/exp-test/coot
#scp -r recovered.mtz recovered.map  fuxingke@fudemac:/Users/fuxingke/Desktop/coot
#prasa_inpmerge_unique1.mtz LABI F=EA SIGF=SIGEA  prasa_5c82_2.86_unique.mtz F=EA SIGF=SIGEA FA=FA SIGFA=SIGFA F=ecalculate.E_sigE.E 5c82_4se_sfall_ecal_unique1.mtz
${path}/../../bin/CF HKLIN  $mtz  MAPOUT recovered_.map HKLOUT recovered_.mtz\
<<eof-DM
TITL charge flipping
LABI F=EA SIGF=SIGEA FA=FA SIGFA=SIGFA
FRAC 0.13
THRE 0.0
SAY $say
WEAK $weak
SCUT $scut
CCTR $cctr
NITR $nitr
NTRY $ntry
#WILS
#NWLS 50
#RWLS 3.5
CONS SE 8
#ECAL
#NECA 100
RESO $reso
END
eof-DM
