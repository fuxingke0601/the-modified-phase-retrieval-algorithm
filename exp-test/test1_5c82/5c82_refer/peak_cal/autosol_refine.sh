#!/bin/bash



if [ $# != 7 ];then
echo "\$1 initial mtz, including I+, I-, OR F+, F-"
echo "\$2 the sequence of protein"
echo "\$3 the pdb of heavy atom"
echo "\$4 the type of heavy atom"
echo "\$5 the wavelength of X-ray"
echo "\$6 the solvent_fraction"
echo "\$7 the ncs_copies"
echo "eg: $0 3k9g-sf.mtz 3K9G.fasta 3k9g_12I.pdb I 1.54180 0.5145 1"
echo "运行时注意sca文件是否使用缩放因子，查看sca中是否有****"
exit
fi

mtzfile=$1
seqfile=$2
pdbha_file=$3
atom_type=$4
lambda=$5
solvent_fraction=$6
ncs_copies=$7

echo "$0 $1 $2 $3 $4 $5 $6 $7" >autosol_run_log.txt

filetitle=$(basename $mtzfile | sed 's/...[az]$//g')
column_name=$( mtzdmp ${mtzfile} | grep -A2 " * Column Labels :" |sed -n "3,3 p")
echo "$column_name"
letter1="I(+)"
letter2="I(-)"

###############
if echo $column_name | grep -q ${letter1} && echo $column_name | grep -q ${letter2} ; then
mtz2various HKLIN "$mtzfile" HKLOUT "${filetitle}.sca" >> autosol_run_log.txt <<END
 OUTPUT SCALEPACK
labin I(+)=I(+) SIGI(+)=SIGI(+) I(-)=I(-) SIGI(-)=SIGI(-)
#SCALE 0.01
end
END
########
else
########
mtz2various HKLIN "$mtzfile" HKLOUT "${filetitle}.sca" <<END
 OUTPUT SCALEPACK
labin F(+)=F(+) SIGF(+)=SIGF(+) F(-)=F(-) SIGF(-)=SIGF(-)
FSQUARED
#SCALE 0.01
end
END

fi
###################
num=$( grep 'HETATM\|ATOM' ${pdbha_file} | wc -l ) #重原子数目
phenix.autosol data=${filetitle}.sca   atom_type=${atom_type} lambda=${lambda} solvent_fraction=${solvent_fraction}  nproc=4 seq_file=${seqfile}   sites=${num} sites_file=${pdbha_file} ncs_copies=${ncs_copies}
