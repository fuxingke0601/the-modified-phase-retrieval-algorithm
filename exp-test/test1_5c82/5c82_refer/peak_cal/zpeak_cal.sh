#!/bin/bash

#从CC最好的mtz中，使用peakmax得到重原子坐标，根据峰值设定占位率
#使用bp3进行相位修正

if [ $# != 6 ];then
echo "\$1 the mtz file of heavy atom after CF"
echo "\$2 the numbers of searched heavy atom"
echo "\$3 the type of heavy atom "
echo "\$4 whether to refine the coordinates of heavy atom"
echo "\$5 the wavelength/A"
echo "\$6 initial mtz, including F(+),SIGF(+),F(-),SIGF(-)"
echo "eg: $0 recover_211.mtz 17 SE 1 0.9785 rcaoy_truncate.mtz"
exit
fi

bestCC_mtz=$1
num=$2
atom_type=$3
flag=$4       #是否使用BP3修正重原子
wavelength=$5
mtz_file=$6

 fft HKLIN ${bestCC_mtz} MAPOUT bestCC_map.map >peak_cal_log.txt<<END
title fft_to_map
xyzlim asu
scale F1 1.0
labin F1=FC PHI=PHIC
end
END

 peakmax MAPIN bestCC_map.map XYZOUT pdb_CF.pdb XYZFRC peaks_ha.txt >>peak_cal_log.txt<<END
threshhold rms 1.2
numpeaks $num
output brookhaven frac
bfactor 20.0 1.0
residue HAT
atname ${atom_type}
chain X
END

awk '{if(NR>2) max=($6>max)?$6:max} END{OFS="\t";}{if(NR>2) $6=$6/max;if(NR>2) print $6}' peaks_ha.txt | awk -v atom="${atom_type}" 'NR==FNR{a[NR]=$1; nr=NR} NR>FNR{if (FNR<=4) print $0;
                           else if (NF>1) printf "%s %6d %3s %4s %s %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11s\n","ATOM",$6,$3,$4,$5,$6,$7,$8,$9,a[FNR-4],20,atom;
                           else print $0}' - pdb_CF.pdb >pdb_CF_format.pdb

######################################################################################
###################使用bp3进行重原子坐标修正##########################################
######################################################################################
if [ ${flag} == 1 ];then
   
   #判断mtz文件中是否有F+,F-,如果没有用ctruncate计算
   column_type=$( mtzdmp ${mtz_file} | grep -A2 "* Column Types :" |sed -n "3,3 p") #得到mtz文件中每一列的数据类型
   letter1="G"
   letter2="L"
   letter3="K"
   letter4="M"
   if echo $column_type | grep -q ${letter1} && echo $column_type | grep -q ${letter2} ; then
       echo "mtz file have F+,F-" >ctruncate.log
   else
       ctruncate -mtzin ${mtz_file} -mtzout "ctruncate.mtz"  -colano "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]" >ctruncate.log
       mtz_file="ctruncate.mtz"
    fi


   fprime1=$(gemmi fprime --wavelength=${wavelength} ${atom_type} |awk '{if(NR==2) printf $4}')
   fprime2=$(gemmi fprime --wavelength=${wavelength} ${atom_type} |awk '{if(NR==2) printf $5}' )

   #echo "$fprime1 $fprime2"

 bp3 HKLIN ${mtz_file} HKLOUT Bp3.mtz << END
TARG SAD
XTAL bp3_phase
DNAM bp3_phase
COLU F+=F(+) SF+=SIGF(+) F-=F(-) SF-=SIGF(-)
FORM SE FP=${fprime1} FPP=${fprime2}
MODL pdb_CF.pdb
LABO FB=BP3_FB PHIB=BP3_PHIB FOM=BP3_FOM HLA=BP3_HLA HLB=BP3_HLB HLC=BP3_HLC HLD=BP3_HLD FDIFF=BP3_DELFAN PDIFF=BP3_PHDELAN
END

grep "HETATM" heavy-1.pdb |awk 'NR==FNR{a[NR]=$6; b[NR]=$7; c[NR]=$8; o[NR]=$9; bb[NR]=$10; nr=NR} NR>FNR{if (FNR<=4) print $0;
                           else if (NF>1) printf "%s %6d %3s %4s %s %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11s\n","ATOM",$2,$3,$4,$5,$6,a[FNR-4],b[FNR-4],c[FNR-4],o[FNR-4],bb[FNR-4],$12;
                           else print $0}' - pdb_CF_format.pdb >heavy-1_format.pdb

grep "HETATM" heavy-oh-1.pdb |awk 'NR==FNR{a[NR]=$6; b[NR]=$7; c[NR]=$8; o[NR]=$9; bb[NR]=$10; nr=NR} NR>FNR{if (FNR<=4) print $0;
                           else if (NF>1) printf "%s %6d %3s %4s %s %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11s\n","ATOM",$2,$3,$4,$5,$6,a[FNR-4],b[FNR-4],c[FNR-4],o[FNR-4],bb[FNR-4],$12;
                           else print $0}' - pdb_CF_format.pdb >heavy-oh-1_format.pdb


fi
