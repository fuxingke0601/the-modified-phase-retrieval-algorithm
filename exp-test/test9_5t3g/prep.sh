#!/bin/bash
# mtz2various HKLIN /Users/fuxingke/Documents/data/SAD/yao_deqiang/caoy.mtz HKLOUT caoy_MTZ2V.hkl << END
#OUTP shelx
#LABI F(+)=F-obs-filtered(+) SIGF(+)=SIGF-obs-filtered(+) F(-)=F-obs-filtered(-) SIGF(-)=SIGF-obs-filtered(-)
#END
#./a.sh /Users/fuxingke/Documents/data/SAD/pdb_lib/4_8i59_14se_P41212/8i59-sf.mtz
#./a.sh /Users/fuxingke/Documents/data/SAD/pdb_lib/11_7e1d_4se_P212121/7e1d-sf.mtz 


file=$1
column_type=$( mtzdmp ${file} | grep -A2 "* Column Types :" |sed -n "3,3 p") #得到mtz文件中每一列的数据类型
filetitle=$(basename $1 | sed 's/...[az]$//g')
sg=$( mtzdmp ${file} -e| awk -F"'" '/Space group/{print $2}')                 #得到空间群信息
sg_1=$( echo ${sg} |tr -d ' ')
cell=$( mtzdmp ${file} | grep -A2 "* Cell Dimensions :" |sed -n "3,3 p")
reso=$( mtzdmp ${file} -e|grep -A2 "*  Resolution Range :" |sed -n "3,3 p" | awk '{printf $(NF-2)}')
printf "$column_type $filetitle $sg $sg_1 $cell $reso \n"
# 要查找的字母
letter1="G"
letter2="L"
letter3="K"
letter4="M"

# 使用grep命令查找字母是否在文件的某一行中
if echo $column_type | grep -q ${letter1} && echo $column_type | grep -q ${letter2} ; then
    #把F+, F-取出来放在一个新的文件中
#    cad HKLIN1 ${file} HKLOUT ${filetitle}_cad.mtz > ${filetitle}_cad.log << END
#VALM NaN NOOUTPUT
#LABI FILE_NUMBER 1 E1=F(+) E2=SIGF(+) E3=F(-) E4=SIGF(-)
#LABO  FILE_NUMBER 1 E1=F(+) E2=SIGF(+) E3=F(-) E4=SIGF(-)
#END
#    mtz2various HKLIN ${filetitle}_cad.mtz HKLOUT ${filetitle}_ctrun-cad_MTZ2V.hkl > ${filetitle}_ctrun-cad_MTZ2V.log<< END
#OUTP shelx
#LABI F(+)=F(+) SIGF(+)=SIGF(+) F(-)=F(-) SIGF(-)=SIGF(-)
#END
mtz2various HKLIN ${file} HKLOUT ${filetitle}_ctrun-cad_MTZ2V.hkl > ${filetitle}_ctrun-cad_MTZ2V.log<< END
OUTP shelx
LABI F(+)=F(+) SIGF(+)=SIGF(+) F(-)=F(-) SIGF(-)=SIGF(-)
END

else 

    ctruncate -mtzin ${file} -mtzout "${filetitle}_ctruncate.mtz"  -colano "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]" >${filetitle}_ctruncate.log  #把I+, I- 转为F+, F-
    mtz2various HKLIN ${filetitle}_ctruncate.mtz HKLOUT ${filetitle}_ctrun-cad_MTZ2V.hkl >${filetitle}_ctrun-cad_MTZ2V.log << END
OUTP shelx
LABI F(+)=F(+) SIGF(+)=SIGF(+) F(-)=F(-) SIGF(-)=SIGF(-)
END

fi

 shelxc SHELX > shelxC.log<< END
CELL ${cell}
SPAG ${sg}
SFAC SE
SAD -f ${filetitle}_ctrun-cad_MTZ2V.hkl
END

 f2mtz HKLIN SHELX_fa.hkl HKLOUT f2mtz.mtz > f2mtz.log<< END
FORM '(3I4,2F8.2,I4)'
SYMM ${sg_1}
CELL ${cell}
NAME crystal crystal
LABO H K L FA SIGFA ALPHA
CTYP H H H F Q P
END

 ecalc HKLIN f2mtz.mtz HKLOUT shelx_ecal.mtz >shelx_ecal.log<< END
LABI DPH=FA SIGDPH=SIGFA
LABO E=EA SIGE=SIGEA
END

cad HKLIN1 shelx_ecal.mtz HKLOUT shelx_ecal_cad.mtz >shelx_ecal_cad.log << END
VALM NaN NOOUTPUT
LABI FILE_NUMBER 1 E1=FA E2=SIGFA E3=EA E4=SIGEA E5=ALPHA
LABO  FILE_NUMBER 1 E1=FA E2=SIGFA E3=EA E4=SIGEA E5=ALPHA
END

uniqueify shelx_ecal_cad.mtz shelx_${filetitle}_ecal_cad_unique.mtz 


