#!/bin/bash

#c=$(ls prasa_*2_2.61_unique.mtz)
#echo $c

if [ ! -e pasa_*2_2.61_unique.mtz ];then
	echo "bad"
fi

#shopt -s nullglob

#ls prasa*2.61_unique.mtz
#c=($(ls prs*2.61_unique.mtz))
#echo ${c[*]}
#if [ ${#c[@]} -eq 0 ]; then
#    echo "bad"
#fi

reso=4.01
if (( $(echo "${reso} > 3.8" | bc) ));then
   nitr=748
else
   nitr=498
fi
echo $nitr
