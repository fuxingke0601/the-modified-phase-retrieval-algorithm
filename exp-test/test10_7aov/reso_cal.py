import numpy as np
import gemmi
import sys

'''考虑了衍射点数目，既不能超过25000，又不能少于3000，还兼顾了d’‘/sig'''
#file_path="/Users/fuxingke/Desktop/caoy_dual_space/prasa_cop/shelxC.log"
file_path=sys.argv[1]
with open(file_path) as file:
    file_contents=file.readlines()

keyword1="Resl."
i=0
for line in file_contents:
    if keyword1 in line:
        #print(line,i)
        break
    i=i+1
#print(file_contents[i],file_contents[i+7])
reso_string=file_contents[i]
dano_sigma_string=file_contents[i+7]
# number_refl_string=file_contents[i+1]
# 使用字符串处理函数将数字提取出来
reso = [float(num) for num in reso_string.split() if num.replace('.', '', 1).isdigit()]
dano_sigma=[float(num) for num in dano_sigma_string.split() if num.replace('.', '', 1).isdigit()]
# number_refl_array=[int(num) for num in number_refl_string.split() if num.replace('.', '', 1).isdigit()]
reso_array=np.array(reso)[::-1]
dano_sigma_array=np.array(dano_sigma)[::-1]

i=0
cut_1_2=0
for dano in dano_sigma_array:
    if dano>=1.2:
        cut_1_2=i
        break
    i+=1
i=0
cut_1_0=0
for dano in dano_sigma_array:
    if dano>=1.0:
        cut_1_0=i
        break
    i+=1
i=0
cut_0_8=0
for dano in dano_sigma_array:
    if dano>=0.8:
        cut_0_8=i
        break
    i+=1
index_dano=np.max([cut_1_2,cut_0_8,cut_1_0])
if index_dano !=0:
    #reso1=reso_array[index-1]
    index_dano=index_dano-1
else:
    #reso1=reso_array[0]
    index_dano=0

#加上衍射点数目判定
mtz_file=sys.argv[2]
mtz_gemmi=gemmi.read_mtz_file(mtz_file)
mtz_array=np.array(mtz_gemmi,copy=True)
hkl_array = mtz_array[:,0:3]
d_spacing = mtz_gemmi.cell.calculate_d_array(hkl_array.astype(int))
d_spacing=np.array(d_spacing)
num_refl_array = np.sum(d_spacing > reso_array[:, np.newaxis], axis=1)

reso2=0
i=0
index_num=0
numbers=np.arange(num_refl_array.shape[0])[::-1]
#print(index)
for dano in numbers:
    if num_refl_array[dano]>25000:
        #reso2=reso_array[dano+1]
        index_num=dano+1
        break
    i+=1

#reso=np.max([])
index=index_dano
if index_num>0:
    #限制衍射点数目小于25000
    index=np.max([index_dano,index_num])
if num_refl_array[index]<=3000:
    if num_refl_array[cut_1_0]<=3000:
        index=cut_0_8
    else:
        index=cut_1_0
reso=reso_array[index]

print(reso)


