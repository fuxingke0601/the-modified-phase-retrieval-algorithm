import  numpy as np
import matplotlib.pyplot as plt
import sys

filename=sys.argv[1]
'''
with open(filename) as file:
    skew_str=file.readlines()

skew_array=np.array(skew_str).astype(float)

plt.plot(skew_array)
plt.show()
'''
with open(filename) as file:
    str_array=file.readlines()

skew_array=[]
cc9_array=[]

for str in str_array:
    # skew_array=str.split
    skew_str=str.split()[0]
    cc9_str=str.split()[1]
#    if skew_str.replace('.', '', 1).lstrip('-+').isdigit() & cc9_str.replace('.', '', 1).lstrip('-+').isdigit():
    if skew_str.replace('.', '', 1).lstrip('-+').isdigit():
        skew_array.append(float(skew_str))
        cc9_array.append(float(cc9_str))

plt.plot(skew_array)
#plt.plot(cc9_array)
plt.show()
