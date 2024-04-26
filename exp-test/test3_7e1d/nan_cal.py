import numpy as np
import matplotlib.pyplot as plt
import sys

def nan_cal(file_path,search_text):
    with open(file_path) as file:
        file_contents = file.readlines()

    i=0; flag=False
    num_CC=[]
    for line in file_contents:
        if flag:
            line_number=line.strip().split()
            # print(line_number,line_number[-1])
            num_CC.append(line_number[-1])
            i+=1
            if i==10:
                break

        if search_text in line:
            flag=True

    # print(num_CC,np.array(num_CC,dtype=float))
    num_CC_array=np.array(num_CC,dtype=float)
    num_count=np.isnan(num_CC_array).sum()
    return num_count
file_path=sys.argv[1]
num_count=nan_cal(file_path,"the number of multiplicity point")
print(num_count)

