import numpy as np
import gemmi
#import matplotlib.pyplot as plt
import sys

mtz_file=sys.argv[1]

mtz=gemmi.read_mtz_file(mtz_file)

mtz_array=np.array(mtz,copy=False)
#print(mtz.column_labels(),mtz_array.shape)


# EA_array=mtz_array[:,7]
# SIGEA_array=mtz_array[:,8]
EA_array=mtz.column_with_label('EA')
SIGEA_array=mtz.column_with_label('SIGEA')
EA_array=np.array(EA_array,copy=True)
SIGEA_array=np.array(SIGEA_array,copy=True)

EA_array=EA_array[~np.isnan(EA_array)]
SIGEA_array=SIGEA_array[~np.isnan(SIGEA_array)]
#print(EA_array.shape,SIGEA_array.shape)


d_thres=2.1
sigma_hist=np.zeros([20001])
#print(sigma_hist.shape)

total=0
for i in range(EA_array.shape[0]):
    if SIGEA_array[i]>0.0001 and EA_array[i]<d_thres and EA_array[i]>0.0001:
        sigma_hist[ min( int(SIGEA_array[i]*10000),20000 ) ] +=1
        total+=1

cumul=0
cumul_thresh=max(500, int(total*0.35))
for i in range(20001):
    cumul+=sigma_hist[i]
    if cumul>=cumul_thresh:
        sigd_thres=(i+1)/10000
        break
#print(sigd_thres,cumul_thresh)
print(sigd_thres)
