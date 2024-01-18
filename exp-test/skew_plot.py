import  numpy as np
import matplotlib.pyplot as plt
import sys

filename=sys.argv[1]

with open(filename) as file:
    skew_str=file.readlines()

skew_array=np.array(skew_str).astype(float)

plt.plot(skew_array)
plt.show()
