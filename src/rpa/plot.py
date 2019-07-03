import numpy as np
from matplotlib import pyplot as plt

w,p,P,W=np.loadtxt("result.txt",unpack=True)

plt.figure()
plt.plot(p,P)
plt.savefig("Pp.png")
