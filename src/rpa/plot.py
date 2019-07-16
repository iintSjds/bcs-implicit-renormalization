import numpy as np
from matplotlib import pyplot as plt

w,p,delta=np.loadtxt("delta.txt",unpack=True)

plt.figure()
plt.plot(p,delta,",")
plt.savefig("delta.png")
