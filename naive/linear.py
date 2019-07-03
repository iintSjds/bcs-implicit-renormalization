import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import sys

T,lam=np.loadtxt(sys.argv[1],unpack=True)

plt.plot(-np.log(T),lam)
plt.savefig(sys.argv[1]+".png")

print(stats.linregress(-np.log(T),lam))

slope,intercept,a,b,c=stats.linregress(-np.log(T),lam)
print np.exp(-(1.0-intercept)/slope)
