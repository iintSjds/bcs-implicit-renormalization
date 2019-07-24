import numpy as np
import h5py
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

f=h5py.File("delta.h5","r")
w=np.array(f["/delta/w"])
q=np.array(f["/delta/q"])
d=np.array(f["/delta/delta"])
f.close()

d.reshape((len(w),len(q)))
ww=np.array([w for i in range(len(q))])
qq=np.array([q for i in range(len(w))])

ww=ww.T

np.savetxt("delta_plot.txt",np.column_stack((ww.ravel(),qq.ravel(),d.ravel())),fmt="%.16f")

fig=plt.figure()
ax=fig.add_subplot(111,projection="3d")

ax.scatter(ww,qq,d,marker=".")
plt.show()
