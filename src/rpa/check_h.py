import numpy as np
import h5py


g=1.0

f=h5py.File("h.h5","r")
w=np.array(f["/h/w"])
q=np.array(f["/h/q"])
h=np.array(f["/h/h1"])
f.close()

ww=np.array([w for i in range(len(q))])
qq=np.array([q for i in range(len(w))])
ww=ww.T

ww=ww.ravel()
qq=qq.ravel()

ha=4*np.pi*g*(np.log(qq)-g/(ww**2+g)*(np.log(qq)-0.5*np.log(qq**2+ww**2+g)))

np.savetxt("h_check.txt",np.column_stack((ww,qq,h,ha,(h-ha)/(h+1e-15))),fmt="%.16f")

