import numpy as np
import h5py

w,q,P=np.loadtxt("pi.txt",unpack=True)
w1,q1,P1=np.loadtxt("pi1.txt",unpack=True)

#Pa=q**2/4/np.pi/(w**2+q**2)


Pa=P1

np.savetxt("pi_check.txt",np.column_stack((w,q,P,Pa)),fmt="%.16f")

# w,k,H=np.loadtxt("h1.txt",unpack=True)
# g=1.0
# Ha=4*np.pi*g*(np.log(k)-g/(w**2+g)*(np.log(k)-0.5*np.log(k**2+w**2+g)))

# np.savetxt("h1_check.txt",np.column_stack((w,k,H,Ha,H-Ha)),fmt="%.16f")

# w,k,p,w0=np.loadtxt("w0.txt",unpack=True)
# k=k+1e-15
# g=1.0
# w0a=4*np.pi*g*(np.log(np.abs((k+p)/(k-p)))-0.5*g/(w**2+g)*(np.log(((k+p)/(k-p))**2)-np.log(np.abs((((k+p)**2+w**2+g)/(((k-p)**2)+w**2+g))))))/k/p
# np.savetxt("w0_check.txt",np.column_stack((w,k,p,w0,w0a)),fmt="%.16f")

# f=h5py.File("w.h5","r+")
# f["/w0/w0"].write_direct(w0a)
# f.close()
