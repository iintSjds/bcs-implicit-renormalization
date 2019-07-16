import numpy as np

w,q,P=np.loadtxt("pi.txt",unpack=True)

Pa=q**2/4/np.pi/(w**2+q**2)

np.savetxt("pi_check.txt",np.column_stack((w,q,P,Pa)),fmt="%.16f")

w,k,H=np.loadtxt("h1.txt",unpack=True)
g=1.0
Ha=4*np.pi*g*(np.log(k)-g/(w**2+g)*(np.log(k)-0.5*np.log(k**2+w**2+g)))

np.savetxt("h1_check.txt",np.column_stack((w,k,H,Ha,H-Ha)),fmt="%.16f")

w,k,p,w0=np.loadtxt("w0.txt",unpack=True)
k=k+0.0000000001
g=1.0
w0a=4*np.pi*g*(np.log(np.abs((k+p)/(k-p)))-0.5*g/(w**2+g)*(np.log(((k+p)/(k-p))**2)-np.log(np.abs((((k+p)**2+w**2+g)/(((k-p)**2)+w**2+g))))))/k/p
np.savetxt("w0_check.txt",np.column_stack((w,k,p,w0,w0a)),fmt="%.16f")
