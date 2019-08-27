import numpy as np
from matplotlib import pyplot as plt

T=0.001
g=1.0

w0,w1,k0,k1,f=np.loadtxt("f.txt",unpack=True)

k=k0
p=k1
w=w0+w1
w0a=4*np.pi*g*(np.log(np.abs((k+p)/(k-p)))-0.5*g/(w**2+g)*(np.log(((k+p)/(k-p))**2)-np.log(np.abs((((k+p)**2+w**2+g)/(((k-p)**2)+w**2+g))))))/k/p
w=np.abs(w0-w1)
w0a+=4*np.pi*g*(np.log(np.abs((k+p)/(k-p)))-0.5*g/(w**2+g)*(np.log(((k+p)/(k-p))**2)-np.log(np.abs((((k+p)**2+w**2+g)/(((k-p)**2)+w**2+g))))))/k/p


fa=w0a*T*p**2/(2*np.pi)**2/(w1**2+(p**2-1)**2)

np.savetxt("f_check.txt",np.column_stack((w0,w1,k0,k1,f,fa)),fmt="%.10e")

plt.figure()
plt.plot(w0,k0,",")
plt.savefig("dist0.png")
plt.figure()
plt.plot(w1,k1,",")
plt.savefig("dist1.png")

err=np.abs(w0-w1)/w1*np.abs(k0-k1)/k1

print [i for i in range(len(err)) if err[i]==err.max()]
print [w1[i] for i in range(len(err)) if err[i]==err.max()]
print [k1[i] for i in range(len(err)) if err[i]==err.max()]
print [err[i] for i in range(len(err)) if err[i]==err.max()]
