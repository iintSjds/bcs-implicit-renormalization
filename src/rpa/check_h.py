import numpy as np
import h5py
from matplotlib import pyplot as plt

g=3.0*1.0421235224
# f=h5py.File("pi_0.h5","r")
# w=np.array(f["/pi/w"])
# q=np.array(f["/pi/q"])
# h=np.array(f["/pi/pi"])
# f.close()

w,q,P=np.loadtxt("pi.txt",unpack=True)

# ww=np.array([w for i in range(len(q))])
# qq=np.array([q for i in range(len(w))])
# ww=ww.T

ww=w
qq=q

h=P

# ww=ww.ravel()
# qq=qq.ravel()

X=len([i for i in qq if i==qq[0]])
Y=len(qq)/X

print X,Y

#ha=4*np.pi*g*(np.log(qq)-g/(ww**2+g)*(np.log(qq)-0.5*np.log(qq**2+ww**2+g)))

#aha=ha-np.array([ha[(i//len(q))*len(q)] for i in range(len(ha))])

np.savetxt("w.txt",np.column_stack((ww,qq,h)),fmt="%.16f")

W=4*np.pi*g/(qq**2+4*np.pi*g*h)

N=0
np.savetxt("w0.txt",np.column_stack((ww[N*Y:(N+1)*Y],qq[N*Y:(N+1)*Y],W[N*Y:(N+1)*Y])),fmt="%.16f")
N=1
np.savetxt("w1.txt",np.column_stack((ww[N*Y:(N+1)*Y],qq[N*Y:(N+1)*Y],W[N*Y:(N+1)*Y],(qq*qq*W)[N*Y:(N+1)*Y])),fmt="%.16f")
N=2
np.savetxt("w2.txt",np.column_stack((ww[N*Y:(N+1)*Y],qq[N*Y:(N+1)*Y],W[N*Y:(N+1)*Y],(qq*qq*W)[N*Y:(N+1)*Y])),fmt="%.16f")

N=0
M=60
q0=qq[N*Y:(N+1)*Y]
plt.plot(qq[N*Y:(N+1)*Y][:M],W[N*Y:(N+1)*Y][:M])
plt.plot(q0[:M],4*np.pi*g/(q0**2+2/np.pi*g)[:M])
plt.savefig("w0.png")

N=2
i=1
M=30

ratio=4*np.pi*g*2./3./np.pi**2

plt.figure()
plt.plot(qq[N*Y:(N+i)*Y][:M],(W*qq*qq)[N*Y:(N+i)*Y][:M])
plt.plot(q0[:M],4*np.pi*g/(1+ratio/(ww[N*Y:(N+i)*Y][0])**2)*np.ones(M))
plt.savefig("w1.png")

print (4*np.pi*g*2./3./np.pi**2),1.629032**2
