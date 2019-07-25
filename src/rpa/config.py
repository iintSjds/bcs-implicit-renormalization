# this is for generating configurations of calc_pi
# w_exp and calc_delta

import numpy as np
import sys

T=0.001
mu=1.0
m=0.5
#rs=1.0/1.0421235224
rs=2.5

Kc=15
wc=15

qc=0.00000001
kc=0.000001

if len(sys.argv)==2:
    rs=float(sys.argv[1])

# calculate grids
w=[0,] # freq grid of pi and w
count=0
multi=1
while w[-1]<2*wc:
    if w[-1]<0.1*mu:
        w.append(w[-1]+2*np.pi*T*multi)
    elif w[-1]<4*mu:
        w.append(w[-1]+((0.1*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    else:
        w.append(w[-1]+((0.5*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    count+=1
    if count==7:
        count=0
        multi*=2

v=[np.pi*T,] #freq grid of delta
count=0
multi=1
while v[-1]<wc:
    if v[-1]<0.1*mu:
        v.append(v[-1]+2*np.pi*T*multi)
    elif v[-1]<0.2*mu: 
        v.append(v[-1]+((0.05*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    elif v[-1]<0.4*mu:
        v.append(v[-1]+((0.02*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    elif v[-1]<2*mu:
        v.append(v[-1]+((0.1*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    else:
        v.append(v[-1]+((0.8*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    count+=1
    if count==5:
        count=0
        multi*=2

q=[qc,]
while q[-1]<2*Kc:
    if q[-1]<0.1*mu:
        step=q[-1]
        for i in range(3):
            q.append(q[-1]+step)
    elif q[-1]<4*mu:
        q.append(q[-1]+0.2*mu)
    elif q[-1]<8*mu:
        q.append(q[-1]+1.0*mu)
    else:
        q.append(q[-1]+2.0*mu)

k=[kc,]
Nlog=20
while k[-1]<Kc:
    if k[-1]<0.2*mu:
        step=k[-1]
        for i in range(3):
            k.append(k[-1]+step)
    elif k[-1]<0.7*mu:
        k.append(k[-1]+0.1*mu)
    elif (mu-k[-1])>kc:
        n=3
        step=(mu-k[-1])/n
        for j in range(n-1):
            k.append(k[-1]+step)
    elif k[-1]<mu:
        k.append(mu)
        k.append(mu+kc)
    elif k[-1]<1.3*mu:
        step=k[-1]-mu
        for i in range(2):
            k.append(k[-1]+step)
    elif k[-1]<4.0*mu:
        k.append(k[-1]+0.2*mu)
    else:
        k.append(k[-1]+1.0*mu)

# write pi.in, containing configs for calc_pi
file=open("pi.in","w")
file.write("set T=%.16f\n" %T)
file.write("set mu=%.16f\n" %mu)
file.write("set m=%.16f\n" %m)
file.write("set freq=grid\n")
for i in w:
    file.write("%.16f\n" %i)
file.write("end grid\n")
file.write("set mmt=grid\n")
for i in q:
    file.write("%.16f\n" %i)
file.write("end grid\n")
file.close()

# write w.in
file=open("w.in","w")
file.write("set T=%.16f\n" %T)
file.write("set mu=%.16f\n" %mu)
file.write("set m=%.16f\n" %m)
file.write("set rs=%.16f\n" %rs)
file.write("set mmt=grid\n")
for i in k:
    file.write("%.16f\n" %i)
file.write("end grid\n")
file.close()

# write delta.in
file=open("delta.in","w")
file.write("set T=%.16f\n" %T)
file.write("set mu=%.16f\n" %mu)
file.write("set m=%.16f\n" %m)
file.write("set rs=%.16f\n" %rs)
file.write("set v=grid\n")
for i in v:
    file.write("%.16f\n" %i)
file.write("end grid\n")
file.close()
