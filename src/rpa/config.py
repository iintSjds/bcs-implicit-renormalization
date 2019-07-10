# this is for generating configurations of calc_pi
# w_exp and calc_delta

import numpy as np
import sys

T=0.001
mu=1.0
m=0.5
rs=1.0

kc=20.0
wc=20.0

qc=0.0001

if len(sys.argv)==2:
    rs=float(sys.argv[1])

# calculate grids
w=[0,]
count=0
multi=1
while w[-1]<wc:
    if w[-1]<0.1*mu:
        w.append(w[-1]+2*np.pi*T*multi)
    elif w[-1]<2*mu:
        w.append(w[-1]+0.1*mu)
    else:
        w.append(w[-1]+0.5*mu)
    count+=1
    if count==7:
        count=0
        multi*=2

q=[qc,]
while q[-1]<kc:
    if q[-1]<0.1*mu:
        q.append(q[-1]*2)
    elif q[-1]<4*mu:
        q.append(q[-1]+0.1*mu)
    else:
        q.append(q[-1]+0.5*mu)

# write pi.in, containing configs for calc_pi
file=open("pi.in","w")
file.write("set T=%f\n" %T)
file.write("set mu=%f\n" %mu)
file.write("set m=%f\n" %m)
file.write("set grid=freq\n")
for i in w:
    file.write("%f\n" %i)
file.write("end grid\n")
file.write("set grid=mmt\n")
for i in q:
    file.write("%f\n" %i)
file.write("end grid\n")
file.close()
