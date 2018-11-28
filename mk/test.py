# -*- coding:utf-8 -*-

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as C
import math

def solver(w, t):
	Oq, COq, q= w
	return [
		
-(+0.10493458023251406*q*0.3-0.10493458023251406*COq-0.19393592065823853*COq*Oq+3.0761508903176893e-55*0.1*q*q)
-(-0.10493458023251406*q*0.3+0.10493458023251406*COq-0.021219747830592432*0.6*q*q+9.923890644120988e-05*Oq*Oq-0.021219747830592432*0.6*q*q+9.923890644120988e-05*Oq*Oq+0.19393592065823853*COq*Oq-3.0761508903176893e-55*0.1*q*q+0.19393592065823853*COq*Oq-3.0761508903176893e-55*0.1*q*q)
,		
+0.10493458023251406*q*0.3-0.10493458023251406*COq-0.19393592065823853*COq*Oq+3.0761508903176893e-55*0.1*q*q
,		
-0.10493458023251406*q*0.3+0.10493458023251406*COq-0.021219747830592432*0.6*q*q+9.923890644120988e-05*Oq*Oq-0.021219747830592432*0.6*q*q+9.923890644120988e-05*Oq*Oq+0.19393592065823853*COq*Oq-3.0761508903176893e-55*0.1*q*q+0.19393592065823853*COq*Oq-3.0761508903176893e-55*0.1*q*q
,		]
t = np.logspace(-20., 15., 1000001)
track = odeint(solver, [0.0,0.0,1.0], t)
print(track)
conlist = []
for i in track[-1]:
	conlist.append(i)
Oq, COq, q= conlist


print(+0.10493458023251406*q*0.3-0.10493458023251406*COq)

print(+0.021219747830592432*0.6*q*q-9.923890644120988e-05*Oq*Oq)

print(+0.19393592065823853*COq*Oq-3.0761508903176893e-55*0.1*q*q)
print(conlist)
plt.plot(t, track[:, 0], 'b', label="O*")
plt.plot(t, track[:, 1], 'g', label="CO*")
plt.plot(t, track[:, 2], 'r', label="*")
plt.xscale('symlog')    
plt.legend(loc='best')    
plt.xlabel('t')    
plt.show()