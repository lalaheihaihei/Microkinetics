# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: ode
@time: 2018/11/27 17:35
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as C
import math

def solver(w, t):
    q, Oq = w
    return [
        - 1.0388318052209737e-06 * q * 0.3 + 1.218567544319086e-06 * (
        1 - Oq - q) - 0.017150364238188037 * 0.6 * q * q + 7.758133527434427e-05 * Oq * Oq + 0.18640526909338634 * (
        1 - Oq - q) * Oq - 1.1686382736183174e-55 * 0.1 * q * q
        ,
        + 0.017150364238188037 * 0.6 * q * q - 7.758133527434427e-05 * Oq * Oq - 0.18640526909338634 * (
        1 - Oq - q) * Oq + 1.1686382736183174e-55 * 0.1 * q * q

    ]

t = np.logspace(-20., 8., 1000001) # 创建时间点

track = odeint(solver, [1., 0.], t)
print(track)

'''
q, Oq = track[-1][0], track[-1][1]
dqdt = - 1.0388318052209737e-06 * q * 0.3 + 1.218567544319086e-06 * (
        1 - Oq - q) - 0.017150364238188037 * 0.6 * q * q + 7.758133527434427e-05 * Oq * Oq + 0.18640526909338634 * (
        1 - Oq - q) * Oq - 1.1686382736183174e-55 * 0.1 * q * q
dOqdt = + 0.017150364238188037 * 0.6 * q * q - 7.758133527434427e-05 * Oq * Oq - 0.18640526909338634 * (
        1 - Oq - q) * Oq + 1.1686382736183174e-55 * 0.1 * q * q
print(dqdt, dOqdt)
'''

plt.plot(t, track[:, 0], 'b', label='q')
plt.plot(t, track[:, 1], 'g', label='Oq')
plt.xscale('symlog')
plt.legend(loc='best')
plt.xlabel('t')
plt.show()


