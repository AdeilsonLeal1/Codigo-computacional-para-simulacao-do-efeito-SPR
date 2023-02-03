import csv

import matplotlib.pyplot as plt
from numpy import *
import pandas as pd
import Setting_Layers as sl
import Reflectivity as ref
import tools

from scipy import *
import sympy as smp


def p(i):
    dc = 0.79
    a = 1/(sqrt(2*pi))

    frac1 = (i - 845)/30
    frac2 = (i - 853)/15
    frac3 = (i - 845)/42

    f = dc + (a/2.7)* exp(-(frac1**2))
    g = dc + (a/2.7)* exp(-(frac2**2))
    h = 1.01 + (a/3)* exp(-(frac3**2))

    a = f*g*h
    return a

def p_2(i):
    dc = 0.632
    a = 1/(sqrt(2*pi))
    b = 30

    
    f = dc + (a/1.1) * exp(-((i - 850)**2)/(b**2))

    a = f
    return a


"""ACF = (pi/180)  # Angle Conversion Factor.

thickness = [1, 50*1E-9, 1]  # Thickness of each layer
material = [1, 13, 16]   # Material of each layer
R_Tm = []       # Stores lists with reflectivity curves for each of the interactions
ref_index = [complex(1.4826,0), complex(0,0), complex(1.33,0)]  # Refractive index of the layers
n_layers=3

lambda_i = arange(400, 1000, 1) * 1e-9 

# Starting and ending angle
a1 = 65.5090 * ACF 
a2 = 73.3877 * ACF

theta_i = linspace(a1, a2, 128)  # Array with angles of incidence

#  Defining the structure

r_TM = []      # Local variable that temporarily stores reflectivity values for each of the interactions

for theta in theta_i:
    ref_index = [complex(1.4826,0), complex(0,0), complex(1.33,0)]
    Num = []
    Den = []
    s1=0
    s2=0
    for lamb in lambda_i:
        ref_index[1] = sl.set_index(material[1], lamb)
        a, b = ref.reflectivity(n_layers, thickness, ref_index, theta, lamb),  p(lamb*1E9)
        Num.append(a * b)
        s1 = s1 + Num[-1]
        Den.append(p(lamb*1E9))
        s2= s2+Den[-1]

    r_TM.append(s1/s2)


R_Tm.append(r_TM)

plt.plot(theta_i/ACF, R_Tm[0])
tools.data_processing(ref_index)"""






"""dados = pd.read_csv("experimento_H2O.csv", encoding='latin1')
pixel = dados['pixel']
#pixel = pixel[::-1]
signal = dados[' signal']
theta = []

for i in pixel:
    theta.append((3.1522 * 1E-5 * i**2) - (0.0661 * i) + 73.4533)

plt.plot(theta, signal, 'o')
plt.show()

"""
"""dset2 = pd.read_csv("Default Dataset 2.csv", encoding='latin1')
dset1 = pd.read_csv("Default Dataset.csv", encoding='latin1')

dset1.dropna()
dset2.dropna()

energy = dset1['e2']/10000
lambda_i = dset1['l'] + dset1['l2']/10000

energy2 = dset2['e2']/10000
lambda_i2 = dset2['l'] + dset2['l2']/10000

y=[]
y2 = []
x = list(linspace(400, 1000, 300))

for i in x:
    a = p(i)
    y2.append(p_2(i))
    y.append(a)

plt.plot(lambda_i, energy, label='base')
plt.plot(x, y, label='aprox')
plt.plot(x, y2, label='aprox2')
plt.legend()
plt.grid()
plt.show()

"""

Idc = float(1.33000)

tau = 200

idcs = [0.0010, 0, 0.0005, 0, 0.0005, 0.0015, 0.004, 0]

carga = [True, False, True, False, True, True, False, False]
count = 0
ri = []
while (count < 1000):
    ri.append(1.33)
    count+=1
count = 0

for dc, load in zip(idcs, carga):
    if load:
        while (count < 1500):
                v = Idc + dc*(1 - exp(-count/tau))
                ri.append(v)
                count+=1
    else:
        while (count < 2500):
                v = (Idc+dc) - dc*(exp(-count/tau))
                ri.append(v)
                count+=1
    
    count = 0

x = linspace (0,len(ri), len(ri) )
y = [1.33 for i in x]
plt.plot(x, ri, '--')
plt.plot(x, y)
plt.yticks(arange(1.329, 1.334, 0.001))
plt.show()