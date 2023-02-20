import csv

import matplotlib.pyplot as plt
from numpy import *
from numpy import polynomial
import pandas as pd
import Setting_Layers as sl
import Reflectivity as ref
import tools

from scipy import *
import sympy as smp
from scipy.interpolate import interp1d 

dset1 = pd.read_csv("Default Dataset.csv", encoding='latin1')
dset1.dropna()

energy = dset1['E_2']/100000
lambda_i = dset1['l'] + dset1['L_2']/100000

p =  interp1d(lambda_i, energy, kind='cubic')

ACF = (pi/180)  # Angle Conversion Factor.

thickness = [1, 50*1E-9, 1]  # Thickness of each layer
material = [1, 13, 16]   # Material of each layer
R_Tm = []       # Stores lists with reflectivity curves for each of the interactions
ref_index = [complex(1.4826,0), complex(0,0), complex(1.33,0)]  # Refractive index of the layers
n_layers=3

lambda_i = arange(750, 935, 1) * 1e-9 

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

plt.subplots(dpi = 200)
plt.plot(theta_i/ACF, R_Tm[0], label ='Com espalhamento')
plt.xlabel("Ângulo de incidância (°)")
plt.ylabel("Reflectância normalizada")
tools.data_processing(ref_index)




"""dados = pd.read_csv("experimento_H2O.csv", encoding='latin1')
pixel = dados['pixel']
#pixel = pixel[::-1]
signal = dados[' signal']
theta = []

for i in pixel:
    theta.append((3.1522 * 1E-5 * i**2) - (0.0661 * i) + 73.4533)

plt.plot(theta, signal, '--')
plt.show()
"""
"""
dset2 = pd.read_csv("Default Dataset 2.csv", encoding='latin1')
dset1 = pd.read_csv("Default Dataset.csv", encoding='latin1')

dset1.dropna()
dset2.dropna()

energy = dset1['e2']/10000
lambda_i = dset1['l'] + dset1['l2']/10000


p_20 = polynomial.Polynomial.fit(lambda_i, energy, 20)


print([p_20])

y=[]
yp = []
y2 = []


x = list(linspace(740, 940, 300))

for i in x:
    a = p(i)
    yp.append(p_20(i))
    y2.append(p(i))
    y.append(a)



fig, ax = plt.subplots(dpi=200)
plt.plot(lambda_i, energy, label='base')
#plt.plot(x, y, label='Aproximação')

plt.plot(x, yp, label='Polinomio grau 16')

plt.plot(x, y2, label='aprox2')
plt.xlabel(f"Comprimento de onda - {chr(955)} - (nm)")
plt.ylabel(f"Intensidade normalizada")
plt.legend()
plt.grid()
plt.show()
"""