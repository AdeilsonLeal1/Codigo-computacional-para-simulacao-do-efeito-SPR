# In the "tools" module, tools are presented to calculate the SPR resonance point and plot the reflectivity curve

from cmath import pi
import csv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import butter, lfilter, freqz, filtfilt


def point_SPR(reflet, ax_x, modo):
    # The "Point_SPR" method returns the resonance point of the curve,
    # either the resonance angle in degrees or the resonance wavelength in nanometers

    # Receives the position of the minimum point of the curve
    c = reflet.index(min(reflet))
    if modo == 1:
        return ax_x[c] * (180 / pi)  # Returns the angle in degrees
    else:
        return ax_x[c] * 1E9  # Returns the wavelength in nanometers


def plot(x_i, R_Tm, resonance_point, modo):

    # Correspondence between AIM and WIM modes
    y = 'Angle' if modo == 1 else 'Wavelength'
    c = '°' if modo == 1 else 'nm'
    n = chr(952) if modo == 1 else chr(955)
    ax_x = x_i * (180 / pi) if modo == 1 else x_i * 1E9

    fig0, ax_TM = plt.subplots(dpi=200)
    ax_TM.plot(ax_x, (R_Tm[0]))
    ax_TM.set_title(f"Reflectivity vs. {y} of Incidence",
                    loc='center', pad='6')
    ax_TM.set(xlabel=f'Incidence {y} ({c})', ylabel='Reflectivity')
    text = f"{n}$_S$$_P$$_R$ = {resonance_point:.4f}{c}"
    ax_TM.annotate(text, (resonance_point, min(R_Tm[0])), xytext=(ax_x[0], 0.1), arrowprops=dict(facecolor='blue', arrowstyle="->"),
                   bbox={'facecolor': 'white', 'edgecolor': 'gray', 'alpha': 0.7})
    ax_TM.grid()
    ax_TM.set_yticks(np.arange(0, 1.20, 0.20))

    plt.show()


def data_processing(ref_index):
    refletivity = pd.read_csv("reflectivity.csv")
    refletivity_processing = pd.read_csv("reflectivity_noise.csv")
    refletivity_processing['MMS'] = refletivity_processing['Reflectivity'].rolling(
        window=10).mean()
    refletivity_processing['MME'] = refletivity_processing['Reflectivity'].ewm(
        span=2).mean()
    refletivity_processing['MME2'] = refletivity_processing['MME'].ewm(
        span=2).mean()
    refletivity_processing['MME4'] = refletivity_processing['MME2'].ewm(
        span=2).mean()
    refletivity_processing['MME8'] = refletivity_processing['MME4'].ewm(
        span=2).mean()
    refletivity_processing['MME16'] = refletivity_processing['MME8'].ewm(
        span=2).mean()
    
    resonance_point = point_SPR(list(refletivity_processing['MME8']),list(refletivity['Angle']*pi/180), 1)
    
    print(resonance_point)
    print(f"Indice indice após processamento{return_analyte(ref_index, resonance_point)}")
    
    signal_filtered = butter_lowpass_filter(data=refletivity_processing['Reflectivity'], cutoff=2, fs=30, order=5 )

    fig1, fig = plt.subplots(dpi=200)
    fig.plot(refletivity['Angle'],
             refletivity['Reflectivity'], label="Theorical")
    fig.plot(refletivity_processing['Angle'],
             refletivity_processing['Reflectivity'], label="Sinal bruto")
    fig.plot(refletivity_processing['Angle'],
             signal_filtered, label="Signal filtered")
    fig.plot(refletivity_processing['Angle'],
             refletivity_processing['MME8'], label="MME8")

    fig.legend()
    fig.grid()
    plt.savefig('Curva SPR.png')
    plt.show()
    plt.close()


def save_csv(name, x, y):
    # 1. create the file
    f = open(name, 'w', newline='', encoding='utf-8')

    # 2. create the recording object
    w = csv.writer(f)

    # 3. record the lines
    w.writerow(['Angle', 'Reflectivity'])

    for i in range(len(y)):
        w.writerow([x[i], y[i]])

    # close the files
    f.close()


def return_analyte(refractive_index, point_SPR):
    print(refractive_index)
    emr = np.real(refractive_index[1])**2 - np.imag(refractive_index[1])**2
    index_analyte = np.sqrt((emr*(refractive_index[0]*np.sin(point_SPR*pi/180))**2)/
        (emr-(refractive_index[0]*np.sin(point_SPR*pi/180))**2))
    
    return np.round(index_analyte,3)


def butter_lowpass_filter(data, cutoff, fs, order):
    normal_cutoff = cutoff / (0.5 * fs)
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y