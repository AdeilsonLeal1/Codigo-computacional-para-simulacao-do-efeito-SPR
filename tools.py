# In the "tools" module, tools are presented to calculate the SPR resonance point and plot the reflectivity curve

from cmath import pi

import matplotlib.pyplot as plt
import numpy as np


def Point_SPR(refletancia, ax_x, modo):
    # The "Point_SPR" method returns the resonance point of the curve,
    # either the resonance angle in degrees or the resonance wavelength in nanometers

    # Receives the position of the minimum point of the curve
    c = refletancia.index(min(refletancia))
    if modo == 1:
        return ax_x[c] * (180 / pi)  # Returns the angle in degrees
    else:
        return ax_x[c] * 1E9  # Returns the wavelength in nanometers

def plot(x_i, R_Tm, resonance_point, a1, modo):
    
    # Correspondence between AIM and WIM modes
    y = 'Angle' if modo == 1 else 'Wavelength'
    c = '°' if modo == 1 else 'nm'
    n = chr(952) if modo == 1 else chr(955)
    ax_x = x_i * (180 / pi) if modo == 1 else x_i * 1E9
    z = (180 / pi) if modo == 1 else 1E9

    fig0, ax_TM = plt.subplots(dpi=200)
    ax_TM.plot(ax_x, (R_Tm[0]))
    ax_TM.set_title(f"Reflectivity vs. {y} of Incidence",
                        loc='center', pad='6')
    ax_TM.set(xlabel=f'Incidence {y} ({c})', ylabel='Reflectivity')
    text = f"{n}$_S$$_P$$_R$ = {resonance_point:.4f}{c}"
    ax_TM.annotate(text, (resonance_point, min(R_Tm[0])), xytext=((a1 * z), 0.1), arrowprops=dict(facecolor='blue', arrowstyle="->"),
                    bbox={'facecolor': 'white', 'edgecolor': 'gray', 'alpha': 0.7})
    ax_TM.grid()
    ax_TM.set_yticks(np.arange(0, 1.20, 0.20))

    plt.show()
