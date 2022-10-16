# In the "tools" module, tools are presented to calculate the SPR resonance point and the sensor sensitivity

from cmath import pi


def Point_SPR(refletancia, axe_x, modo):
        # The "Point_SPR" method returns the resonance point of the curve, 
        # either the resonance angle in degrees or the resonance wavelength in nanometers

        c = refletancia.index(min(refletancia))  # Receives the position of the minimum point of the curve
        if modo == 1:
            return axe_x[c] * (180 / pi)  # Returns the angle in degrees
        else:
            return axe_x[c] * 1E9  # Returns the wavelength in nanometers