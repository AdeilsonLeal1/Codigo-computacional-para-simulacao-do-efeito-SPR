##########     Reflectivity     ##########
# module that calculates the Fresnel coefficient

from cmath import cos, pi, sin, sqrt
from numpy import array


def reflectivity(nLayers, d,  index, theta_i, wavelenght):

    j = complex(0, 1)  # Simplification for the complex number "j"
    k0 = (2 * pi) / wavelenght  # Wave number

    b = []  # b_j -> Phase shift in each layer
    q_TM = []  # q_TM_j -> Admittance in TM polarization
    M_Tm = []  # M_Tm_j -> Transfer matrix between each layer - TM polarization

    for layer in range(nLayers):
        y = sqrt((index[layer] ** 2) - ((index[0] * sin(theta_i)) ** 2))

        b.append(k0 * d[layer] * y)
        q_TM.append(y / index[layer] ** 2)

        # Total Transfer Matrix
        if layer < (nLayers - 1):
            M_Tm.append(array([[cos(b[layer]), (-j / q_TM[layer]) * sin(b[layer])],
                               [-j * q_TM[layer] * sin(b[layer]), cos(b[layer])]]))

    Mt_TM = M_Tm[0]  # Mt_TM -> Total Transfer Matrix - TM polarization

    for k in range(nLayers - 2):
        Mt_TM = Mt_TM @ M_Tm[k + 1]

    num_TM = (Mt_TM[0][0] + Mt_TM[0][1] * q_TM[nLayers - 1]) * q_TM[0] - (
        Mt_TM[1][0] + Mt_TM[1][1] * q_TM[nLayers - 1])
    den_TM = (Mt_TM[0][0] + Mt_TM[0][1] * q_TM[nLayers - 1]) * q_TM[0] + (
        Mt_TM[1][0] + Mt_TM[1][1] * q_TM[nLayers - 1])

    r_TM = num_TM / den_TM  # 'r_TM'-> Fresnel reflection coefficient - TM polarization

    return (abs(r_TM) ** 2)  # Reflectance - TM polarization
