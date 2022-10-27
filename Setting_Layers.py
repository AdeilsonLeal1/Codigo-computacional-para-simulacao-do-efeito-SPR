##########     Setting_Layers     ##########
# module that includes the functions needed to build the layers of the structure
# setLayers -> sets the material for each layer
from cmath import sqrt
from numpy import interp, real, imag


def setLayers(layer):

    if layer == 0:
        print("\n=====================================================\n"
                "==========    Set Layers Characteristics   ==========\n"
                "=====================================================\n\n"
                "==========  1st Layer - Optical substrate  ==========\n")
        d = 1
        material = int(input(f"\n1 - BK7   2 - Silica   3 - N-F2   4 - Synthetic sapphire(Al2O3)"
                                    f"\n5 - SF10  6 - FK51A    7 - N-SF14 8 - Acrylic SUVT "
                                    f"\n9 - Other - (Custom only in AIM mode)  "
                                    f"\n\nMaterial -> "))
    else:
        print(
            f"\n==========            {layer+1}st Layer            ==========\n")
        material= int(input(f"\n 1 - BK7     2 - Silica      3 - N-F2       4 - Synthetic sapphire(Al2O3) "
                                    f"\n 5 - SF10    6 - FK51A       7 - N-SF14     8 - Acrylic SUVT"
                                    f"\n 9 - PVA    10 - Glycerin   11 - Quartz    12 - Aluminium"
                                    f"\n13 - Gold   14 - Silver     15 - Copper    16 - Water (RI = 1.33)"
                                    f"\n17 - Air    18 - LiF        19 - Cytop     20 - Other - (Custom only in AIM mode)\n\nMaterial -> "))
        
        d = float(input("Thickness (nm): ")) * 1e-9

    return d, material


def setAnalyte(n, analyte):
    analyte_list = []   # Array with the variations in the analyte
    print("\n=====================================================\n"
          "==========  Define Analyte Characteristics  =========\n"
          "=====================================================\n")

    d_na = float(
        input("\nRefractive index variation step (eg. delta_na = 0.005 RIU): "))

    for i in range(n):
        analyte_list.append((analyte + i*d_na))

    return analyte_list


def set_index(material, wi):
    B1, B2, B3, C1, C2, C3, X, n, k_index = 0, 0, 0, 0, 0, 0, [
    ], [], []  # Initialization of variables
    Lambda_i = wi * 1e6  # Incidence Wavelength in micrometers
    j = 0 + 1j
    # Materials that compose the prisms, modeled by the Sellmeier equation according to the RefractiveIndex.info
    # https://refractiveindex.info/
    if 0 < material <= 8:
        if material == 1:  # BK7
            B1, B2, B3, C1, C2, C3 = 1.03961212, 2.31792344E-1, 1.01046945, 6.00069867E-3, 2.00179144E-2, 103.560653

        elif material == 2:  # Silica
            B1, B2, B3, C1, C2, C3 = 0.6961663, 0.4079426, 0.8974794, 4.6791E-3, 1.35121E-2, 97.934003

        elif material == 3:  # N-F2
            B1, B2, B3, C1, C2, C3 = 1.39757037, 1.59201403E-1, 1.2686543, 9.95906143E-3, 5.46931752E-2, 119.2483460

        elif material == 4:  # Synthetic Sapphire(Al2O3)
            B1, B2, B3, C1, C2, C3 = 1.4313493, 0.65054713, 5.3414021, 0.00527993, 0.0142383, 325.01783

        elif material == 5:  # SF10
            B1, B2, B3, C1, C2, C3 = 1.62153902, 0.256287842, 1.64447552, 0.0122241457, 0.0595736775, 147.468793

        elif material == 6:  # FK51A
            B1, B2, B3, C1, C2, C3 = 0.971247817, 0.216901417, 0.904651666, 0.00472301995, 0.0153575612, 168.68133

        elif material == 7:  # N-SF14
            B1, B2, B3, C1, C2, C3 = 1.69022361, 0.288870052, 1.704518700, 0.01305121130, 0.0613691880, 149.5176890

        elif material == 8:  # Acrylic SUVT
            B1, B2, B3, C1, C2, C3 = 0.59411, 0.59423, 0, 0.010837, 0.0099968, 0

        # Sellmeier equation
        n = sqrt(1 + ((B1 * Lambda_i ** 2) / (Lambda_i ** 2 - C1)) + ((B2 * Lambda_i ** 2) / (Lambda_i ** 2 - C2))
                 + ((B3 * Lambda_i ** 2) / (Lambda_i ** 2 - C3)))

    # Glycerol e PVA according to the RefractiveIndex.info: https://refractiveindex.info/
    elif 8 < material <= 10:
        if material == 9:  # PVA
            B1, B2, B3, C1, C2, C3 = 1.460, 0.00665, 0, 0, 0, 0
        elif material == 10:  # Glycerin/glycerol
            B1, B2, B3, C1, C2, C3 = 1.45797, 0.00598, -0.00036, 0, 0, 0
        # Equation that models the refractive index as a function of wavelength
        n = B1 + (B2 / Lambda_i ** 2) + (B3 / Lambda_i ** 4)

    # Quartz according to the RefractiveIndex.info: https://refractiveindex.info/
    elif 10 < material <= 11:
        B1, B2, B3, C1, C2, C3 = 2.356764950, - \
            1.139969240E-2, 1.087416560E-2, 3.320669140E-5, 1.086093460E-5, 0
        n = sqrt(B1 + (B2 * Lambda_i ** 2) + (B3 / Lambda_i ** 2) +
                 (C1 / Lambda_i ** 4) + (C2 / Lambda_i ** 6))

    # Metals
    elif 11 < material <= 15:
        """
        X - Wavelength in micrometers,
        n - Real part of the refractive index e k_index - Imaginary part of the refractive index
        according to Johnson and Christy, 1972
        """
        if material == 12:  # Refractive index of Aluminum according to the Drude model
            LambdaP, LambdaC = 1.0657E-7, 2.4511E-5
            n = sqrt(1 - (((wi ** 2) * LambdaC) /
                     ((LambdaC + (j * wi)) * (LambdaP ** 2))))
        elif material == 13:  # Gold
            X = [0.1879, 0.1916, 0.1953, 0.1993, 0.2033, 0.2073, 0.2119, 0.2164, 0.2214, 0.2262, 0.2313, 0.2371,
                 0.2426, 0.2490, 0.2551, 0.2616, 0.2689, 0.2761, 0.2844, 0.2924, 0.3009, 0.3107, 0.3204, 0.3315,
                 0.3425, 0.3542, 0.3679, 0.3815, 0.3974, 0.4133, 0.4305, 0.4509, 0.4714, 0.4959, 0.5209, 0.5486,
                 0.5821, 0.6168, 0.6595, 0.7045, 0.756, 0.8211, 0.892, 0.984, 1.088, 1.216, 1.393, 1.61, 1.937, 3.5]

            n = [1.28, 1.32, 1.34, 1.33, 1.33, 1.30, 1.30, 1.30, 1.30, 1.31, 1.30, 1.32, 1.32,
                 1.33, 1.33, 1.35, 1.38, 1.43, 1.47, 1.49, 1.53, 1.53, 1.54,
                 1.48, 1.48, 1.50, 1.48, 1.46, 1.47, 1.46, 1.45, 1.38, 1.31, 1.04,
                 0.62, 0.43, 0.29, 0.21, 0.14, 0.13, 0.14, 0.16, 0.17, 0.22, 0.27, 0.35, 0.43, 0.56, 0.92, 1.8]

            k_index = [1.188, 1.203, 1.226, 1.251, 1.277, 1.304, 1.35, 1.387, 1.427, 1.46, 1.497, 1.536, 1.577,
                       1.631, 1.688, 1.749, 1.803, 1.847, 1.869, 1.878, 1.889, 1.893, 1.898, 1.883, 1.871, 1.866, 1.895,
                       1.933, 1.952, 1.958, 1.948, 1.914, 1.849, 1.833, 2.081, 2.455, 2.863, 3.272, 3.697, 4.103, 4.542,
                       5.083, 5.663, 6.35, 7.15, 8.145, 9.519, 11.21, 13.78, 25]

        elif material == 14:  # Silver
            X = [0.1879, 0.1916, 0.1953, 0.1993, 0.2033, 0.2073, 0.2119, 0.2164, 0.2214, 0.2262, 0.2313,
                 0.2371, 0.2426, 0.249, 0.2551, 0.2616, 0.2689, 0.2761, 0.2844, 0.2924, 0.3009, 0.3107,
                 0.3204, 0.3315, 0.3425, 0.3542, 0.3679, 0.3815, 0.3974, 0.4133, 0.4305, 0.4509, 0.4714,
                 0.4959, 0.5209, 0.5486, 0.5821, 0.6168, 0.6595, 0.7045, 0.756, 0.8211, 0.892, 0.984, 1.088,
                 1.216, 1.393, 1.61, 1.937, 5]

            n = [1.07, 1.1, 1.12, 1.14, 1.15, 1.18, 1.2, 1.22, 1.25, 1.26, 1.28, 1.28, 1.3, 1.31, 1.33, 1.35,
                 1.38, 1.41, 1.41, 1.39, 1.34, 1.13, 0.81, 0.17, 0.14, 0.1, 0.07, 0.05, 0.05, 0.05, 0.04, 0.04,
                 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04, 0.03, 0.04, 0.04, 0.04, 0.04, 0.09, 0.13, 0.15,
                 0.24, 2]

            k_index = [1.212, 1.232, 1.255, 1.277, 1.296, 1.312, 1.325, 1.336, 1.342, 1.344, 1.357, 1.367, 1.378,
                       1.389, 1.393, 1.387, 1.372, 1.331, 1.264, 1.161, 0.964, 0.616, 0.392, 0.829, 1.142, 1.419,
                       1.657, 1.864, 2.07, 2.275, 2.462, 2.657, 2.869, 3.093, 3.324, 3.586, 3.858, 4.152, 4.483,
                       4.838, 5.242, 5.727, 6.312, 6.992, 7.795, 8.828, 10.1, 11.85, 14.08, 35]

        elif material == 15:  # Copper
            X = [0.1879, 0.1916, 0.1953, 0.1993, 0.2033, 0.2073, 0.2119, 0.2164, 0.2214, 0.2262, 0.2313, 0.2371,
                 0.2426, 0.249, 0.2551, 0.2616, 0.2689, 0.2761, 0.2844, 0.2924, 0.3009, 0.3107, 0.3204, 0.3315,
                 0.3425, 0.3542, 0.3679, 0.3815, 0.3974, 0.4133, 0.4305, 0.4509, 0.4714, 0.4959, 0.5209, 0.5486,
                 0.5821, 0.6168, 0.6595, 0.7045, 0.756, 0.8211, 0.892, 0.984, 1.088, 1.216, 1.393, 1.61, 1.937, 5]

            n = [0.94, 0.95, 0.97, 0.98, 0.99, 1.01, 1.04, 1.08, 1.13, 1.18, 1.23, 1.28, 1.34, 1.37, 1.41, 1.41,
                 1.45, 1.46, 1.45, 1.42, 1.4, 1.38, 1.38, 1.34, 1.36, 1.37, 1.36, 1.33, 1.32, 1.28, 1.25, 1.24, 1.25,
                 1.22, 1.18, 1.02, 0.7, 0.3, 0.22, 0.21, 0.24, 0.26, 0.3, 0.32, 0.36, 0.48, 0.6, 0.76, 1.09, 2.5]

            k_index = [1.337, 1.388, 1.44, 1.493, 1.55, 1.599, 1.651, 1.699, 1.737, 1.768, 1.792, 1.802, 1.799,
                       1.783, 1.741, 1.691, 1.668, 1.646, 1.633, 1.633, 1.679, 1.729, 1.783, 1.821, 1.864, 1.916, 1.975, 2.045, 2.116,
                       2.207, 2.305, 2.397, 2.483, 2.564, 2.608, 2.577, 2.704, 3.205, 3.747, 4.205, 4.665, 5.18, 5.768, 6.421, 7.217,
                       8.245, 9.439, 11.12, 13.43, 35]

        # Method that calculates the refractive indices of the metals by linear interpolation
        # using the points contained in X, n e k_index previously described
        n_interp = interp(Lambda_i, X, n)
        k_interp = interp(Lambda_i, X, k_index)
        n = complex(n_interp, k_interp)

    # Refractive index of the Water
    elif material == 16:
        n = 1.33

    # Refractive index of the Air
    elif material == 17:
        n = 1.0000

    # Refractive index of the LiF (Lithium Fluoride) according to the RefractiveIndex.info:
    # https://refractiveindex.info/
    elif material == 18:
        B1, B2, B3, C1, C2, C3 = 0.92549, 6.96747, 0, 5.4405376E-3, 1075.1841, 0
        n = sqrt(1 + ((B1 * Lambda_i ** 2) / (Lambda_i ** 2 - C1)) + ((B2 * Lambda_i ** 2) / (Lambda_i ** 2 - C2))
                 + ((B3 * Lambda_i ** 2) / (Lambda_i ** 2 - C3)))

    elif material == 19:  # Cytop
        # According to the AGC chemicals company. Available in:
        # https://www.agc-chemicals.com/jp/en/fluorine/products/cytop/download/index.html
        X = [0.2, 2]

        n_cy = [1.34, 1.34]

        # Method that calculates the refractive indices of the metals by linear interpolation
        # using the points contained in X, n e k_index previously described
        n = complex(interp(Lambda_i, X, n_cy))

    n0 = round(real(n), 5)  # Rounded to five decimal places
    k0 = round(imag(n), 5)  # Rounded to five decimal places

    # Returns the complex refractive index for each material
    return (n0 + k0 * j)


def set_index_custom():
    id_real = float(input(f"Custom refractive index:\n    * Real part: -> "))
    id_imaginary = float(input(f"    * Imaginary part: -> "))
                 
    return complex(id_real, id_imaginary)