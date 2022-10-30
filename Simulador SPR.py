# This computational code allows calculating and plotting the reflectivity curve of a SPR (Surface Plasmon Resonance)
# sensor based on the Kretschmann configuration, operating both in angular interrogation mode (AIM) and
# in wavelength interrogation mode (WIM).

from cmath import nan, pi
import numpy as np
import os
import Setting_Layers as sl
import Reflectivity as ref
import tools

ACF = (pi/180)  # Angle Conversion Factor.

mod_int = nan   # Interrogation mode
lambda_i = nan  # Incident wavelength
theta_i = nan   # Angle of incidence
thickness = []  # Thickness of each layer
material = []   # Material of each layer
R_Tm = []       # Stores lists with reflectivity curves for each of the interactions
resonance_point = nan   # Variable that stores the SPR resonance point
ref_index = []  # Refractive index of the layers

# screen cleaning
if os.name == 'nt':
    os.system('cls')
else:
    os.system('clear')

# introduction
print("\nThis computacional code allows calculating calculating and plotting the reflectivity curve of a SPR sensor \n"
      "based on the Kretschmann configuration, operating both in angular interrogation mode (AIM) and "
      "in wavelength \ninterrogation mode (WIM).\n")
x = input("Press to continue... \n>>>")
if os.name == 'nt':
    os.system('cls')
else:
    os.system('clear')

# Selecting interrogation mode
mod_int = int(input(
    "Interrogation mode: \n     1 - Angular Interrogation Mode (AIM)\n     2 - Wavelength Interrogation Mode (WIM)\n >> "))

if mod_int == 1:
    #  Setting Angular interrogation Mode:
    #       - Incident wavelength and angular range (ex. 60 deg to 80 deg);
    print("## Angular Interrogation Mode\n"
          "\n## Setting angular interrogation mode:\n\n ")

    # Incident wavelength
    lambda_i = float(input("Incident wavelength (nm): ")) * 1e-9

    # Starting and ending angle
    a1 = float(input("Starting angle value (degress): "))*ACF
    a2 = float(input("Ending angle value (degress): "))*ACF

    theta_i = np.arange(a1, a2, 0.001*ACF)  # Array with angles of incidence

    #  Defining the structure
    n_layers = int(
        input("\nSet the number of layers in your structure:\n     N = "))

    for layer in range(n_layers):
        d, m = sl.setLayers(layer)  # Sets the thickness and material for each layer
        thickness.append(d)
        material.append(m)

        if layer == 0:
            if material[layer] == 9:
                ref_index.append(sl.set_index_custom())
            else:
                ref_index.append(sl.set_index(material[layer],lambda_i))
        else:
            if material[layer] == 20:
                ref_index.append(sl.set_index_custom())
            else:
                ref_index.append(sl.set_index(material[layer],lambda_i))


    #n_subs = int(input("\nAnalyze sensitivity for how many interactions?\n >> "))
    n_subs = 1

    #list_analyte = sl.setAnalyte(n_subs,ref_index[-1])

    for i in range(n_subs):
        r_TM = []      # Local variable that temporarily stores reflectivity values for each of the interactions

        #ref_index[-1] = complex(list_analyte[i])

        for t in range(len(theta_i)):
            r_TM.append(ref.Reflectivity(
                n_layers, thickness, ref_index, theta_i[t], lambda_i,))

        resonance_point = tools.point_SPR(r_TM, theta_i, mod_int)

        R_Tm.append(r_TM)

        print(f"\n\nResonance Angle = {resonance_point:.4f}°\n\n")
    tools.plot(theta_i, R_Tm, resonance_point, mod_int)

elif mod_int == 2:
    #  Setting Wavelength interrogation Mode:
    #       - Angle of incidence and spectral range (ex. 400 nm to 1500 nm);
    print("## Wavelength Interrogation Mode\n"
          "\n## Setting wavelength interrogation mode:\n\n ")

    # Angle of incidence
    theta_i = float(input("Angle of incidence (degress): "))*ACF

    # Starting and ending wavelength
    a1 = float(input("Starting wavelength value (nm): ")) * 1e-9
    a2 = float(input("Ending wavelength value (nm): ")) * 1e-9

    lambda_i = np.arange(a1, a2, 1e-11)  # Array with wavelength of incidence

    #  Defining the structure
    n_layers = int(
        input("\nSet the number of layers in your structure:\n     N = "))
    
    for layer in range(n_layers):
        d, m = sl.setLayers(layer)  # Sets the thickness and material for each layer
        thickness.append(d)
        material.append(m)

    #n_subs = int(input("\nAnalyze sensitivity for how many interactions?\n >> "))
    n_subs = 1

    #list_analyte = sl.setAnalyte(n_subs,ref_index[-1])

    for i in range(n_subs):
        r_TM = []      # Local variable that temporarily stores reflectivity values for each of the interactions
        #ref_index[-1] = complex(list_analyte[i])
        for t in range(len(lambda_i)):
            ref_index = []
            for layer in range(n_layers):
                ref_index.append(sl.set_index(material[layer],lambda_i[t]))

            r_TM.append(ref.Reflectivity(
                n_layers, thickness, ref_index, theta_i, lambda_i[t],))

        resonance_point = tools.point_SPR(r_TM, lambda_i, mod_int)

        R_Tm.append(r_TM)

        print(f"\n\nResonance Wavelength = {resonance_point:.4f} nm\n\n")

    tools.plot(lambda_i, R_Tm, resonance_point, mod_int)


else:
    print("Invalid!!")
