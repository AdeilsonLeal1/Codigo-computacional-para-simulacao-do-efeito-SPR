# This computational code allows calculating and plotting the reflectivity curve of a SPR (Surface Plasmon Resonance)
# sensor based on the Kretschmann configuration, operating both in angular interrogation mode (AIM) and
# in wavelength interrogation mode (WIM).

from cmath import nan, pi
import numpy as np
import os
import Setting_Layers as sl
import Reflectivity as ref
import matplotlib.pyplot as plt
import tools

ACF = (pi/180)  # Angle Conversion Factor.

mod_int = nan   # Interrogation mode
lambda_i = nan  # Incident wavelength
theta_i = nan   # Angle of incidence
R_Tm = []       # Stores lists with reflectivity curves for each of the interactions

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

    d, ref_index = sl.setLayers(n_layers, lambda_i)

    n_subs = int(
        input("\nHow many substances to analyze in the recipe?\n >> "))

    list_analyte = sl.setAnalyte(n_subs)

    for i in range(n_subs):
        r_TM = []      # Local variable that temporarily stores reflectivity values for each of the interactions

        ref_index[-1] = complex(list_analyte[i])

        for t in range(len(theta_i)):
            r_TM.append(ref.Reflectivity(n_layers, d, ref_index, theta_i[t], lambda_i,))




        R_Tm.append(r_TM)

        print(f"d = {d}\nRefractive Index = {ref_index}")

    y = 'Angle'
    c = 'Degrees'
    n = chr(952)
    ax_x = theta_i / ACF

    fig0, ax_TM = plt.subplots()
    ax_TM.plot(ax_x, (R_Tm[0]))
    ax_TM.set_title("Reflectivity vs. Angle of Incidence", loc='center', pad='6')
    ax_TM.set(xlabel=f'Incidence {y} ({c})', ylabel='Reflectivity')
    ax_TM.grid()
    ax_TM.set_yticks(np.arange(0, 1.20, 0.20))

    fig1, ax_TM2 = plt.subplots()
    for i in range(n_subs):
        ax_TM2.plot(ax_x, (R_Tm[i]))
    ax_TM2.set(xlabel=f'Incidence {y} ({c})', ylabel='Reflectivity')
    ax_TM2.grid()
    ax_TM.set_yticks(np.arange(0, 1.20, 0.20))

    #plt.yticks(np.arange(0, 1.20, 0.20))
    plt.show()

elif mod_int == 2:
    #  Setting Wavelength interrogation Mode:
    #       - Angle of incidence and spectral range (ex. 400 nm to 1500 nm);
    print("## Wavelength Interrogation Mode\n"
          "\n## Setting wavelength interrogation mode:\n\n ")

    # Angle of incidence
    theta_i = float(input("Angle of incidence (degress): "))*ACF

    # Starting and ending wavelength
    a1 = float(input("Starting wavelength value (nm): "))
    a2 = float(input("Ending wavelength value (nm): "))

    lambda_i = np.arange(a1, a2, 0.1)  # Array with wavelength of incidence

    print(
        f"len(lambda_i) = {len(lambda_i)}\n lambda_i= {lambda_i} \nAngle of incidence= {theta_i/ACF} deg. \n a1 = {a1} nm and a2 = {a2} nm")


else:
    print("Invalid!!")
