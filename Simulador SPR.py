# This computational code allows calculating and plotting the reflectivity curve of a SPR (Surface Plasmon Resonance) 
# sensor based on the Kretschmann configuration, operating both in angular interrogation mode (AIM) and 
# in wavelength interrogation mode (WIM).

from cmath import nan, pi
import numpy as np 
import os

ACF = (pi/180)  # Angle Conversion Factor. 

mod_int = nan
lambda_i = nan
theta_i = nan

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
mod_int = int(input("Interrogation mode: \n     1 - Angular Interrogation Mode (AIM)\n     2 - Wavelength Interrogation Mode (WIM)\n >> "))

if mod_int == 1:
    #  Setting Angular interrogation Mode:
    #       - Incident wavelength and angular range (ex. 60 deg to 80 deg);
    print("## Angular Interrogation Mode\n"
          "\n## Setting angular interrogation mode:\n\n ")

    # Incident wavelength
    lambda_i = float(input("Incident wavelength (nm): "))
    
    # Starting and ending angle
    a1 = float(input("Starting angle value (degress): "))*ACF
    a2 = float(input("Ending angle value (degress): "))*ACF
    
    theta_i = np.arange(a1, a2, 0.001*ACF) # Array with angles of incidence

    print(f"len(theta_i) = {len(theta_i)}\n theta_i = {theta_i} \nIncidet wavelength = {lambda_i} nm \n a1 = {a1/ACF} deg. and a2 = {a2/ACF} deg.")
    
    
    #  Defining the structure
    n_layers = int(input("Set the number of layers in your structure:\n N = "))






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
    
    lambda_i = np.arange(a1, a2, 0.1) # Array with wavelength of incidence

    print(f"len(lambda_i) = {len(lambda_i)}\n lambda_i= {lambda_i} \nAngle of incidence= {theta_i/ACF} deg. \n a1 = {a1} nm and a2 = {a2} nm")


else:
    print("Invalid!!")
