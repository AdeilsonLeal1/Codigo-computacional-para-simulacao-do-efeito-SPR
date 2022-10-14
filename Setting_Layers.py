##########     Setting_Layers     ##########
# module that includes the functions needed to build the layers of the structure
# setLayers -> sets the material for each layer


d = []          # Layers thickness
material = []   # The composition of the layers


def setLayers(n):
    # Choosing layer material
    for layer in range(n-1):
        if layer == 0:
            print("\n=====================================================\n"
                    "==========    Set Layers Characteristics   ==========\n"
                    "=====================================================\n\n"
                    "==========  1st Layer - Optical substrate  ==========\n")
            d.append(1)
            material.append((int(input(f"\n1 - BK7   2 - Sílica   3 - N-F2   4 - Safira sintética(Al2O3)"
                                    f"\n5 - SFL6  6 - FK51A    7 - N-SF14 8 - Acrilico SUVT   "
                                    f"\n\nMaterial -> "))))
        else:
            print(f"\n==========            {layer+1}st Layer            ==========\n")
            material.append((int(input(f"\n 1 - BK7     2 - Sílica      3 - N-F2       4 - Safira sintética("
                                    f"Al2O3) "
                                    f"\n 5 - SFL6    6 - FK51A       7 - N-SF14     8 - Acrilico SUVT"
                                    f"\n 9 - PVA    10 - Glicerina  11 - Quartzo   12 - Aluminio"
                                    f"\n13 - Ouro   14 - Prata      15 - Cobre     16 - Água "
                                    f"\n17 - Ar    \n\nMaterial -> "))))

           
            d.append(float(input("Thickness (nm): ")) * 1e-9)

    return d, material
