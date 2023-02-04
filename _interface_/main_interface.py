from random import gauss, uniform
from scipy.integrate import *
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QWidget
import pandas as pd
from Simulador_Spreeta import Ui_Form
from numpy import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure
import Setting_Layers as sl
import icons_
import sys
import matplotlib.pyplot as plt
import time

REF = 0

class MainWindow(QWidget, Ui_Form):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.setWindowIcon(QtGui.QIcon('icons/LOGO.png'))
        self.setWindowTitle("Simulador Spreeta")

        # Declaração da área de plotagem das figuras:
        
        # Channel sample
        self.figure_sample = plt.figure(dpi=250)
        self.canvas_sample = FigureCanvas(self.figure_sample )
        self.toolbar_sample = NavigationToolbar2QT(self.canvas_sample , self)
        self.verticalLayout_5.addWidget(self.canvas_sample )
        self.gridLayout_2.addWidget(self.toolbar_sample , 2, 0, 1, 4)

        # SPR curve
        self.figure_spr= plt.figure(dpi=250)
        self.canvas_spr = FigureCanvas(self.figure_spr)
        self.toolbar_spr = NavigationToolbar2QT(self.canvas_spr, self)
        self.verticalLayout_3.addWidget(self.canvas_spr)
        self.gridLayout_4.addWidget(self.toolbar_spr, 2, 0, 1,2)


        # Sensorgram
        self.figure_sensorgram = plt.figure(dpi=250)
        self.canvas_sensorgram = FigureCanvas(self.figure_sensorgram )
        self.toolbar_sensorgram = NavigationToolbar2QT(self.canvas_sensorgram , self)
        self.verticalLayout.addWidget(self.canvas_sensorgram )
        self.gridLayout_3.addWidget(self.toolbar_sensorgram , 2, 0, 1, 2)

        # Raw signal
        self.figure_raw = plt.figure(dpi=250)
        self.canvas_raw = FigureCanvas(self.figure_raw)
        self.toolbar_raw = NavigationToolbar2QT(self.canvas_raw , self)
        self.verticalLayout_4.addWidget(self.canvas_raw )
        self.gridLayout_5.addWidget(self.toolbar_raw, 2, 0, 1, 3)

        #Eventos

        # Calibração 
        self.btn_Dry_flow_cell.clicked.connect(self.dryCellCalibration)
        self.btn_Wet_flow_cell.clicked.connect(self.wetCellCalibration)

        # Atualização do indice de refração
        self.refractive_index_analyte.valueChanged.connect(lambda: self.changeSlider(self.refractive_index_analyte, self.refractive_index_slider))
        self.refractive_index_slider.valueChanged.connect(lambda: self.changeSpinBox(self.refractive_index_analyte, self.refractive_index_slider))
        self.cbox_format.currentIndexChanged.connect(self.changeFormat)


    def dryCellCalibration(self):
        self.btn_Dry_flow_cell.setEnabled(False)
        self.btn_Dry_flow_cell.setStyleSheet("background-color: rgb(190, 190, 190); color:rgb(90,90,90)")
        self.setSignalReference()
        
        for i in range(100):
            time.sleep(0.01)
            self.pgsBar_dry_cell.setValue(i+1)
    
    def wetCellCalibration(self):
        self.btn_Wet_flow_cell.setEnabled(False)
        self.btn_Wet_flow_cell.setStyleSheet("background-color: rgb(190, 190, 190); color:rgb(90,90,90)")
        self.setSignalReferenceWet()
        for i in range(100):
            time.sleep(0.01)
            self.pgsBar_wet_cell.setValue(i+1)

    def changeSpinBox(self, spin, slider):
        spin.setValue(float(slider.value()/1E6))
    
    def changeSlider(self, spin, slider):
        slider.setValue(int(spin.value()*1E6))

    def setSignalReference(self):
        dados = pd.read_csv("experimento_H2O.csv", encoding='latin1')
        dados4 = pd.read_csv("signal reference.csv", encoding='latin1')

        pixel = dados['pixel']
        signal_ref = dados4['Reference']

        mn = min(signal_ref)
        mx = max(signal_ref)

        noise_ref = []
        noise_term = []

        theta_i = []

        for i in pixel:
            theta_i.append((3.1522 * 1E-5 * i**2) - (0.0661 * i) + 73.4533)
            a, b = uniform(mn,mx), gauss(mean(0), 0.05)
            noise_ref.append(a)
            noise_term.append(b)

        # Declaração da função gaussiana que modela o espalhamento do Led
        dc = 0.785
        a = 1/(sqrt(2*pi))

        p = lambda w_i: (dc + (a/2.7)* exp(-(((w_i - 845)/30)**2))) * (dc + (a/2.7)* exp(-(((w_i - 853)/15)**2))) * (1.01 + (a/3)* exp(-(((w_i - 845)/42)**2)))

        # Definição de função para calcular o indice de refração
        def set_index_3(lambda_i, material):
            
            ref_index = (sl.set_index(material, lambda_i*1E-9))

            return ref_index

        idx = lambda w_i, material: set_index_3(w_i, material)

        r = lambda theta, lambda_i: self.reflectivity(3, [1, 50*1E-9, 1], [idx(lambda_i, 21), idx(lambda_i, 13), idx(lambda_i,17)], theta*pi/180, lambda_i*1e-9)
        
        def f(w_i, theta):
            return (r(theta, w_i)*p(w_i))

        integrals = [[a, quad(f,  745, 935, args=(a))[0]] for a in linspace(73.6872,65.80896, 128)]

        reflectancia= array(integrals).T[1] / quad(p, 745, 935)[0]

        reflectancia_norm = reflectancia/max(reflectancia)

        reflectancia_norm = (reflectancia_norm * noise_ref)+ noise_term

        arquivo = open('Reference_data.csv','w')
        arquivo.write("Pixel #,Angle,Signal,Noise_Ref")
        for i in range(len(theta_i)):
            arquivo.write(f"\n{pixel[i]},{theta_i[i]},{reflectancia_norm[i]},{noise_ref[i]}")
        arquivo.close()
        
        self.plotReference()
    
    def setSignalReferenceWet(self):
        dados = pd.read_csv("experimento_H2O.csv", encoding='latin1')
        dados4 = pd.read_csv("signal reference.csv", encoding='latin1')

        pixel = dados['pixel']
        signal_ref = dados4['Reference']

        mn = min(signal_ref)
        mx = max(signal_ref)

        noise_ref = []
        noise_term = []

        theta_i = []

        for i in pixel:
            theta_i.append((3.1522 * 1E-5 * i**2) - (0.0661 * i) + 73.4533)
            a, b = uniform(mn,mx), gauss(mean(0), 0.05)
            noise_ref.append(a)
            noise_term.append(b)

        # Declaração da função gaussiana que modela o espalhamento do Led
        dc = 0.785
        a = 1/(sqrt(2*pi))

        p = lambda w_i: (dc + (a/2.7)* exp(-(((w_i - 845)/30)**2))) * (dc + (a/2.7)* exp(-(((w_i - 853)/15)**2))) * (1.01 + (a/3)* exp(-(((w_i - 845)/42)**2)))

        # Definição de função para calcular o indice de refração
        def set_index_3(lambda_i, material):
            
            ref_index = (sl.set_index(material, lambda_i*1E-9))

            return ref_index

        idx = lambda w_i, material: set_index_3(w_i, material)

        r = lambda theta, lambda_i: self.reflectivity(3, [1, 50*1E-9, 1], [idx(lambda_i, 21), idx(lambda_i, 13), idx(lambda_i,16)], theta*pi/180, lambda_i*1e-9)
        
        def f(w_i, theta):
            return (r(theta, w_i)*p(w_i))

        integrals = [[a, quad(f,  745, 935, args=(a))[0]] for a in linspace(73.6872,65.80896, 128)]

        reflectancia= array(integrals).T[1] / quad(p, 745, 935)[0]

        reflectancia_norm = reflectancia/max(reflectancia)

        reflectancia_norm = (reflectancia_norm * noise_ref)+ noise_term

        arquivo = open('Reference_data_wet.csv','w')
        arquivo.write("Pixel #,Angle,Signal,Noise_Ref")
        for i in range(len(theta_i)):
            arquivo.write(f"\n{pixel[i]},{theta_i[i]},{reflectancia_norm[i]},{noise_ref[i]}")
        arquivo.close()
        
        self.plotReferenceWet()

    def reflectivity(self, nLayers, d,  index, theta_i, wavelenght):

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

        a = (abs(r_TM) ** 2)
        
        return a # Reflectance - TM polarization

    def plotReference(self):
        global REF
        REF = 1
        
        font=dict(size=4, family="Sans-Serif")
        plt.rc('font', **font)

        plt.subplots_adjust(left=0.125,
            bottom=0.150, 
            right=0.960, 
            top=0.950, 
            wspace=0.1, 
            hspace=0.2)

        signal_ref_txt = pd.read_csv('Reference_data.csv', encoding='latin1')
        pixel_ = signal_ref_txt['Pixel #']
        angle_ = signal_ref_txt['Angle']
        signal_ = signal_ref_txt['Signal']
        noise_ref = signal_ref_txt['Noise_Ref']
        
        if self.cbox_format.currentText() == "Signal vs. Pixel":
            self.figure_raw.clear()
            plt.plot(pixel_, signal_, label='Signal', linewidth=0.5)
            plt.grid(alpha=0.5)
        else:
            self.figure_raw.clear()
            plt.plot(angle_, (signal_/noise_ref), label= 'Signal Ref TXT', linewidth=0.5)
            plt.grid(alpha=0.5)
        
        self.canvas_raw.draw()
    
    def plotReferenceWet(self):
        global REF
        REF = 2
        
        font=dict(size=4, family="Sans-Serif")
        plt.rc('font', **font)

        plt.subplots_adjust(left=0.125,
            bottom=0.150, 
            right=0.960, 
            top=0.950, 
            wspace=0.1, 
            hspace=0.2)

        signal_ref_txt = pd.read_csv('Reference_data_wet.csv', encoding='latin1')
        pixel_ = signal_ref_txt['Pixel #']
        angle_ = signal_ref_txt['Angle']
        signal_ = signal_ref_txt['Signal']
        noise_ref = signal_ref_txt['Noise_Ref']
        
        if self.cbox_format.currentText() == "Signal vs. Pixel":
            self.figure_raw.clear()
            plt.plot(pixel_, signal_, label='Signal', linewidth=0.5)
            plt.grid(alpha=0.5)
        else:
            self.figure_raw.clear()
            plt.plot(angle_, (signal_/noise_ref), label= 'Signal Ref TXT', linewidth=0.5)
            plt.grid(alpha=0.5)
        
        self.canvas_raw.draw()


    def changeFormat(self):
        if REF == 1:
            self.plotReference()
        elif REF == 2:
            self.plotReferenceWet()
        else:
            pass


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    Widget = MainWindow()
    Widget.show()
    app.exec()