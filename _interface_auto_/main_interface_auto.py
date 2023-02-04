from scipy.signal import butter, lfilter, freqz, filtfilt
from scipy.integrate import *
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QWidget
import pandas as pd
from Simulador_Spreeta_auto import Ui_Form
from numpy import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
import Setting_Layers as sl
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
        font=dict(size=4, family="Sans-Serif")
        plt.rc('font', **font)

        # Sample
        self.figure_sample= plt.figure(dpi=250)
        self.canvas_sample = FigureCanvas(self.figure_sample)
        self.toolbar_sample = NavigationToolbar2QT(self.canvas_sample, self)
        self.verticalLayout_5.addWidget(self.canvas_sample)
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

        # Atualização do formato
        self.cbox_format.currentIndexChanged.connect(self.changeFormat)

        # Iniciando monitoramento
        self.btn_start.clicked.connect(self.startMonitoring)


    def dryCellCalibration(self):
        self.btn_Dry_flow_cell.setEnabled(False)
        self.btn_Dry_flow_cell.setStyleSheet("background-color: rgb(190, 190, 190); color:rgb(90,90,90)")
        self.setSignalReference()
        
        for i in range(100):
            time.sleep(0.02)
            self.pgsBar_dry_cell.setValue(i+1)
    
    def wetCellCalibration(self):
        self.btn_Wet_flow_cell.setEnabled(False)
        self.btn_Wet_flow_cell.setStyleSheet("background-color: rgb(190, 190, 190); color:rgb(90,90,90)")
        self.setSignalReferenceWet()
        for i in range(100):
            time.sleep(0.02)
            self.pgsBar_wet_cell.setValue(i+1)

    def setSignalReference(self):
        self.plotReference()
    
    def setSignalReferenceWet(self):      
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

        plt.subplots_adjust(top=0.939,
                            bottom=0.218,
                            left=0.125,
                            right=0.969,
                            hspace=0.2,
                            wspace=0.2)
        
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
            plt.plot(angle_, (signal_/noise_ref), label='Signal', linewidth=0.5)
            plt.grid(alpha=0.5)
        
        self.canvas_raw.draw()
    
    def plotReferenceWet(self):
        global REF
        REF = 2

        plt.subplots_adjust(top=0.939,
                            bottom=0.218,
                            left=0.125,
                            right=0.969,
                            hspace=0.2,
                            wspace=0.2)
        
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
            plt.plot(angle_, (signal_/noise_ref), label='Signal', linewidth=0.5)
            plt.grid(alpha=0.5)
        
        self.canvas_raw.draw()

    def changeFormat(self):
        if REF == 1:
            self.plotReference()
        elif REF == 2:
            self.plotReferenceWet()
        else:
            pass

    def startMonitoring(self):

        signal_ref_txt = pd.read_csv('Reference_data_wet.csv', encoding='latin1')
        pixel_ = signal_ref_txt['Pixel #']
        angle_ = signal_ref_txt['Angle']
        signal_ = signal_ref_txt['Signal']
        noise_ref = signal_ref_txt['Noise_Ref']
        
        plt.subplots_adjust(top=0.96,
                            bottom=0.16,
                            left=0.125,
                            right=0.97,
                            hspace=0.2,
                            wspace=0.2)
        
        self.figure_sample.clear()
        self.ax = self.figure_sample.add_subplot()
        self.ax.plot(angle_, (signal_/noise_ref), label= 'Signal Ref TXT', linewidth=0.5)
        self.ax.grid(alpha=0.5)
            
        self.canvas_sample.draw()
        time.sleep(2)
        self.figure_sample.clear()

        signal_sensor = pd.read_csv('Sensor_data.txt', encoding='latin1', delimiter='\t')
        angle_ = []
        signal_ = []
        
        for i in range(len(signal_sensor['time(s)'])):
            a = signal_sensor['Angle'][i].replace('[', '').replace(']', '').replace(' ', '')
            s = signal_sensor['SPR Curve'][i].replace('[', '').replace(']', '').replace(' ', '')

            a = a.split(',')
            angle_ = [float(j) for j in a]
            s = s.split(',')
            signal_ = [float(j) for j in s]
            
            signal_filtered = self.butter_lowpass_filter(data=signal_, cutoff=2, fs=30, order=5 )

            self.plotSprCurve(angle_, signal_filtered)
            time.sleep(2)

    
    
    def butter_lowpass_filter(self, data, cutoff, fs, order):
        normal_cutoff = cutoff / (0.5 * fs)
        # Get the filter coefficients 
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        y = filtfilt(b, a, data)
        return y
    
    def plotSprCurve(self,_angle, _signal):       
        self.figure_spr.clear()
        plt.subplots_adjust(top=0.96,
                            bottom=0.16,
                            left=0.125,
                            right=0.97,
                            hspace=0.2,
                            wspace=0.2)
        self.ax1 = self.figure_spr.add_subplot()
        self.ax1.plot(_angle, _signal, linewidth=0.5)
        self.ax1.grid(alpha=0.5)
            
        self.canvas_spr.draw()
        


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    Widget = MainWindow()
    Widget.show()
    app.exec()