# PyQt5 modules
from PyQt5.QtWidgets import QMainWindow, QInputDialog, QListWidgetItem
from PyQt5.QtCore import Qt

# Project modules
from src.ui.mainwindow import Ui_MainWindow
from src.Filters2 import FilterClass
import scipy.signal as ss
import matplotlib.pyplot as plt
import numpy as np
from src.plottingClasses import BodePlot
from src.PlantillaClass import PlantillaClass

class MainWindow(QMainWindow, Ui_MainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)

        self.listaFiltros = []
        self.listaPlantillas = []

        # Configuración gráfica

        # Gráficos de atenuación
        self.atenMagPlot = BodePlot(parent = self.atenMagPlotBox)
        self.atenMagPlotBox.layout().addWidget(self.atenMagPlot.navToolBar)
        self.atenMagPlotBox.layout().addWidget(self.atenMagPlot) #importante setear tamaños mínimos desde QTDesigner, sino no anda

        self.atenFasePlot = BodePlot(parent=self.atenFasePlotBox)
        self.atenFasePlotBox.layout().addWidget(self.atenFasePlot.navToolBar)
        self.atenFasePlotBox.layout().addWidget(self.atenFasePlot)

        # Gráficos de ganancia
        self.ganMagPlot = BodePlot(parent = self.ganMagPlotBox)
        self.ganMagPlotBox.layout().addWidget(self.ganMagPlot.navToolBar)
        self.ganMagPlotBox.layout().addWidget(self.ganMagPlot)

        self.ganFasePlot = BodePlot(parent=self.ganFasePlotBox)
        self.ganFasePlotBox.layout().addWidget(self.ganFasePlot.navToolBar)
        self.ganFasePlotBox.layout().addWidget(self.ganFasePlot)

        # Gráficos de polos y ceros

        # Configuración de las pestañas y clicks
        self.crearPlantillaButton.clicked.connect(self.crearPlantilla)
        self.crearFiltroButton.clicked.connect(self.crearFiltro)


        self.borrarFiltrosButton.clicked.connect(self.deleteFiltros)
        self.borrarTodosFiltrosButton.clicked.connect(self.deleteAllFiltros)
        self.borrarPlantillasButton.clicked.connect(self.deletePlantillas)
        self.borrarTodasPlantillasButton.clicked.connect(self.deleteAllPlantillas)


    def crearPlantilla(self):
        newPlantilla = PlantillaClass()

        text = self.tipoFiltroComboBox.currentText()

        if text == "Pasa Bajos":
            fp = float(self.fpLP.toPlainText())
            fa = float(self.faLP.toPlainText())
            ap= float(self.apVal.toPlainText())
            aa = float(self.aaVal.toPlainText())


            newPlantilla.crearPasaBajos(fp, fa, ap, aa)

        elif text == "Pasa Altos":
            fp = float(self.fpHP.toPlainText())
            fa = float(self.faHP.toPlainText())
            ap = float(self.apVal.toPlainText())
            aa = float(self.aaVal.toPlainText())

            newPlantilla.crearPasaAltos(fp, fa, ap, aa)

        elif text == "Pasa Banda":
            tipo = self.BPTypeComboBox.currentText()
            if tipo == "Anchos de Banda":
                dfp = float(self.dFpBP.toPlainText())
                dfa = float(self.dFaBP.toPlainText())
                ap = float(self.apVal.toPlainText())
                aa = float(self.aaVal.toPlainText())

                newPlantilla.crearPasaBandaBW(dfp, dfa, ap, aa)

            elif tipo == "Frecuencias":
                fo = float(self.f0BP.toPlainText())
                fpx = float(self.fpXBP.toPlainText())
                fpy = float(self.fpYBP.toPlainText())
                fax = float(self.faXBP.toPlainText())
                fay = float(self.faYBP.toPlainText())
                ap = float(self.apVal.toPlainText())
                aa = float(self.aaVal.toPlainText())

                newPlantilla.crearPasaBandaFreq(fo, fpx, fpy, fax, fay, ap, aa)

        elif text == "Rechaza Banda":
            tipo = self.BSTypeComboBox.currentText()
            if tipo == "Anchos de Banda":
                dfp = float(self.dFpBP.toPlainText())
                dfa = float(self.dFaBP.toPlainText())
                ap = float(self.apVal.toPlainText())
                aa = float(self.aaVal.toPlainText())

                newPlantilla.crearRechazaBandaBW(dfp, dfa, ap, aa)

            elif tipo == "Frecuencias":
                fo = float(self.f0BS.toPlainText())
                fpx = float(self.fpXBS.toPlainText())
                fpy = float(self.fpYBS.toPlainText())
                fax = float(self.faXBS.toPlainText())
                fay = float(self.faYBS.toPlainText())
                ap = float(self.apVal.toPlainText())
                aa = float(self.aaVal.toPlainText())

                newPlantilla.crearRechazaBandaFreq(fo, fpx, fpy, fax, fay, ap, aa)


        plantillaName, ok = QInputDialog.getText(self, "Agregar Plantilla", 'Nombre la plantilla: ')
        if not ok:
            return
        if len(plantillaName) < 1:
            plantillaName = "Plantilla " + str(len(self.listaPlantillas))

        newPlantilla.crearNombre(plantillaName)

        self.listaPlantillas.append(newPlantilla)

        item = QListWidgetItem(plantillaName)
        item.setCheckState(Qt.Checked)

        self.listPlantillasWidget.addItem(item)

    def crearFiltro(self):
        text = self.tipoAproxComboBox.currentText()

        if text == "Butterworth":
            text = "butter"
        elif text == "Chebyshev I":
            text = "cheby1"
        elif text == "Chebyshev II (inverso)":
            text = "cheby2"
        elif text == "Cauer (elíptico)":
            text = "ellip"

        for i in range(len(self.listaPlantillas)):
            if self.listPlantillasWidget.item(i).checkState() == 2:
                filtroActual = FilterClass()
                p = self.listaPlantillas[i]

                if p.tipoPlantilla == "LP":
                    filtroActual.getLPTransferFunction(p.fp, p.fa, p.ap, p.aa, text, 0)

                elif p.tipoPlantilla == "HP":
                    filtroActual.getHPTransferFunction(p.fp, p.fa, p.ap, p.aa, text, 0)

                elif p.tipoPlantilla == "BPF":
                    filtroActual.getBPTransferFunctionFreq([p.fpx,p.fpy], [p.fax, p.fay], p.ap, p.aa, text, 0)

                elif p.tipoPlantilla == "BPBW":
                    filtroActual.getBPTransferFunctionBW(p.f0, p.dfp, p.dfa, p.ap, p.aa, text, 0)

                elif p.tipoPlantilla == "BSF":
                    filtroActual.getBSTransferFunctionFreq([p.fpx, p.fpy], [p.fax, p.fay], p.ap, p.aa, text, 0)

                elif p.tipoPlantilla == "BSBW":
                    filtroActual.getBSTransferFunctionBW(p.f0, p.dfp, p.dfa, p.ap, p.aa, text, 0)



                filtroName, ok = QInputDialog.getText(self, "Agregar Filtro", 'Nombre del filtro: ')
                if not ok:
                    return
                if len(filtroName) < 1:
                    filtroName = "Filtro " + str(len(self.listaFiltros))

                filtroActual.crearNombre(filtroName)

                self.listaFiltros.append(filtroActual)

                item = QListWidgetItem(filtroName)
                item.setCheckState(Qt.Checked)

                self.listFiltrosWidget.addItem(item)

    def deletePlantillas(self):
        numberOfFunctions = len(self.listaPlantillas)
        for i in range(numberOfFunctions):
            if self.listPlantillasWidget.item(numberOfFunctions-1-i).checkState() == 2:
                self.listPlantillasWidget.takeItem(numberOfFunctions-1-i)
                self.listaPlantillas.pop(numberOfFunctions-1-i)

    def deleteAllPlantillas(self):
        numberOfFunctions = len(self.listaPlantillas)
        for i in range(numberOfFunctions):
            self.listPlantillasWidget.takeItem(numberOfFunctions - 1 - i)
            self.listaPlantillas.pop(numberOfFunctions - 1 - i)

    def deleteFiltros(self):
        numberOfFunctions = len(self.listaFiltros)
        for i in range(numberOfFunctions):
            if self.listFiltrosWidget.item(numberOfFunctions-1-i).checkState() == 2:
                self.listFiltrosWidget.takeItem(numberOfFunctions-1-i)
                self.listaFiltros.pop(numberOfFunctions-1-i)

    def deleteAllFiltros(self):
        numberOfFunctions = len(self.listaFiltros)
        for i in range(numberOfFunctions):
            self.listFiltrosWidget.takeItem(numberOfFunctions - 1 - i)
            self.listaFiltros.pop(numberOfFunctions - 1 - i)











