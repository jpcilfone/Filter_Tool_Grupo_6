# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'designer\mainwindow.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(954, 755)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName("tabWidget")
        self.tab_5 = QtWidgets.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.tab_5)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.tabWidget_2 = QtWidgets.QTabWidget(self.tab_5)
        self.tabWidget_2.setMinimumSize(QtCore.QSize(500, 0))
        self.tabWidget_2.setObjectName("tabWidget_2")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.tab)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.atenMagPlotBox = QtWidgets.QGroupBox(self.tab)
        self.atenMagPlotBox.setMinimumSize(QtCore.QSize(300, 300))
        self.atenMagPlotBox.setTitle("")
        self.atenMagPlotBox.setObjectName("atenMagPlotBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.atenMagPlotBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.gridLayout_4.addWidget(self.atenMagPlotBox, 0, 0, 1, 1)
        self.atenFasePlotBox = QtWidgets.QGroupBox(self.tab)
        self.atenFasePlotBox.setMinimumSize(QtCore.QSize(300, 300))
        self.atenFasePlotBox.setTitle("")
        self.atenFasePlotBox.setObjectName("atenFasePlotBox")
        self.gridLayout_10 = QtWidgets.QGridLayout(self.atenFasePlotBox)
        self.gridLayout_10.setObjectName("gridLayout_10")
        self.gridLayout_4.addWidget(self.atenFasePlotBox, 1, 0, 1, 1)
        self.tabWidget_2.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.tab_2)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.ganMagPlotBox = QtWidgets.QGroupBox(self.tab_2)
        self.ganMagPlotBox.setTitle("")
        self.ganMagPlotBox.setObjectName("ganMagPlotBox")
        self.gridLayout_8 = QtWidgets.QGridLayout(self.ganMagPlotBox)
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.gridLayout_6.addWidget(self.ganMagPlotBox, 0, 0, 1, 1)
        self.ganFasePlotBox = QtWidgets.QGroupBox(self.tab_2)
        self.ganFasePlotBox.setTitle("")
        self.ganFasePlotBox.setObjectName("ganFasePlotBox")
        self.gridLayout_11 = QtWidgets.QGridLayout(self.ganFasePlotBox)
        self.gridLayout_11.setObjectName("gridLayout_11")
        self.gridLayout_6.addWidget(self.ganFasePlotBox, 1, 0, 1, 1)
        self.tabWidget_2.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.tab_3)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.polosCerosPlot = QtWidgets.QGroupBox(self.tab_3)
        self.polosCerosPlot.setTitle("")
        self.polosCerosPlot.setObjectName("polosCerosPlot")
        self.gridLayout_12 = QtWidgets.QGridLayout(self.polosCerosPlot)
        self.gridLayout_12.setObjectName("gridLayout_12")
        self.gridLayout_7.addWidget(self.polosCerosPlot, 0, 0, 1, 1)
        self.tabWidget_2.addTab(self.tab_3, "")
        self.gridLayout_3.addWidget(self.tabWidget_2, 0, 3, 10, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(self.tab_5)
        self.groupBox_2.setMaximumSize(QtCore.QSize(200, 16777215))
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label = QtWidgets.QLabel(self.groupBox_2)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)
        self.tipoFiltroComboBox = QtWidgets.QComboBox(self.groupBox_2)
        self.tipoFiltroComboBox.setObjectName("tipoFiltroComboBox")
        self.tipoFiltroComboBox.addItem("")
        self.tipoFiltroComboBox.addItem("")
        self.tipoFiltroComboBox.addItem("")
        self.tipoFiltroComboBox.addItem("")
        self.gridLayout_2.addWidget(self.tipoFiltroComboBox, 0, 1, 1, 1)
        self.stackedWidget_2 = QtWidgets.QStackedWidget(self.groupBox_2)
        self.stackedWidget_2.setMaximumSize(QtCore.QSize(16777215, 250))
        self.stackedWidget_2.setObjectName("stackedWidget_2")
        self.page_3 = QtWidgets.QWidget()
        self.page_3.setObjectName("page_3")
        self.gridLayout_16 = QtWidgets.QGridLayout(self.page_3)
        self.gridLayout_16.setObjectName("gridLayout_16")
        self.label_53 = QtWidgets.QLabel(self.page_3)
        self.label_53.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_53.setObjectName("label_53")
        self.gridLayout_16.addWidget(self.label_53, 0, 0, 1, 1)
        self.fpLP = QtWidgets.QPlainTextEdit(self.page_3)
        self.fpLP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.fpLP.setObjectName("fpLP")
        self.gridLayout_16.addWidget(self.fpLP, 0, 1, 1, 1)
        self.label_27 = QtWidgets.QLabel(self.page_3)
        self.label_27.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_27.setObjectName("label_27")
        self.gridLayout_16.addWidget(self.label_27, 1, 0, 1, 1)
        self.faLP = QtWidgets.QPlainTextEdit(self.page_3)
        self.faLP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.faLP.setObjectName("faLP")
        self.gridLayout_16.addWidget(self.faLP, 1, 1, 1, 1)
        self.stackedWidget_2.addWidget(self.page_3)
        self.page_4 = QtWidgets.QWidget()
        self.page_4.setObjectName("page_4")
        self.gridLayout_17 = QtWidgets.QGridLayout(self.page_4)
        self.gridLayout_17.setObjectName("gridLayout_17")
        self.label_55 = QtWidgets.QLabel(self.page_4)
        self.label_55.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_55.setObjectName("label_55")
        self.gridLayout_17.addWidget(self.label_55, 0, 0, 1, 1)
        self.fpHP = QtWidgets.QPlainTextEdit(self.page_4)
        self.fpHP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.fpHP.setObjectName("fpHP")
        self.gridLayout_17.addWidget(self.fpHP, 0, 1, 1, 1)
        self.label_54 = QtWidgets.QLabel(self.page_4)
        self.label_54.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_54.setObjectName("label_54")
        self.gridLayout_17.addWidget(self.label_54, 1, 0, 1, 1)
        self.faHP = QtWidgets.QPlainTextEdit(self.page_4)
        self.faHP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.faHP.setObjectName("faHP")
        self.gridLayout_17.addWidget(self.faHP, 1, 1, 1, 1)
        self.stackedWidget_2.addWidget(self.page_4)
        self.page_7 = QtWidgets.QWidget()
        self.page_7.setObjectName("page_7")
        self.gridLayout_18 = QtWidgets.QGridLayout(self.page_7)
        self.gridLayout_18.setObjectName("gridLayout_18")
        self.BPTypeComboBox = QtWidgets.QComboBox(self.page_7)
        self.BPTypeComboBox.setObjectName("BPTypeComboBox")
        self.BPTypeComboBox.addItem("")
        self.BPTypeComboBox.addItem("")
        self.gridLayout_18.addWidget(self.BPTypeComboBox, 0, 0, 1, 2)
        self.stackedWidget_4 = QtWidgets.QStackedWidget(self.page_7)
        self.stackedWidget_4.setObjectName("stackedWidget_4")
        self.page_8 = QtWidgets.QWidget()
        self.page_8.setObjectName("page_8")
        self.gridLayout_19 = QtWidgets.QGridLayout(self.page_8)
        self.gridLayout_19.setObjectName("gridLayout_19")
        self.label_63 = QtWidgets.QLabel(self.page_8)
        self.label_63.setObjectName("label_63")
        self.gridLayout_19.addWidget(self.label_63, 0, 0, 1, 1)
        self.f0BP = QtWidgets.QPlainTextEdit(self.page_8)
        self.f0BP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.f0BP.setObjectName("f0BP")
        self.gridLayout_19.addWidget(self.f0BP, 0, 1, 1, 1)
        self.label_62 = QtWidgets.QLabel(self.page_8)
        self.label_62.setObjectName("label_62")
        self.gridLayout_19.addWidget(self.label_62, 1, 0, 1, 1)
        self.dFpBP = QtWidgets.QPlainTextEdit(self.page_8)
        self.dFpBP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.dFpBP.setObjectName("dFpBP")
        self.gridLayout_19.addWidget(self.dFpBP, 1, 1, 1, 1)
        self.label_61 = QtWidgets.QLabel(self.page_8)
        self.label_61.setObjectName("label_61")
        self.gridLayout_19.addWidget(self.label_61, 2, 0, 1, 1)
        self.dFaBP = QtWidgets.QPlainTextEdit(self.page_8)
        self.dFaBP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.dFaBP.setObjectName("dFaBP")
        self.gridLayout_19.addWidget(self.dFaBP, 2, 1, 1, 1)
        self.stackedWidget_4.addWidget(self.page_8)
        self.page_9 = QtWidgets.QWidget()
        self.page_9.setObjectName("page_9")
        self.gridLayout_20 = QtWidgets.QGridLayout(self.page_9)
        self.gridLayout_20.setObjectName("gridLayout_20")
        self.label_58 = QtWidgets.QLabel(self.page_9)
        self.label_58.setObjectName("label_58")
        self.gridLayout_20.addWidget(self.label_58, 0, 0, 1, 1)
        self.fpXBP = QtWidgets.QPlainTextEdit(self.page_9)
        self.fpXBP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.fpXBP.setObjectName("fpXBP")
        self.gridLayout_20.addWidget(self.fpXBP, 0, 1, 1, 1)
        self.label_60 = QtWidgets.QLabel(self.page_9)
        self.label_60.setObjectName("label_60")
        self.gridLayout_20.addWidget(self.label_60, 1, 0, 1, 1)
        self.fpYBP = QtWidgets.QPlainTextEdit(self.page_9)
        self.fpYBP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.fpYBP.setObjectName("fpYBP")
        self.gridLayout_20.addWidget(self.fpYBP, 1, 1, 1, 1)
        self.label_64 = QtWidgets.QLabel(self.page_9)
        self.label_64.setObjectName("label_64")
        self.gridLayout_20.addWidget(self.label_64, 2, 0, 1, 1)
        self.faXBP = QtWidgets.QPlainTextEdit(self.page_9)
        self.faXBP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.faXBP.setObjectName("faXBP")
        self.gridLayout_20.addWidget(self.faXBP, 2, 1, 1, 1)
        self.label_59 = QtWidgets.QLabel(self.page_9)
        self.label_59.setObjectName("label_59")
        self.gridLayout_20.addWidget(self.label_59, 3, 0, 1, 1)
        self.faYBP = QtWidgets.QPlainTextEdit(self.page_9)
        self.faYBP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.faYBP.setObjectName("faYBP")
        self.gridLayout_20.addWidget(self.faYBP, 3, 1, 1, 1)
        self.stackedWidget_4.addWidget(self.page_9)
        self.gridLayout_18.addWidget(self.stackedWidget_4, 1, 0, 1, 2)
        self.stackedWidget_2.addWidget(self.page_7)
        self.page_10 = QtWidgets.QWidget()
        self.page_10.setObjectName("page_10")
        self.formLayout = QtWidgets.QFormLayout(self.page_10)
        self.formLayout.setObjectName("formLayout")
        self.BSTypeComboBox = QtWidgets.QComboBox(self.page_10)
        self.BSTypeComboBox.setObjectName("BSTypeComboBox")
        self.BSTypeComboBox.addItem("")
        self.BSTypeComboBox.addItem("")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.SpanningRole, self.BSTypeComboBox)
        self.stackedWidget_5 = QtWidgets.QStackedWidget(self.page_10)
        self.stackedWidget_5.setObjectName("stackedWidget_5")
        self.page_11 = QtWidgets.QWidget()
        self.page_11.setObjectName("page_11")
        self.gridLayout_21 = QtWidgets.QGridLayout(self.page_11)
        self.gridLayout_21.setObjectName("gridLayout_21")
        self.label_67 = QtWidgets.QLabel(self.page_11)
        self.label_67.setObjectName("label_67")
        self.gridLayout_21.addWidget(self.label_67, 0, 0, 1, 1)
        self.f0BS = QtWidgets.QPlainTextEdit(self.page_11)
        self.f0BS.setMaximumSize(QtCore.QSize(16777215, 50))
        self.f0BS.setObjectName("f0BS")
        self.gridLayout_21.addWidget(self.f0BS, 0, 1, 1, 1)
        self.label_66 = QtWidgets.QLabel(self.page_11)
        self.label_66.setObjectName("label_66")
        self.gridLayout_21.addWidget(self.label_66, 1, 0, 1, 1)
        self.dFpBS = QtWidgets.QPlainTextEdit(self.page_11)
        self.dFpBS.setMaximumSize(QtCore.QSize(16777215, 50))
        self.dFpBS.setObjectName("dFpBS")
        self.gridLayout_21.addWidget(self.dFpBS, 1, 1, 1, 1)
        self.label_65 = QtWidgets.QLabel(self.page_11)
        self.label_65.setObjectName("label_65")
        self.gridLayout_21.addWidget(self.label_65, 2, 0, 1, 1)
        self.dFaBS = QtWidgets.QPlainTextEdit(self.page_11)
        self.dFaBS.setMaximumSize(QtCore.QSize(16777215, 50))
        self.dFaBS.setObjectName("dFaBS")
        self.gridLayout_21.addWidget(self.dFaBS, 2, 1, 1, 1)
        self.stackedWidget_5.addWidget(self.page_11)
        self.page_12 = QtWidgets.QWidget()
        self.page_12.setObjectName("page_12")
        self.gridLayout_22 = QtWidgets.QGridLayout(self.page_12)
        self.gridLayout_22.setObjectName("gridLayout_22")
        self.label_69 = QtWidgets.QLabel(self.page_12)
        self.label_69.setObjectName("label_69")
        self.gridLayout_22.addWidget(self.label_69, 0, 0, 1, 1)
        self.fpXBS = QtWidgets.QPlainTextEdit(self.page_12)
        self.fpXBS.setMaximumSize(QtCore.QSize(16777215, 50))
        self.fpXBS.setObjectName("fpXBS")
        self.gridLayout_22.addWidget(self.fpXBS, 0, 1, 1, 1)
        self.label_71 = QtWidgets.QLabel(self.page_12)
        self.label_71.setObjectName("label_71")
        self.gridLayout_22.addWidget(self.label_71, 1, 0, 1, 1)
        self.fpYBS = QtWidgets.QPlainTextEdit(self.page_12)
        self.fpYBS.setMaximumSize(QtCore.QSize(16777215, 50))
        self.fpYBS.setObjectName("fpYBS")
        self.gridLayout_22.addWidget(self.fpYBS, 1, 1, 1, 1)
        self.label_70 = QtWidgets.QLabel(self.page_12)
        self.label_70.setObjectName("label_70")
        self.gridLayout_22.addWidget(self.label_70, 2, 0, 1, 1)
        self.faXBS = QtWidgets.QPlainTextEdit(self.page_12)
        self.faXBS.setMaximumSize(QtCore.QSize(16777215, 50))
        self.faXBS.setObjectName("faXBS")
        self.gridLayout_22.addWidget(self.faXBS, 2, 1, 1, 1)
        self.label_68 = QtWidgets.QLabel(self.page_12)
        self.label_68.setObjectName("label_68")
        self.gridLayout_22.addWidget(self.label_68, 3, 0, 1, 1)
        self.faYBS = QtWidgets.QPlainTextEdit(self.page_12)
        self.faYBS.setMaximumSize(QtCore.QSize(16777215, 50))
        self.faYBS.setObjectName("faYBS")
        self.gridLayout_22.addWidget(self.faYBS, 3, 1, 1, 1)
        self.stackedWidget_5.addWidget(self.page_12)
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.SpanningRole, self.stackedWidget_5)
        self.stackedWidget_2.addWidget(self.page_10)
        self.gridLayout_2.addWidget(self.stackedWidget_2, 1, 0, 1, 2)
        self.label_73 = QtWidgets.QLabel(self.groupBox_2)
        self.label_73.setMaximumSize(QtCore.QSize(100, 50))
        self.label_73.setObjectName("label_73")
        self.gridLayout_2.addWidget(self.label_73, 2, 0, 2, 1)
        self.aaVal = QtWidgets.QPlainTextEdit(self.groupBox_2)
        self.aaVal.setMaximumSize(QtCore.QSize(16777215, 50))
        self.aaVal.setObjectName("aaVal")
        self.gridLayout_2.addWidget(self.aaVal, 4, 1, 3, 1)
        self.crearPlantillaButton = QtWidgets.QPushButton(self.groupBox_2)
        self.crearPlantillaButton.setObjectName("crearPlantillaButton")
        self.gridLayout_2.addWidget(self.crearPlantillaButton, 8, 0, 3, 2)
        self.label_72 = QtWidgets.QLabel(self.groupBox_2)
        self.label_72.setMaximumSize(QtCore.QSize(100, 50))
        self.label_72.setObjectName("label_72")
        self.gridLayout_2.addWidget(self.label_72, 4, 0, 3, 1)
        self.apVal = QtWidgets.QPlainTextEdit(self.groupBox_2)
        self.apVal.setMaximumSize(QtCore.QSize(16777215, 50))
        self.apVal.setObjectName("apVal")
        self.gridLayout_2.addWidget(self.apVal, 2, 1, 2, 1)
        self.listPlantillasWidget = QtWidgets.QListWidget(self.groupBox_2)
        self.listPlantillasWidget.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.listPlantillasWidget.setObjectName("listPlantillasWidget")
        self.gridLayout_2.addWidget(self.listPlantillasWidget, 11, 0, 1, 2)
        self.borrarPlantillasButton = QtWidgets.QPushButton(self.groupBox_2)
        self.borrarPlantillasButton.setObjectName("borrarPlantillasButton")
        self.gridLayout_2.addWidget(self.borrarPlantillasButton, 13, 0, 2, 1)
        self.borrarTodasPlantillasButton = QtWidgets.QPushButton(self.groupBox_2)
        self.borrarTodasPlantillasButton.setObjectName("borrarTodasPlantillasButton")
        self.gridLayout_2.addWidget(self.borrarTodasPlantillasButton, 13, 1, 2, 1)
        self.graficarPlantillasButton = QtWidgets.QPushButton(self.groupBox_2)
        self.graficarPlantillasButton.setObjectName("graficarPlantillasButton")
        self.gridLayout_2.addWidget(self.graficarPlantillasButton, 15, 0, 1, 2)
        self.gridLayout_5.addLayout(self.gridLayout_2, 0, 0, 1, 1)
        self.gridLayout_3.addWidget(self.groupBox_2, 0, 0, 10, 1)
        self.groupBox = QtWidgets.QGroupBox(self.tab_5)
        self.groupBox.setMaximumSize(QtCore.QSize(200, 16777215))
        self.groupBox.setTitle("")
        self.groupBox.setObjectName("groupBox")
        self.gridLayout_25 = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout_25.setObjectName("gridLayout_25")
        self.gridLayout_24 = QtWidgets.QGridLayout()
        self.gridLayout_24.setObjectName("gridLayout_24")
        self.borrarFiltrosButton = QtWidgets.QPushButton(self.groupBox)
        self.borrarFiltrosButton.setObjectName("borrarFiltrosButton")
        self.gridLayout_24.addWidget(self.borrarFiltrosButton, 5, 0, 1, 1)
        self.cambiarOrdenFiltroButton = QtWidgets.QPushButton(self.groupBox)
        self.cambiarOrdenFiltroButton.setObjectName("cambiarOrdenFiltroButton")
        self.gridLayout_24.addWidget(self.cambiarOrdenFiltroButton, 8, 0, 1, 1)
        self.listFiltrosWidget = QtWidgets.QListWidget(self.groupBox)
        self.listFiltrosWidget.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.listFiltrosWidget.setObjectName("listFiltrosWidget")
        self.gridLayout_24.addWidget(self.listFiltrosWidget, 4, 0, 1, 2)
        self.crearFiltroButton = QtWidgets.QPushButton(self.groupBox)
        self.crearFiltroButton.setObjectName("crearFiltroButton")
        self.gridLayout_24.addWidget(self.crearFiltroButton, 2, 0, 2, 2)
        self.rangoLP = QtWidgets.QPlainTextEdit(self.groupBox)
        self.rangoLP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.rangoLP.setObjectName("rangoLP")
        self.gridLayout_24.addWidget(self.rangoLP, 9, 1, 1, 1)
        self.ordenLP = QtWidgets.QTextBrowser(self.groupBox)
        self.ordenLP.setMaximumSize(QtCore.QSize(16777215, 50))
        self.ordenLP.setObjectName("ordenLP")
        self.gridLayout_24.addWidget(self.ordenLP, 8, 1, 1, 1)
        self.graficarFiltrosButton = QtWidgets.QPushButton(self.groupBox)
        self.graficarFiltrosButton.setObjectName("graficarFiltrosButton")
        self.gridLayout_24.addWidget(self.graficarFiltrosButton, 6, 0, 2, 2)
        self.tipoAproxComboBox = QtWidgets.QComboBox(self.groupBox)
        self.tipoAproxComboBox.setObjectName("tipoAproxComboBox")
        self.tipoAproxComboBox.addItem("")
        self.tipoAproxComboBox.addItem("")
        self.tipoAproxComboBox.addItem("")
        self.tipoAproxComboBox.addItem("")
        self.gridLayout_24.addWidget(self.tipoAproxComboBox, 0, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBox)
        self.label_6.setObjectName("label_6")
        self.gridLayout_24.addWidget(self.label_6, 0, 0, 1, 1)
        self.borrarTodosFiltrosButton = QtWidgets.QPushButton(self.groupBox)
        self.borrarTodosFiltrosButton.setObjectName("borrarTodosFiltrosButton")
        self.gridLayout_24.addWidget(self.borrarTodosFiltrosButton, 5, 1, 1, 1)
        self.cambiarRangoFiltroButton = QtWidgets.QPushButton(self.groupBox)
        self.cambiarRangoFiltroButton.setObjectName("cambiarRangoFiltroButton")
        self.gridLayout_24.addWidget(self.cambiarRangoFiltroButton, 9, 0, 2, 1)
        self.gridLayout_25.addLayout(self.gridLayout_24, 0, 0, 1, 1)
        self.gridLayout_3.addWidget(self.groupBox, 0, 1, 1, 1)
        self.tabWidget.addTab(self.tab_5, "")
        self.tab_6 = QtWidgets.QWidget()
        self.tab_6.setObjectName("tab_6")
        self.tabWidget.addTab(self.tab_6, "")
        self.gridLayout.addWidget(self.tabWidget, 0, 2, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 954, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(0)
        self.stackedWidget_2.setCurrentIndex(0)
        self.stackedWidget_4.setCurrentIndex(0)
        self.stackedWidget_5.setCurrentIndex(0)
        self.tipoFiltroComboBox.currentIndexChanged['int'].connect(self.stackedWidget_2.setCurrentIndex) # type: ignore
        self.BPTypeComboBox.currentIndexChanged['int'].connect(self.stackedWidget_4.setCurrentIndex) # type: ignore
        self.BSTypeComboBox.currentIndexChanged['int'].connect(self.stackedWidget_5.setCurrentIndex) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab), _translate("MainWindow", "Curvas de atenuación"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_2), _translate("MainWindow", "Curvas de ganancia"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_3), _translate("MainWindow", "Polos y ceros"))
        self.label.setText(_translate("MainWindow", "Tipo de Filtro"))
        self.tipoFiltroComboBox.setItemText(0, _translate("MainWindow", "Pasa Bajos"))
        self.tipoFiltroComboBox.setItemText(1, _translate("MainWindow", "Pasa Altos"))
        self.tipoFiltroComboBox.setItemText(2, _translate("MainWindow", "Pasa Banda"))
        self.tipoFiltroComboBox.setItemText(3, _translate("MainWindow", "Rechaza Banda"))
        self.label_53.setText(_translate("MainWindow", "Fp"))
        self.label_27.setText(_translate("MainWindow", "Fa"))
        self.label_55.setText(_translate("MainWindow", "Fp"))
        self.label_54.setText(_translate("MainWindow", "Fa"))
        self.BPTypeComboBox.setItemText(0, _translate("MainWindow", "Anchos de Banda"))
        self.BPTypeComboBox.setItemText(1, _translate("MainWindow", "Frecuencias"))
        self.label_63.setText(_translate("MainWindow", "Fo"))
        self.label_62.setText(_translate("MainWindow", "ΔFp"))
        self.label_61.setText(_translate("MainWindow", "ΔFa"))
        self.label_58.setText(_translate("MainWindow", "Fp-"))
        self.label_60.setText(_translate("MainWindow", "Fp+"))
        self.label_64.setText(_translate("MainWindow", "Fa-"))
        self.label_59.setText(_translate("MainWindow", "Fa+"))
        self.BSTypeComboBox.setItemText(0, _translate("MainWindow", "Anchos de Banda"))
        self.BSTypeComboBox.setItemText(1, _translate("MainWindow", "Frecuencias"))
        self.label_67.setText(_translate("MainWindow", "Fo"))
        self.label_66.setText(_translate("MainWindow", "ΔFp"))
        self.label_65.setText(_translate("MainWindow", "ΔFa"))
        self.label_69.setText(_translate("MainWindow", "Fp-"))
        self.label_71.setText(_translate("MainWindow", "Fp+"))
        self.label_70.setText(_translate("MainWindow", "Fa-"))
        self.label_68.setText(_translate("MainWindow", "Fa+"))
        self.label_73.setText(_translate("MainWindow", "Ap"))
        self.crearPlantillaButton.setText(_translate("MainWindow", "Crear Plantilla"))
        self.label_72.setText(_translate("MainWindow", "Aa"))
        self.borrarPlantillasButton.setText(_translate("MainWindow", "Borrar selec."))
        self.borrarTodasPlantillasButton.setText(_translate("MainWindow", "Borrar todas"))
        self.graficarPlantillasButton.setText(_translate("MainWindow", "Graficar plantillas selec."))
        self.borrarFiltrosButton.setText(_translate("MainWindow", "Borrar selec."))
        self.cambiarOrdenFiltroButton.setText(_translate("MainWindow", "Cambiar Orden"))
        self.crearFiltroButton.setText(_translate("MainWindow", "Crear Filtro"))
        self.graficarFiltrosButton.setText(_translate("MainWindow", "Graficar filtros"))
        self.tipoAproxComboBox.setItemText(0, _translate("MainWindow", "Butterworth"))
        self.tipoAproxComboBox.setItemText(1, _translate("MainWindow", "Chebyshev I"))
        self.tipoAproxComboBox.setItemText(2, _translate("MainWindow", "Chebyshev II (inverso)"))
        self.tipoAproxComboBox.setItemText(3, _translate("MainWindow", "Cauer (elíptico)"))
        self.label_6.setText(_translate("MainWindow", "Tipo de filtro"))
        self.borrarTodosFiltrosButton.setText(_translate("MainWindow", "Borrar todos"))
        self.cambiarRangoFiltroButton.setText(_translate("MainWindow", "Cambiar Rango"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5), _translate("MainWindow", "Aproximaciones"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_6), _translate("MainWindow", "División de etapas"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
