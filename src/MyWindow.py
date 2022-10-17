# PyQt5 modules
from PyQt5.QtWidgets import QMainWindow

# Project modules
from src.ui.mainwindow import Ui_MainWindow
from src.Filters2 import ButterworthApprox,ChebyshevApproxI, ChebyshevApproxII
import src.Filters

class MainWindow(QMainWindow, Ui_MainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)


        #f1 = src.Filters.ChebyshevApproxI(1,20,1,1.1)
        f = ButterworthApprox(1, 20, 1, 1.1)


