# PyQt5 modules
from PyQt5.QtWidgets import QMainWindow

# Project modules
from src.ui.mainwindow import Ui_MainWindow
from src.Filters import ButterworthApprox,ChebyshevApproxI, ChebyshevApproxII


class MainWindow(QMainWindow, Ui_MainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)

        f = ChebyshevApproxII(1, 20, 1, 1.1)


