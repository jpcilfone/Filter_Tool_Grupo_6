# PyQt5 modules
from PyQt5.QtWidgets import QMainWindow

# Project modules
from src.ui.mainwindow import Ui_MainWindow
from src.Filters2 import FilterClass
import scipy.signal as ss
import matplotlib.pyplot as plt
import numpy as np

class MainWindow(QMainWindow, Ui_MainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)

        tipo = "HP"
        Fp = 10
        Fa = 3
        Ap = 5
        Aa = 80
        f = FilterClass()
        tF = f.getHPTransferFunction(Fp, Fa, Ap, Aa, "ellip", 1)


        if tipo == "LP":
            w = np.logspace(np.log10(Fp *2*np.pi / 10), np.log10(Fa *2*np.pi* 10), 1000)
            w, m, p = ss.bode(tF, w, n=1000)
            plt.semilogx(w / (2 * np.pi), m)

            plt.fill([Fp/10, Fp, Fp, Fp/10], [0, 0, -Ap, -Ap], '0.9', lw=0)
            plt.fill([Fa, Fa*10, Fa*10, Fa], [-Aa, -Aa, -Aa -10, -Aa -10], '0.9', lw=0)

        elif tipo == "HP":
            w = np.logspace(np.log10(Fa * 2 * np.pi / 10), np.log10(Fp * 2 * np.pi * 10), 1000)
            w, m, p = ss.bode(tF, w, n=1000)
            plt.semilogx(w / (2 * np.pi), m)

            plt.fill([Fp, Fp*10, Fp*10, Fp], [0, 0, -Ap, -Ap], '0.9', lw=0)
            plt.fill([Fa/10, Fa, Fa, Fa/10], [-Aa, -Aa, -Aa*2, -Aa*2], '0.9', lw=0)

        elif tipo == "BP":
            w = np.logspace(np.log10(Fa[0] * 2 * np.pi / 10), np.log10(Fa[1] * 2 * np.pi*2), 1000)
            w, m, p = ss.bode(tF, w, n=1000)
            plt.semilogx(w / (2 * np.pi), m)

            plt.fill([Fp[0], Fp[1], Fp[1], Fp[0]], [0, 0, -Ap, -Ap], '0.9', lw=0)
            plt.fill([Fa[0]/10, Fa[0], Fa[0], Fa[0]/10], [-Aa, -Aa, -Aa*2, -Aa*2], '0.9', lw=0)
            plt.fill([Fa[1], Fa[1]*2, Fa[1]*2, Fa[1]], [-Aa, -Aa, -Aa * 2, -Aa * 2], '0.9', lw=0)

        elif tipo == "BS":
            w = np.logspace(np.log10(Fa[0] * 2 * np.pi / 10), np.log10(Fa[1] * 2 * np.pi * 10), 1000)
            w, m, p = ss.bode(tF, w, n=1000)
            plt.semilogx(w / (2 * np.pi), m)

            plt.fill([Fa[0], Fa[1], Fa[1], Fa[0]], [-Aa, -Aa, -Aa-10, -Aa-10], '0.9', lw=0)
            plt.fill([Fp[0]/10, Fp[0], Fp[0], Fp[0]/10], [0, 0, -Ap, -Ap], '0.9', lw=0)
            plt.fill([Fp[1], Fp[1]*10, Fp[1]*10, Fp[1]], [0, 0, -Ap, -Ap], '0.9', lw=0)

        plt.show()

        N = f.n +1
        tF = f.changeFilterOrder(N)

        if tipo == "LP":
            w = np.logspace(np.log10(Fp * 2 * np.pi / 10), np.log10(Fa * 2 * np.pi * 10), 1000)
            w, m, p = ss.bode(tF, w, n=1000)
            plt.semilogx(w / (2 * np.pi), m)

            plt.fill([Fp / 10, Fp, Fp, Fp / 10], [0, 0, -Ap, -Ap], '0.9', lw=0)
            plt.fill([Fa, Fa * 10, Fa * 10, Fa], [-Aa, -Aa, -Aa - 10, -Aa - 10], '0.9', lw=0)

        elif tipo == "HP":
            w = np.logspace(np.log10(Fa * 2 * np.pi / 10), np.log10(Fp * 2 * np.pi * 10), 1000)
            w, m, p = ss.bode(tF, w, n=1000)
            plt.semilogx(w / (2 * np.pi), m)

            plt.fill([Fp, Fp * 10, Fp * 10, Fp], [0, 0, -Ap, -Ap], '0.9', lw=0)
            plt.fill([Fa / 10, Fa, Fa, Fa / 10], [-Aa, -Aa, -Aa * 2, -Aa * 2], '0.9', lw=0)

        elif tipo == "BP":
            w = np.logspace(np.log10(Fa[0] * 2 * np.pi / 10), np.log10(Fa[1] * 2 * np.pi * 2), 1000)
            w, m, p = ss.bode(tF, w, n=1000)
            plt.semilogx(w / (2 * np.pi), m)

            plt.fill([Fp[0], Fp[1], Fp[1], Fp[0]], [0, 0, -Ap, -Ap], '0.9', lw=0)
            plt.fill([Fa[0] / 10, Fa[0], Fa[0], Fa[0] / 10], [-Aa, -Aa, -Aa * 2, -Aa * 2], '0.9', lw=0)
            plt.fill([Fa[1], Fa[1] * 2, Fa[1] * 2, Fa[1]], [-Aa, -Aa, -Aa * 2, -Aa * 2], '0.9', lw=0)

        elif tipo == "BS":
            w = np.logspace(np.log10(Fa[0] * 2 * np.pi / 10), np.log10(Fa[1] * 2 * np.pi * 10), 1000)
            w, m, p = ss.bode(tF, w, n=1000)
            plt.semilogx(w / (2 * np.pi), m)

            plt.fill([Fa[0], Fa[1], Fa[1], Fa[0]], [-Aa, -Aa, -Aa - 10, -Aa - 10], '0.9', lw=0)
            plt.fill([Fp[0] / 10, Fp[0], Fp[0], Fp[0] / 10], [0, 0, -Ap, -Ap], '0.9', lw=0)
            plt.fill([Fp[1], Fp[1] * 10, Fp[1] * 10, Fp[1]], [0, 0, -Ap, -Ap], '0.9', lw=0)

        plt.show()










