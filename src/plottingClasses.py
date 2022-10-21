from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure()
        self.axes = self.fig.add_subplot()
        self.fig.set_tight_layout(True)
        super().__init__(self.fig)
        self.navToolBar = NavigationToolbar(self, parent=parent)

class BodePlot(MplCanvas):
    def __init__(self, parent=None):
        if parent is not None:
            super().__init__(parent)

    def plotMag(self, w, mag, labels):
        for i in range(len(w)):
            self.axes.semilogx(w[i], mag[i], label = labels[i])

        self.axes.set_xlabel('Frecuencia [Hz]')
        self.axes.set_ylabel('|H| [dB]')

        self.fig.canvas.draw()

    def plotFase(self, w, fase, labels):
        for i in range(len(w)):
            self.axes.semilogx(w[i], fase[i], label = labels[i])

        self.axes.set_xlabel('Frecuencia [Hz]')
        self.axes.set_ylabel('Fase (°)')

        self.fig.canvas.draw()


class PolosCerosPlot(MplCanvas):
    def __init__(self, parent=None):
        if parent is not None:
            super().__init__(parent)

    def plotPolosCeros(self, polesZeros, labels):
        for i in range(len(polesZeros)):
            poles = polesZeros[i][0]
            zeros = polesZeros[i][1]

            if len(poles) != 0:
                self.axesPolesZeros.scatter(np.real(poles), np.imag(poles),
                                                        marker='x', label='Polos ' + labels [i])

            if len(zeros) != 0:
                self.axesPolesZeros.scatter(np.real(zeros), np.imag(zeros),
                                                             marker='o', label='Ceros ' + labels[i])

        self.axes.ticklabel_format(axis='both', style='sci')
        self.axes.grid(visible=True)
        self.axes.axhline(0, color='black', linewidth=1)
        self.axes.axvline(0, color='black', linewidth=1)  # marcamos los ejes

        self.axes.set_ylabel('jω')
        self.axes.set_xlabel('σ')
        self.axes.legend(loc=0)

        self.fig.canvas.draw()

