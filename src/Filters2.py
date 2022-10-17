import numpy as np
import sympy as sp
from sympy.abc import s,w
import scipy.signal as ss
from sympy.functions.special import polynomials as poly
import matplotlib.pyplot as plt

MIN = 1e-10

class FilterClass:
    def __init__(self):

        self.transferFunction = None    #transferencia final
        self.tfLP = None                #transferencia del pasabajos


class ButterworthApprox(FilterClass):
    def __init__(self):
        super().__init__()

        def getLPTransferFunction(self, Fp, Fa, Ap, Aa):
            Wp = Fp * 2 * np.pi()
            Wa = Fa * 2 * np.pi()
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b, a = ss.butter(N, Wn, 'lowpass', True)

            return ss.TransferFunction(b,a)

        def getHPTransferFunction(self, Fp, Fa, Ap, Aa):
            Wp = 1
            Wa = 2 * np.pi() * Fp / (Fa * 2 * np.pi())

            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', True)
            b, a = ss.lp2hp(b1,a1, Fp * 2 * np.pi())

            return ss.TransferFunction(b, a)

        def getBPTransferFunction(self, Fo, dFp, dFa, Ap, Aa):

            Fa = dFa / 2 + np.sqrt((dFa**2 + 4 * Fo**2)) / 2        # Hallo el valor de la frecuencia de atenuación más alta

            Wp = 1
            Wa = np.abs((-(2 * np.pi() * Fa) ** 2 + (Fo * 2 * np.pi()) ** 2) /
                        (- (dFp * 2 *np.pi()) * (2 * np.pi() * Fa)))

            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', True)
            b, a = ss.lp2bp(b1, a1, Fo * 2 * np.pi(), dFp * 2 *np.pi())

            return ss.TransferFunction(b, a)

        def getBPTransferFunction(self, Fps, Fas, Ap, Aa):

            Fo = np.sqrt(Fps[0]*Fps[1])
            Fp = Fps[1]

            dFp = Fps[1] - Fps[0]

            if (Fas[1] / Fo < Fo / Fas[0]):
                Fa = Fas[1]

            else:
                Fa = 2*Fo - Fas[0]

            Wp = 1
            Wa = np.abs((-(2 * np.pi() * Fa) ** 2 + (Fo * 2 * np.pi()) ** 2) /
                        (- (dFp * 2 * np.pi()) * (2 * np.pi() * Fa)))

            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', True)
            b, a = ss.lp2bp(b1, a1, Fo * 2 * np.pi(), dFp * 2 *np.pi())

            return ss.TransferFunction(b, a)

        def getBSTransferFunction(self, Fo, dFp, dFa, Ap, Aa):

            Fa = dFa / 2 + np.sqrt((dFa ** 2 + 4 * Fo ** 2)) / 2  # Hallo el valor de la frecuencia de atenuación más alta

            Wp = 1
            Wa = np.abs((-(2 * np.pi() * Fa) ** 2 + (Fo * 2 * np.pi()) ** 2) /
                        (- (dFp * 2 * np.pi()) * (2 * np.pi() * Fa)))

            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', True)
            b, a = ss.lp2hp(b1, a1, Fo * 2 * np.pi())

            return ss.TransferFunction(b, a, Fo * 2 * np.pi())

        def getBSTransferFunction(self, Fps, Fas, Ap, Aa):

            Fo = np.sqrt(Fps[0] * Fps[1])
            Fp = Fps[1]

            dFp = Fps[1] - Fps[0]

            if (Fas[1] / Fo < Fo / Fas[0]):
                Fa = Fas[1]

            else:
                Fa = 2 * Fo - Fas[0]

            Wp = 1
            Wa = np.abs((-(2 * np.pi() * Fa) ** 2 + (Fo * 2 * np.pi()) ** 2) /
                        (- (dFp * 2 * np.pi()) * (2 * np.pi() * Fa)))

            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', True)
            b, a = ss.lp2hp(b1, a1, Fo * 2 * np.pi())

            return ss.TransferFunction(b, a, Fp * 2 * np.pi())


class ChebyshevApproxI(FilterClass):
    def __init__(self, Ap, Aa, Fp, Fa):
        super().__init__(Ap, Aa)

        self.Fa = Fa
        self.Fp = Fp
        self.Wp = 1
        self.Wa = Fa / Fp

        N, Wn = ss.cheb1ord(self.Wp, self.Wa, Ap, Aa, True)
        self.tfLP = ss.cheby1(N, Ap, Wn, 'lowpass', True)

        wp, m, p = ss.bode(self.tfLP, np.linspace(0, 2, 100))

        plt.plot(wp, 10 ** (m / 20))
        plt.show()

class ChebyshevApproxII(FilterClass):
    def __init__(self, Ap, Aa, Fp, Fa):
        super().__init__(Ap, Aa)

        self.Fa = Fa
        self.Fp = Fp
        self.Wp = 1
        self.Wa = Fa / Fp

        N, Wn = ss.cheb2ord(self.Wp, self.Wa, Ap, Aa, True)
        self.tfLP = ss.cheby2(N, Wn, 'lowpass', True)

