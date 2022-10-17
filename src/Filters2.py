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

    '''

    '''
    def getLPTransferFunction(self, Fp, Fa, Ap, Aa, btype):
        Wp = Fp * 2 * np.pi
        Wa = Fa * 2 * np.pi

        a = [1]
        b = [1]

        if (btype == "butter"):
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b, a = ss.butter(N, Wn, 'lowpass', analog = True)

        elif (btype == "cheby1"):
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b, a = ss.cheby1(N, Ap, Wn, 'lowpass', analog = True)

        elif (btype == "cheby2"):
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b, a = ss.cheby2(N, Aa, Wn, 'lowpass', analog = True)

        elif (btype == "ellip"):
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b, a = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)

        self.transferFunction = ss.TransferFunction(b,a)

        return self.transferFunction

    '''

    '''
    def getHPTransferFunction(self, Fp, Fa, Ap, Aa, btype):
        Wp = 1
        Wa = 2 * np.pi * Fp / (Fa * 2 * np.pi)

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)


        b, a = ss.lp2hp(b1,a1, Fp * 2 * np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

    '''

    '''
    def getBPTransferFunctionBW(self, Fo, dFp, dFa, Ap, Aa, btype):

        Fa = dFa / 2 + np.sqrt((dFa**2 + 4 * Fo**2)) / 2        # Hallo el valor de la frecuencia de atenuación más alta

        Wp = 1
        Wa = np.abs((-(2 * np.pi * Fa) ** 2 + (Fo * 2 * np.pi) ** 2) /
                    (- (dFp * 2 *np.pi) * (2 * np.pi * Fa)))

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)

        b, a = ss.lp2bp(b1, a1, Fo * 2 * np.pi, dFp * 2 *np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

    '''
    
    '''
    def getBPTransferFunctionFreq(self, Fps, Fas, Ap, Aa, btype):

        Fo = np.sqrt(Fps[0]*Fps[1])
        Fp = Fps[1]

        dFp = Fps[1] - Fps[0]

        if (Fas[1] / Fo < Fo / Fas[0]):
            Fa = Fas[1]

        else:
            Fa = 2*Fo - Fas[0]

        Wp = 1
        Wa = np.abs((-(2 * np.pi * Fa) ** 2 + (Fo * 2 * np.pi) ** 2) /
                    (- (dFp * 2 * np.pi) * (2 * np.pi * Fa)))

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)


        b, a = ss.lp2bp(b1, a1, Fo * 2 * np.pi, dFp * 2 *np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

    '''
    
    '''
    def getBSTransferFunctionBW(self, Fo, dFp, dFa, Ap, Aa, btype):

        Fa = dFa / 2 + np.sqrt((dFa ** 2 + 4 * Fo ** 2)) / 2  # Hallo el valor de la frecuencia de atenuación más alta

        Wp = 1
        Wa = np.abs(((dFp * 2 * np.pi) * (2 * np.pi * Fa)) /
            (-(2 * np.pi * Fa) ** 2 + (Fo * 2 * np.pi) ** 2))

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)

        b, a = ss.lp2bs(b1, a1, Fo * 2 * np.pi, dFp * 2 * np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

    '''

    '''
    def getBSTransferFunctionFreq(self, Fps, Fas, Ap, Aa, btype):

        Fo = np.sqrt(Fps[0] * Fps[1])
        Fp = Fps[1]

        dFp = Fps[1] - Fps[0]

        if (Fas[1] / Fo < Fo / Fas[0]):
            Fa = Fas[1]

        else:
            Fa = 2 * Fo - Fas[0]

        Wp = 1
        Wa = np.abs(((dFp * 2 * np.pi) * (2 * np.pi * Fa)) /
                    (-(2 * np.pi * Fa) ** 2 + (Fo * 2 * np.pi) ** 2))

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)

        b, a = ss.lp2bs(b1, a1, Fo * 2 * np.pi, dFp * 2 * np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

