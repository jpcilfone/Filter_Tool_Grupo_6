import numpy as np
import scipy.signal as ss
import sympy as sp
from sympy.abc import s,w
from sympy.functions.special import polynomials as poly

Epsilon = 1e-10

class FilterClass:
    def __init__(self):

        self.transferFunction = None    #transferencia final
        self.tfLP = None                #transferencia del pasabajos

        self.filterType = None          #Tipo de filtro
        self.filterName = None
        self.n = None
        self.Wp = None                  #Frecuencias del pasa bajos normalizado
        self.Wa = None
        self.Ap = None
        self.Aa = None

        self.Fp = None
        self.Fo = None
        self.dFp = None

        self.currentN = None            #nuevo N
        self.Wan = None              # nuevo Wan

    '''
        Low Pass
    '''
    def getLPTransferFunction(self, Fp, Fa, Ap, Aa, btype):

        self.filterType = "LP"

        Wp = 1
        Wa = Fa / Fp

        self.Wp = Wp
        self.Wa = Wa
        self.Ap = Ap
        self.Aa = Aa
        self.Fp = Fp

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            self.filterName = "butter"
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog = True)

        elif (btype == "cheby1"):
            self.filterName = "cheby1"
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog = True)

        elif (btype == "cheby2"):
            self.filterName = "cheby2"
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog = True)

        elif (btype == "ellip"):
            self.filterName = "ellip"
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)

        b, a = ss.lp2lp(b1, a1, Fp * 2 * np.pi)
        self.transferFunction = ss.TransferFunction(b,a)

        return self.transferFunction

    '''
        High Pass
    '''
    def getHPTransferFunction(self, Fp, Fa, Ap, Aa, btype):
        self.filterType = "HP"
        Wp = 1
        Wa = 2 * np.pi * Fp / (Fa * 2 * np.pi)

        self.Wp = Wp
        self.Wa = Wa
        self.Ap = Ap
        self.Aa = Aa
        self.Fp = Fp

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            self.filterName = "butter"
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            self.filterName = "cheby1"
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            self.filterName = "cheby2"
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            self.filterName = "ellip"
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)


        b, a = ss.lp2hp(b1,a1, Fp * 2 * np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

    '''
        Band Pass con ANCHOS DE BANDA
    '''
    def getBPTransferFunctionBW(self, Fo, dFp, dFa, Ap, Aa, btype):
        self.filterType = "BP"

        #Fa = dFa / 2 + np.sqrt((dFa**2 + 4 * Fo**2)) / 2        # Hallo el valor de la frecuencia de atenuación más alta

        Wp = 1
        Wa = dFa / dFp

        self.Wp = Wp
        self.Wa = Wa
        self.Ap = Ap
        self.Aa = Aa

        self.Fo = Fo
        self.dFp = dFp

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            self.filterName = "butter"
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            self.filterName = "cheby1"
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            self.filterName = "cheby2"
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            self.filterName = "ellip"
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)

        b, a = ss.lp2bp(b1, a1, Fo * 2 * np.pi, dFp * 2 *np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

    '''
        Band Pass con FRECUENCIAS CRÍTICAS
    '''
    def getBPTransferFunctionFreq(self, Fps, Fas, Ap, Aa, btype):
        self.filterType = "BP"

        Fo = np.sqrt(Fps[0]*Fps[1])
        Fp = Fps[1]

        dFp = Fps[1] - Fps[0]
        self.dFp = dFp

        Fa = 1

        if (Fas[1] / Fo < Fo / Fas[0]):
            Fa = Fas[1]
            dFa = Fa - Fo ** 2 / Fa

        else:
            Fa = Fas[0]
            dFa = Fo ** 2 / Fa - Fa

        Wp = 1
        Wa = dFa / dFp

        self.Wp = Wp
        self.Wa = Wa
        self.Ap = Ap
        self.Aa = Aa

        self.Fo = Fo
        self.dFp = dFp

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            self.filterName = "butter"
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            self.filterName = "cheby1"
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            self.filterName = "cheby2"
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            self.filterName = "ellip"
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)


        b, a = ss.lp2bp(b1, a1, Fo * 2 * np.pi, dFp * 2 * np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

    '''
        Band Stop con ANCHOS DE BANDA
    '''
    def getBSTransferFunctionBW(self, Fo, dFp, dFa, Ap, Aa, btype):
        self.filterType = "BS"
        #Fa = dFa / 2 + np.sqrt((dFa ** 2 + 4 * Fo ** 2)) / 2  # Hallo el valor de la frecuencia de atenuación más alta

        Wp = 1
        Wa = dFp / dFa

        self.Wp = Wp
        self.Wa = Wa
        self.Ap = Ap
        self.Aa = Aa

        self.Fo = Fo
        self.dFp = dFp

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            self.filterName = "butter"
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            self.filterName = "cheby1"
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            self.filterName = "cheby2"
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            self.filterName = "ellip"
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)

        b, a = ss.lp2bs(b1, a1, Fo * 2 * np.pi, dFp * 2 * np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

    '''
        Band Stop con FRECUENCIAS CRÍTICAS
    '''
    def getBSTransferFunctionFreq(self, Fps, Fas, Ap, Aa, btype):
        self.filterType = "BS"

        Fo = np.sqrt(Fps[0] * Fps[1])
        Fp = Fps[1]

        dFp = Fps[1] - Fps[0]

        if (Fas[1] / Fo > Fo / Fas[0]):
            Fa = Fas[1]
            dFa = Fa - Fo ** 2 / Fa

        else:
            Fa = Fas[0]
            dFa = Fo ** 2 / Fa - Fa

        Wp = 1
        Wa = dFp / dFa

        self.Wp = Wp
        self.Wa = Wa
        self.Ap = Ap
        self.Aa = Aa

        self.Fo = Fo
        self.dFp = dFp

        b1 = [1]
        a1 = [1]

        if (btype == "butter"):
            self.filterName = "butter"
            N, Wn = ss.buttord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.butter(N, Wn, 'lowpass', analog=True)

        elif (btype == "cheby1"):
            self.filterName = "cheby1"
            N, Wn = ss.cheb1ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby1(N, Ap, Wn, 'lowpass', analog=True)

        elif (btype == "cheby2"):
            self.filterName = "cheby2"
            N, Wn = ss.cheb2ord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.cheby2(N, Aa, Wn, 'lowpass', analog=True)

        elif (btype == "ellip"):
            self.filterName = "ellip"
            N, Wn = ss.ellipord(Wp, Wa, Ap, Aa, True)
            self.n = N
            b1, a1 = ss.ellip(N, Ap, Aa, Wn, 'lowpass', analog=True)

        b, a = ss.lp2bs(b1, a1, Fo * 2 * np.pi, dFp * 2 * np.pi)

        self.transferFunction = ss.TransferFunction(b, a)

        return self.transferFunction

    '''
        Cambiar el orden del filtro
    '''
    def changeFilterOrder(self, N):

        b1 = [1]
        a1 = [1]

        if self.filterName == "butter":
            Wan1 = 10 ** (np.log10((10 ** (self.Aa / 10) - 1) / (10 ** (self.Ap / 10) - 1)) / (2 * N))
            Wan2 = 10 ** (np.log10((10 ** (self.Aa / 10) - 1) / (10 ** (self.Ap / 10) - 1)) / (2 * (N-1)))

            Wan = (Wan1 + Wan2) / 2

            self.Wan = Wan

            N1, Wn = ss.buttord(self.Wp, Wan, self.Ap, self.Aa, True)
            self.currentN = N1                                      #N1 debería ser igual a N
            b1, a1 = ss.butter(N1, Wn, 'lowpass', analog=True)

        elif self.filterName == "cheby1":
            Wan1 = np.cosh(np.arccosh(np.sqrt((10 ** (self.Aa / 10) - 1) / (10 ** (self.Ap / 10) - 1))) / (N))
            Wan2 = np.cosh(np.arccosh(np.sqrt((10 ** (self.Aa / 10) - 1) / (10 ** (self.Ap / 10) - 1))) / (N-1))

            Wan = (Wan1 + Wan2) / 2

            self.Wan = Wan

            N1, Wn = ss.cheb1ord(self.Wp, Wan, self.Ap, self.Aa, True)
            self.currentN = N1  # N1 debería ser igual a N
            b1, a1 = ss.cheby1(N1, self.Ap, Wn, 'lowpass', analog=True)

        elif self.filterName == "cheby2":
            Wan1 = np.cosh(np.arccosh(np.sqrt((10 ** (self.Aa / 10) - 1) / (10 ** (self.Ap / 10) - 1))) / (N))
            Wan2 = np.cosh(np.arccosh(np.sqrt((10 ** (self.Aa / 10) - 1) / (10 ** (self.Ap / 10) - 1))) / (N - 1))

            Wan = (Wan1 + Wan2) / 2

            self.Wan = Wan

            N1, Wn = ss.cheb2ord(self.Wp, Wan, self.Ap, self.Aa, True)
            self.currentN = N1  # N1 debería ser igual a N
            b1, a1 = ss.cheby2(N1, self.Aa, Wn, 'lowpass', analog=True)

        elif self.filterName == "ellip":

            self.currentN = N  # N1 debería ser igual a N
            b1, a1 = ss.ellip(N, self.Ap, self.Aa, self.Wp, 'lowpass', analog=True) #debería funcionar


        if self.filterType == "LP":
            b, a = ss.lp2lp(b1, a1, self.Fp * 2 * np.pi)
            self.transferFunction = ss.TransferFunction(b, a)

        elif self.filterType == "HP":
            b, a = ss.lp2hp(b1, a1, self.Fp * 2 * np.pi)

            self.transferFunction = ss.TransferFunction(b, a)

        elif self.filterType == "BP":
            b, a = ss.lp2bp(b1, a1, self.Fo * 2 * np.pi, self.dFp * 2 * np.pi)

            self.transferFunction = ss.TransferFunction(b, a)

        elif self.filterType == "BS":
            b, a = ss.lp2bs(b1, a1, self.Fo * 2 * np.pi, self.dFp * 2 * np.pi)

            self.transferFunction = ss.TransferFunction(b, a)


        return self.transferFunction









