import numpy as np
import sympy as sp
from sympy.abc import s,w
import scipy.signal as ss
from sympy.functions.special import polynomials as poly

MIN = 1e-10

class FilterClass:
    def __init__(self, Ap, Aa):
        self.Ap = Ap            #ganancia de paso en dB (positiva)
        self.Aa = Aa            #ganancia de atenuación en dB (positiva)

        self.transferFunction = None    #transferencia final
        self.tfLP = None
        self.tfLP = None
        self.polinomial = None
        self.Zw = None
        self.Fs = None
        self.poles = None
        self.zeros= None

    def findPolesZeros(self):
        num = np.array(sp.Poly(sp.fraction(self.Fs)[0], s).all_coeffs(), dtype=float) #coef num
        den = np.array(sp.Poly(sp.fraction(self.Fs)[1], s).all_coeffs(), dtype=float) #coef den

        pz = ss.tf2zpk(num, den)    #polos y ceros

        usefulPoles = [p for p in pz[1] if np.real(p) < 0]
        useful_z = [p for p in pz[0] if np.real(p) < MIN]

        self.zeros = useful_z
        self.poles = usefulPoles

    def findLPTransferFunction(self):
        gain = np.cumprod(self.poles)[-1] / np.cumprod(self.zeros)[-1]      #se normaliza la función
        self.tfLP = ss.ZerosPolesGain(self.zeros, self.poles, gain)

class ButterworthApprox(FilterClass):
    def __init__(self, Ap, Aa, Fp, Fa):
        super().__init__(Ap, Aa)

        self.Fa = Fa
        self.Fp = Fp
        self.Wp = 1
        self.Wa = Fa / Fp

        self.xi = np.sqrt(10 ** (self.Ap / 10) - 1)  # se halla el xi

        self.n = int(np.ceil(np.log10( (10**(self.Aa/10)-1) / (self.xi**2)) /
                         (2 * np.log10(self.Wa) )))

        self.setTransferFunctionOrderN(self.n)

    def setTransferFunctionOrderN(self, n):
        self.polinomial = w ** n

        self.Zw = sp.Lambda(w, 1 / (1 + (self.polinomial * self.xi) ** 2))
        self.Fs = sp.simplify(self.Zw(s / sp.I))

        self.findPolesZeros()
        self.findLPTransferFunction()

class ChebyshevApproxI(FilterClass):
    def __init__(self, Ap, Aa, Fp, Fa):
        super().__init__(Ap, Aa)

        self.Fa = Fa
        self.Fp = Fp
        self.Wp = 1
        self.Wa = Fa / Fp

        self.xi = np.sqrt(10 ** (self.Ap / 10) - 1)  # se halla el xi igual que antes

        self.n = int(np.ceil( np.arccosh( (np.sqrt(10 ** (self.Aa / 10) - 1))  / (np.sqrt(10 ** (self.Ap / 10) - 1)))  /
                             (np.arccosh(self.Wa))))

        self.setTransferFunctionOrderN(self.n)

    def setTransferFunctionOrderN(self, n):
        self.polinomial = poly.chebyshevt(n, w)

        self.Zw = sp.Lambda(w, 1 / (1 + (self.polinomial * self.xi) ** 2))
        self.Fs = sp.simplify(self.Zw(s / sp.I))

        self.findPolesZeros()
        self.findLPTransferFunction()

class ChebyshevApproxII(FilterClass):
    def __init__(self, Ap, Aa, Fp, Fa):
        super().__init__(Ap, Aa)

        self.Fa = Fa
        self.Fp = Fp
        self.Wp = 1
        self.Wa = Fa / Fp

        self.xi = np.sqrt(10 ** (self.Ap / 10) - 1)  # se halla el xi igual que antes

        self.n = int(np.ceil( np.arccosh( (np.sqrt(10 ** (self.Aa / 10) - 1))  / (np.sqrt(10 ** (self.Ap / 10) - 1)))  /
                             (np.arccosh(self.Wa))))

        self.setTransferFunctionOrderN(self.n)

    def setTransferFunctionOrderN(self, n):
        self.polinomial = poly.chebyshevt(n, w)

        self.Zw = sp.Lambda(w, 1 / (1 + (self.polinomial * self.xi) ** 2))
        self.Fs = sp.simplify(self.Zw(s / sp.I))

        self.findPolesZeros()
        self.findLPTransferFunction()
