import numpy as np
import sympy as sp
from sympy.abc import s,w,n,xi,Z
import scipy.signal as ss

class FilterClass:
    def __init__(self, Ap, Aa):
        self.Ap = Ap            #ganancia de paso en dB (positiva)
        self.Aa = Aa            #ganancia de atenuación en dB (positiva)

        self.transferFunction = None    #transferencia final
        self.polinomial = None
        self.Zw = None
        self.Fs = None
        self.poles = None
        self.zeros= None

    def findPolesZeros(self):
        num = np.array(sp.Poly(sp.fraction(self.Fs)[0]).all_coeffs(), dtype=float) #coef num
        den = np.array(sp.Poly(sp.fraction(self.Fs)[1]).all_coeffs(), dtype=float) #coef den

        pz = ss.tf2zpk(num, den)    #polos y ceros

        usefulPoles = [p for p in pz[1] if np.real(p) < 0]

        self.zeros = pz[0]
        self.poles = usefulPoles

    def findTransferFunction(self):
        self.transferFunction = ss.ZerosPolesGain(self.zeroes, self.poles, 1)

class ButterworthApprox(FilterClass):
    def __init__(self, Ap, Aa, Fp, Fa):
        super().__init__(Ap, Aa)

        self.Fa = Fa
        self.Fp = Fp
        self.Wp = 1
        self.Wa = Fa / Fp

        self.xi = np.sqrt(10 ** (self.Ap / 10) - 1)  # se halla el xi

        self.n = np.ceil(np.log10( (10**(self.Aa/10)-1) / (self.xi**2)) /
                         (2 * np.log10(self.Wa) ))

        self.polinomial = w ** self.n

        self.Zw = sp.Lambda(w, 1 / (1 + (self.polinomial * self.xi)**2))
        self.Fs = self.Zw(s/sp.I)

        self.findPolesZeros()
        self.findTransferFunction()

        print(self.Fs)

class ChebyshevApproxI(FilterClass):
    def __init__(self, Ap, Aa, Fp, Fa):
        super().__init__(Ap, Aa)

        self.Fa = Fa
        self.Fp = Fp
        self.Wp = 1
        self.Wa = Fa / Fp

        self.xi = np.sqrt(10 ** (self.Ap / 10) - 1)  # se halla el xi

        self.n = np.ceil(np.log10( (10**(self.Aa/10)-1) / (self.xi**2)) /
                         (2 * np.log10(self.Wa) ))

        self.polinomial = w ** self.n

        self.Zw = sp.Lambda(w, 1 / (1 + (self.polinomial * self.xi)**2))
        self.Fs = sp.simplify(self.Zw(s/1j))

        print(self.Fs)
