
class PlantillaClass:
    def __init__(self):
        self.nombre = None
        self.tipoPlantilla = None
        self.ap = None
        self.aa = None
        self.fp = None
        self.fa = None
        self.f0 = None
        self.fpx = None
        self.fpy = None
        self.fax = None
        self.fay = None
        self.dfp = None
        self.dfa = None

    def crearPasaBajos(self, Fp, Fa, Ap, Aa):
        self.tipoPlantilla = "LP"
        self.fp = Fp
        self.fa = Fa
        self.ap = Ap
        self.aa = Aa

    def crearNombre(self, nombre):
        self.nombre = nombre

    def crearPasaAltos(self, Fp, Fa, Ap, Aa):
        self.tipoPlantilla = "HP"
        self.fp = Fp
        self.fa = Fa
        self.ap = Ap
        self.aa = Aa

    def crearPasaBandaFreq(self, Fo, Fpx, Fpy, Fax, Fay, Ap, Aa):
        self.tipoPlantilla = "BPF"
        self.fo = Fo
        self.fpx = Fpx
        self.fpy = Fpy
        self.fax = Fax
        self.fay = Fay
        self.ap = Ap
        self.aa = Aa

    def crearPasaBandaBW(self, dFp, dFa, Ap, Aa):
        self.tipoPlantilla = "BPBW"
        self.dfp = dFp
        self.dfa = dFa
        self.ap = Ap
        self.aa = Aa

    def crearRechazaBandaFreq(self, Fo, Fpx, Fpy, Fax, Fay, Ap, Aa):
        self.tipoPlantilla = "BSF"
        self.fo = Fo
        self.fpx = Fpx
        self.fpy = Fpy
        self.fax = Fax
        self.fay = Fay
        self.ap = Ap
        self.aa = Aa

    def crearRechazaBandaBW(self, dFp, dFa, Ap, Aa):
        self.tipoPlantilla = "BSBW"
        self.dfp = dFp
        self.dfa = dFa
        self.ap = Ap
        self.aa = Aa