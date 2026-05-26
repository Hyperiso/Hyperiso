from pyhyperiso.phyperiso.pyhyperiso.core import QCDConstants as _CppQCDConstants

class QCDConstants:
    def __init__(self, cpp_obj: _CppQCDConstants):
        self._cpp_obj = cpp_obj

    @property
    def Nc(self):
        return self._cpp_obj.Nc

    @property
    def C_F(self):
        return self._cpp_obj.C_F

    @property
    def C_A(self):
        return self._cpp_obj.C_A

    @property
    def beta(self):
        return self._cpp_obj.beta

    @property
    def gamma(self):
        return self._cpp_obj.gamma
    
    def __repr__(self):
        return (f"QCDConstants(Nc={self.Nc}, C_F={self.C_F}, C_A={self.C_A}, "
                f"beta={self.beta}, gamma={self.gamma})")