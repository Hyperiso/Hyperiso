from pyhyperiso.phyperiso.pyhyperiso.observable import Estimate as _CppEstimate
from pyhyperiso.core.Math.scalar import Scalar

class Estimate:
    def __init__(self, cpp_obj):
        """Ne pas instancier manuellement. Utilisé en interne."""
        self._cpp_obj = cpp_obj

    @property
    def central_value(self) -> Scalar:
        return Scalar.from_cpp(self._cpp_obj.central_value)

    @property
    def stat_std(self) -> Scalar:
        return Scalar.from_cpp(self._cpp_obj.stat_std)

    @property
    def syst_std(self) -> Scalar:
        return Scalar.from_cpp(self._cpp_obj.syst_std)

    @property
    def combined_std(self) -> Scalar:
        return Scalar.from_cpp(self._cpp_obj.combined_std())

    def __repr__(self):
        return (f"Estimate(central = {self.central_value.value}, "
                f"stat = {self.stat_std.value}, "
                f"syst = {self.syst_std.value}, "
                f"combined = {self.combined_std.value})")

    @staticmethod
    def from_cpp(cpp_obj) -> "Estimate":
        return Estimate(cpp_obj)