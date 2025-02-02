from Python.Phyperiso import MemoryManager
from Python.Phyperiso import Parameters
from Python.Phyperiso import WilsonManager
from Python.Phyperiso import Model
from Python.Phyperiso import BCoefficientGroup, BPrimeCoefficientGroup
from Python.Phyperiso import ObservableInterface
from Python.Phyperiso import Observables
from Python.Phyperiso import ParameterType
import os

mm = MemoryManager()
print("mmh")
print(os.getcwd() +"/Test/InputFiles/testInput.flha")
print(os.path.exists(os.getcwd() +"/Test/InputFiles/testInput.flha"))
mm.init("Test/InputFiles/testInput.flha")

pa = Parameters()
print("MASS TEST FOR TOP MASS", pa("MASS", 6))
# print("ALPHA_S at 33 GeV", pa.alpha_s(33))
# print("Runinng mass test", pa.running_mass(4.18, 4.18, 33))
print("check existence of value (bottom mass)", pa.exists("MASS", 5))
# print("try getting mt_mt", pa.get_qcd_mass("mt_mt"))
# pa.set_block_value("MASS", 3, 42)
print("setting charm mass to 42, after changes: ", pa("MASS", 3))

wilManag = WilsonManager()

print("trying to set a group")

coeffgroup = BCoefficientGroup()

print("register a group in the manager")
wilManag.register_coefficient_group("BCoefficientGroup", coeffgroup)

print("setting Q_match")
wilManag.set_q_match("BCoefficientGroup", 81)

print("setting matching coefficients")
wilManag.set_matching_coefficient("BCoefficientGroup", "LO")

print("setting running value")
wilManag.set_group_scale("BCoefficientGroup", 4)

print("setting running coefficient")
wilManag.set_run_coefficient("BCoefficientGroup", "LO")

print("Matching C7", wilManag.get_matching_coefficient("BCoefficientGroup", "C7", "LO"))

print("Running C8", wilManag.get_run_coefficient("BCoefficientGroup", "C8", "LO"))

# obsInterface = ObservableInterface()

# print("obs truc", obsInterface.compute_observable(Observables.ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA))

# print("list of blocks : ", mm.get_blocks_list())

# print("info of block mass", mm.get_block_infos("MASS"))

# for block in mm.get_blocks_list():
#     print(block, mm.get_block_infos(block))
# print("info of block mass", mm.get_block_infos("MASS"))

print(pa("MASS", 5))

mm.switch_lha("Test/InputFiles/testinput_thdm.lha", Model.SM)
print(pa.exists("MASS", 5))
pa = Parameters()

print(pa("MASS",5))

mm.switch_lha("Test/InputFiles/testInput.slha", Model.SM)
pa = Parameters()
print(pa("MASS", 5))

print(pa.exists("MASS", 5))
print(mm.switch_model(Model.SUSY))
print(mm.get_blocks_list(ParameterType.SUSY))

coeffprimegroup = BPrimeCoefficientGroup()

print("register a prime group in the manager")
wilManag.register_coefficient_group("BPrimeCoefficientGroup", coeffprimegroup)
wilManag.set_q_match("BPrimeCoefficientGroup", 81)
wilManag.set_matching_coefficient("BPrimeCoefficientGroup", "LO")
wilManag.set_group_scale("BPrimeCoefficientGroup", 81)
wilManag.set_run_coefficient("BPrimeCoefficientGroup", "LO")
print(wilManag.get_matching_coefficient("BPrimeCoefficientGroup", "CP1", "LO"))
wilManag.set_q_match("BPrimeCoefficientGroup", 81)
wilManag.set_matching_coefficient("BPrimeCoefficientGroup", "LO")
print(wilManag.get_matching_coefficient("BPrimeCoefficientGroup", "CP1", "LO"))
# print(wilManag.get_state())
