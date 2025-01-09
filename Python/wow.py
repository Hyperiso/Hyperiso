from Python.Phyperiso import MemoryManager
from Python.Phyperiso import Parameters
from Python.Phyperiso import WilsonManager
from Python.Phyperiso import Model
from Python.Phyperiso import BCoefficientGroup
from Python.Phyperiso import ObservableInterface
from Python.Phyperiso import Observables

import os

mm = MemoryManager()
print("mmh")
print(os.getcwd() +"/Test/InputFiles/testInput.flha")
print(os.path.exists(os.getcwd() +"/Test/InputFiles/testInput.flha"))
mm.init("Test/InputFiles/testInput.flha")

pa = Parameters()
print("MASS TEST FOR TOP MASS", pa("MASS", 6))
print("ALPHA_S at 33 GeV", pa.alpha_s(33))
print("Runinng mass test", pa.running_mass(4.18, 4.18, 33))
print("check existence of value (bottom mass)", pa.exists("MASS", 5))
print("try getting mt_mt", pa.get_qcd_mass("mt_mt"))
pa.set_block_value("MASS", 3, 42)
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

mm.switch_lha("DataBase/lha/testInputv2.slha", Model.SM)
print(pa.exists("MASS", 5))
pa = Parameters()

print(pa("MASS",5))

mm.switch_lha("DataBase/lha/testInputv3.slha", Model.SM)
pa = Parameters()
print(pa("MASS", 5))

print(pa.exists("MASS", 5))

print(mm.get_blocks_list())
