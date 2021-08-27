#!/usr/bin/env python3

from stel_repo import get_stel
from util import runqsc, runVMEC, runBOOZXFORM, runNEO, runSPEC, runREGCOIL
from simsopt_driver import optimize
from pathlib import Path
import os

results_folder = 'Results'
executables_folder = 'Executables'
plotting_folder = 'Plotting'

stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP = get_stel(7, nphi=151)

iota_target = 0.4
nIterations = 30
abs_step_array = [1e-2,1e-3,1e-4,1e-5,1e-6]
rel_step_array = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6]
# abs_step_array = [1e-2]
# rel_step_array = [1e-1]

# Create folder for the results
Path(results_folder+'/'+name).mkdir(parents=True, exist_ok=True)
# Obtain folder paths'
main_path = str(Path(__file__).parent.resolve())
results_path = str(Path(results_folder+'/'+name).resolve())
executables_path = str(Path(executables_folder).resolve())
plotting_path = str(Path(plotting_folder).resolve())

# Go to results folder
os.chdir(results_path)


##### LOOK AT EQS 3.12 and 3.14 from higher order paper
##### and calculate X3, Y3 and optimize for those

# optimize(stel,iota_target,rel_step_array,abs_step_array,nIterations,grad=True)
# stel.plot(r=r_edge,fieldlines=True)
# stel.B_contour(r=0.05)
# stel.plot_axis()
# runqsc(stel,name,r_edge,executables_path,plotting_path)

# ## USE WRAPPERS - VMEC, BOOZ_XFORM, SPEC
# runVMEC(name,executables_path,plotting_path)


##### Look at C++ version of Booz_Xform in Hidden Symmetries rpo
##### and get plotting routine from there to show both
##### QH and QA
####### DEVELOP BRANCH, symplot(b, helical_detail=True)
## https://github.com/hiddenSymmetries/booz_xform/blob/develop/src/booz_xform/plots.py#L97
# runBOOZXFORM(name,executables_path,plotting_path)
# runNEO(name,executables_path,plotting_path)
runREGCOIL(name,executables_path,plotting_path,coilSeparation = coilSeparation,targetValue = targetValue,nCoilsPerNFP = nCoilsPerNFP)
# runSPEC(name,executables_path,plotting_path)

# Go back to main
os.chdir(main_path)