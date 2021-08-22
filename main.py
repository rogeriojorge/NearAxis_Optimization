#!/usr/bin/env python3

from stel_repo import get_stel
from util import runqsc, runVMEC, runBOOZXFORM, runNEO, runSPEC, runREGCOIL
from simsopt_driver import optimize
from pathlib import Path
import os

results_folder = 'Results'
executables_folder = 'Executables'
plotting_folder = 'Plotting'

stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP = get_stel(1, nphi=251)

iota_target = 0.41
nIterations = 20
abs_step_array = [1e-1,1e-2,1e-3,1e-4,1e-6]
rel_step_array = [1e-1,1e-2,1e-3,1e-4]
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

# optimize(stel,iota_target,rel_step_array,abs_step_array,nIterations,grad=True)
# stel.plot(r=r_edge,fieldlines=True)
runqsc(stel,name,r_edge,executables_path,plotting_path)
runVMEC(name,executables_path,plotting_path)
runBOOZXFORM(name,executables_path,plotting_path)
runNEO(name,executables_path,plotting_path)
runREGCOIL(name,executables_path,plotting_path,coilSeparation = coilSeparation,targetValue = targetValue,nCoilsPerNFP = nCoilsPerNFP)
# runSPEC(name,executables_path,plotting_path)

# Go back to main
os.chdir(main_path)