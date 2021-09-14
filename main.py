#!/usr/bin/env python3

from stel_repo import get_stel
from util import runqsc, runVMEC, runBOOZXFORM, runNEO, runSPEC, runREGCOIL
from simsopt_driver import optimize
from pathlib import Path
import os

results_folder = 'Results'
executables_folder = 'Executables'
plotting_folder = 'Plotting'

stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP = get_stel(23, nphi=251)


#####
# Chaning weights of least_squares_problem
# example: LeastSquaresProblem.terms[0].weight = 4.0
#####

iota_target = 0.42
nIterations = 50
abs_step_array = [1e-1,1e-2,1e-3,1e-4,1e-6]
rel_step_array = [1e-1,1e-2,1e-3]
# abs_step_array = [1e-2]
# rel_step_array = [1e-1]
max_fourier_coefficients = 6

# Create folder for the results
Path(results_folder+'/'+name).mkdir(parents=True, exist_ok=True)
# Obtain folder paths'
main_path = str(Path(__file__).parent.resolve())
results_path = str(Path(results_folder+'/'+name).resolve())
executables_path = str(Path(executables_folder).resolve())
plotting_path = str(Path(plotting_folder).resolve())

# Go to results folder
os.chdir(results_path)

# optimize(stel,iota_target,rel_step_array,abs_step_array,nIterations,grad=False,max_fourier_coefficients=max_fourier_coefficients)
# stel.plot(r=r_edge,fieldlines=True)
# stel.B_contour(r=0.05)
# stel.plot_axis()
# runqsc(stel,name,r_edge,executables_path,plotting_path)
# stel.to_vmec('input.'+name,r=r_edge)

# runVMEC(name,executables_path,plotting_path)
# runBOOZXFORM(name)
# runNEO(name,executables_path,plotting_path)
# runREGCOIL(name,executables_path,plotting_path,coilSeparation = coilSeparation,targetValue = targetValue,nCoilsPerNFP = nCoilsPerNFP)
runSPEC(name,executables_path,plotting_path,stel,r_edge)

# Go back to main
os.chdir(main_path)