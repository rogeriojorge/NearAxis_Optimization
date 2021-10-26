#!/usr/bin/env python3

from stel_repo import get_stel
from util import runqsc, runVMEC, runBOOZXFORM, runNEO, runSPEC, runREGCOIL
from simsopt_driver import optimize
from pathlib import Path
import os

results_folder = 'Results'
executables_folder = 'Executables'
plotting_folder = 'Plotting'

stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP = get_stel(24, nphi=131)
## Optimize for magnetic well d2_volume_d_psi2

#####
# Chaning weights of least_squares_problem
# example: LeastSquaresProblem.terms[0].weight = 4.0
#####

iota_target = 0.42
nIterations = 30
abs_step_array = [1e-1,1e-2,1e-3,1e-4,1e-6]
rel_step_array = [1e-1,1e-2]
# abs_step_array = [1e-2]
# rel_step_array = [1e-1]
max_fourier_coefficients = 6

Optimize = False
# Optimize = True

# Create folder for the results
Path(results_folder+'/'+name).mkdir(parents=True, exist_ok=True)
# Obtain folder paths'
main_path = str(Path(__file__).parent.resolve())
results_path = str(Path(results_folder+'/'+name).resolve())
executables_path = str(Path(executables_folder).resolve())
plotting_path = str(Path(plotting_folder).resolve())

# Go to results folder
os.chdir(results_path)

if Optimize:
    try:
        rel_step_array
        abs_step_array
        optimize(stel,iota_target,nIterations=nIterations,rel_step_array=rel_step_array,abs_step_array=abs_step_array,grad=True,max_fourier_coefficients=max_fourier_coefficients)
    except:
        optimize(stel,iota_target,nIterations,max_fourier_coefficients=max_fourier_coefficients)

## runqsc(stel,name,r_edge,executables_path,plotting_path) # DEPRECATED
# stel.plot(savefig='pyQSC_out.'+name+'.params')
# stel.plot_boundary(r=r_edge,fieldlines=True,savefig='pyQSC_out.'+name+'.boundary')
# stel.B_contour(r=0.06)

# stel.to_vmec('input.'+name,r=r_edge, ntheta=26, 
#         params={"ns_array": [16, 49, 101, 151, 201, 251],
#                 "ftol_array": [1e-17,1e-16,1e-15,1e-14,1e-14,1e-13],
#                 "niter_array": [2000,2000,2000,3000,4000,6000]})
# runVMEC(name,executables_path,plotting_path)

# runBOOZXFORM(name)
# runNEO(name,executables_path,plotting_path)
# runSPEC(name,executables_path,plotting_path,stel,r_edge)
# runREGCOIL(name,executables_path,plotting_path,coilSeparation = coilSeparation,targetValue = targetValue,nCoilsPerNFP = nCoilsPerNFP)

# Go back to main
os.chdir(main_path)