#!/usr/bin/env python3

from stel_repo import get_stel
from util import runqsc, runVMEC, runBOOZXFORM, runNEO, runSPEC, runREGCOIL
from simsopt_driver import optimize
from pathlib import Path
import os
import input

# If not optimizing, use a fine resolution
nphi_refined = max(input.nphi, 251)
try:
    if input.Optimize:
        nphi = input.nphi
    else:
        nphi = nphi_refined
except:
    nphi = nphi_refined

# Get stellarator from the repository
stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP = get_stel(input.ind, nphi=nphi)

## Folders operations
# Set the name of important folders
try:
    input.results_folder
except:
    results_folder = 'Results'
try:
    input.executables_folder
except:
    executables_folder = 'Executables'
try:
    input.plotting_folder
except:
    plotting_folder = 'Plotting'
# Create folder for the results
Path(results_folder+'/'+name).mkdir(parents=True, exist_ok=True)
# Obtain folder paths'
main_path = str(Path(__file__).parent.resolve())
results_path = str(Path(results_folder+'/'+name).resolve())
executables_path = str(Path(executables_folder).resolve())
plotting_path = str(Path(plotting_folder).resolve())

# Go to results folder
os.chdir(results_path)

# Run Optimization
try:
    if input.Optimize:
        # if rel_step_array is specified, perform gradient based optimization
        try:
            input.rel_step_array
            input.abs_step_array
            optimize(stel,input.iota_target,nIterations=input.nIterations,rel_step_array=input.rel_step_array,abs_step_array=input.abs_step_array,grad=True,max_fourier_coefficients=input.max_fourier_coefficients)
        except:
            optimize(stel,input.iota_target,nIterations=input.nIterations,max_fourier_coefficients=input.max_fourier_coefficients)
except:
    input.Optimize = False

# Do the plotting
try:
    if input.Plot:
        # runqsc(stel,name,r_edge,executables_path,plotting_path) # DEPRECATED
        stel.plot(savefig='pyQSC_out.'+name+'.params', show=False)
        stel.B_contour(r=r_edge, savefig='pyQSC_out.'+name, ncontours=25, show=False)
        stel.plot_axis(savefig='pyQSC_out.'+name+'axis', show=False)
        stel.B_fieldline(r=r_edge, savefig='pyQSC_out.'+name, show=False)
        stel.plot_boundary(r=r_edge, fieldlines=True, savefig='pyQSC_out.'+name+'.boundary', show=False)
except:
    input.Plot = False

# Run VMEC
try:
    if input.VMEC:
        stel.to_vmec('input.'+name,r=r_edge,
                params={"ns_array": [16, 49, 101, 151],
                        "ftol_array": [1e-17,1e-16,1e-15,1e-14],
                        "niter_array": [2000,2000,2000,3000]})
        runVMEC(name,executables_path,plotting_path)
except:
    input.VMEC = False

# Run BOOZ_XFORM
try:
    if input.BOOZ_XFORM:
        runBOOZXFORM(name)
except:
    input.BOOZ_XFORM = False

# Run NEO
try:
    if input.NEO:
        runNEO(name,executables_path,plotting_path)
except:
    input.NEO = False

# Run SPEC
try:
    if input.SPEC:
        runSPEC(name,executables_path,plotting_path,stel,r_edge)
except:
    input.SPEC = False

# Run REGCOIL
try:
    if input.REGCOIL:
        runREGCOIL(name,executables_path,plotting_path,coilSeparation = coilSeparation,targetValue = targetValue,nCoilsPerNFP = nCoilsPerNFP)
except:
    input.REGCOIL = False

# Go back to main
os.chdir(main_path)