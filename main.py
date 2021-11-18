#!/usr/bin/env python3

from stel_repo import get_stel
from util import runqsc, runVMEC, runBOOZXFORM, runNEO, runSPEC, runREGCOIL
from simsopt_driver import optimize
from pathlib import Path
import os
import Input

print('Starting Near-Axis Optimization')

# If not optimizing, use a fine resolution
nphi_refined = max(Input.nphi, 301)
try:
    if Input.Optimize:
        nphi = Input.nphi
    else:
        nphi = nphi_refined
except Exception as e:
    nphi = nphi_refined

# Get stellarator from the repository
stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP = get_stel(Input.ind, nphi=nphi)

## Folders operations
# Set the name of important folders
try:
    Input.results_folder
except Exception as e:
    results_folder = 'Results'
try:
    Input.executables_folder
except Exception as e:
    executables_folder = 'Executables'
try:
    Input.plotting_folder
except Exception as e:
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

# Check if user specified ftol
try:
    ftol = Input.ftol
except Exception as e:
    ftol = 1e-4

# Run Optimization
try:
    if Input.Optimize:
        # if rel_step_array is specified, perform gradient based optimization
        try:
            Input.rel_step_array
            Input.abs_step_array
            optimize(stel,Input.iota_target,nIterations=Input.nIterations,rel_step_array=Input.rel_step_array,abs_step_array=Input.abs_step_array,grad=True,max_fourier_coefficients=Input.max_fourier_coefficients,ftol=ftol)
        except Exception as e:
            print(e)
            optimize(stel,Input.iota_target,nIterations=Input.nIterations,max_fourier_coefficients=Input.max_fourier_coefficients,ftol=ftol)
except Exception as e:
    print(e)
    Input.Optimize = False

# Check if user specified r_edge
try:
    r_edge = Input.r_edge
except Exception as e:
    r_edge = r_edge

# Do the plotting
try:
    if Input.Plot:
        print('Plotting...')
        # runqsc(stel,name,r_edge,executables_path,plotting_path) # DEPRECATED
        print('  plot()')
        stel.plot(savefig='pyQSC_out.'+name+'.params', show=False)
        print('  B_contour()')
        stel.B_contour(r=r_edge, savefig='pyQSC_out.'+name, ncontours=25, show=False)
        # print('  plot_axis()')
        # stel.plot_axis(savefig='pyQSC_out.'+name+'axis', show=False)
        print('  B_fieldline()')
        stel.B_fieldline(r=r_edge, savefig='pyQSC_out.'+name, show=False)
        print('  plot_boundary()')
        # stel.plot_boundary(r=r_edge, fieldlines=True, savefig='pyQSC_out.'+name+'.boundary', show=False, ntheta=120, nphi=int(120*stel.nfp), ntheta_fourier=30)
        stel.plot_boundary(r=r_edge, fieldlines=False, savefig='pyQSC_out.'+name+'.boundary', show=False, ntheta=120, nphi=int(120*stel.nfp), ntheta_fourier=30)
except Exception as e:
    print(e)
    Input.Plot = False

# Run VMEC
try:
    if Input.VMEC:
        if Input.Optimize:
            input("Copy new optimized configuration into repo")
        print('Outputing to VMEC...')
        stel.to_vmec('Input.'+name,r=r_edge,
                params={"ns_array": [16, 49, 101, 151],
                        "ftol_array": [1e-17,1e-16,1e-15,1e-14],
                        "niter_array": [3000,3000,4000,5000]})
        print('Running VMEC...')
        runVMEC(name,executables_path,plotting_path)
except Exception as e:
    print(e)
    Input.VMEC = False

# Run BOOZ_XFORM
try:
    if Input.BOOZ_XFORM:
        print('Running BOOZ_XFORM...')
        runBOOZXFORM(name)
except Exception as e:
    print(e)
    Input.BOOZ_XFORM = False

# Run NEO
try:
    if Input.NEO:
        print('Running NEO...')
        runNEO(name,executables_path,plotting_path)
except Exception as e:
    print(e)
    Input.NEO = False

# Run SPEC
try:
    if Input.SPEC:
        print('Running SPEC...')
        runSPEC(name,executables_path,plotting_path,stel,r_edge)
except Exception as e:
    print(e)
    Input.SPEC = False

# Run REGCOIL
try:
    if Input.REGCOIL:
        print('Running REGCOIL...')
        runREGCOIL(name,executables_path,plotting_path,coilSeparation = coilSeparation,targetValue = targetValue,nCoilsPerNFP = nCoilsPerNFP)
except Exception as e:
    Input.REGCOIL = False

# Go back to main
os.chdir(main_path)

print('Done')