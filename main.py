#!/usr/bin/env python3

from stel_repo import get_stel
from util import runqsc, runVMEC, runBOOZXFORM, runNEO, runSPEC, runREGCOIL, runSTAGE2, runVMECfree, runBEAMS3D, runVMECrescale
from simsopt_driver import optimize
from pathlib import Path
import os, sys
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
    # print(e)
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
        # print('  B_contour()')
        # stel.B_contour(r=r_edge, savefig='pyQSC_out.'+name, ncontours=25, show=False)
        # print('  B_fieldline()')
        # stel.B_fieldline(r=r_edge, savefig='pyQSC_out.'+name, show=False)
        # print('  plot_boundary()')
        # stel.plot_boundary(r=r_edge, fieldlines=True, savefig='pyQSC_out.'+name+'.boundary', show=False, ntheta=120, nphi=int(120*stel.nfp), ntheta_fourier=30)
        # stel.plot_boundary(r=r_edge, fieldlines=False, savefig='pyQSC_out.'+name+'.boundary', show=False, ntheta=120, nphi=int(120*stel.nfp), ntheta_fourier=30)
        # print('  plot_axis()')
        # stel.plot_axis(frenet_factor=0.15, savefig='pyQSC_out_axis.'+name, show=False)
        # sys.path.insert(1, plotting_path)
        # import Simple3Dplot
        # Simple3Dplot.main(name, stel, r_edge, show=False)
except Exception as e:
    # print(e)
    Input.Plot = False

# Run VMEC
try:
    if Input.VMEC:
        if Input.Optimize:
            input("Copy new optimized configuration into repo")
        print('Outputing to VMEC...')
        stel.to_vmec('input.'+name,r=r_edge,
                params={"ns_array": [16, 49, 101, 151],
                        "ftol_array": [1e-11,1e-12,1e-13,1e-14],
                        "niter_array": [1000, 2000, 2000, 3000]})
        print('Running VMEC...')
        runVMEC(name,stel,executables_path,plotting_path)
except Exception as e:
    # print(e)
    Input.VMEC = False

# Run BOOZ_XFORM
try:
    if Input.BOOZ_XFORM:
        print('Running BOOZ_XFORM...')
        runBOOZXFORM(name)
        stel_from_boozxform = stel.from_boozxform('boozmn_'+name+'.nc', max_s_for_fit = 0.4, N_phi = stel.nphi,
                                max_n_to_plot = 2, show=False, vmec_file='wout_'+name+'.nc', input_stel=stel, nNormal=stel.iotaN-stel.iota)
except Exception as e:
    # print(e)
    Input.BOOZ_XFORM = False

# Run NEO
try:
    if Input.NEO:
        print('Running NEO...')
        runNEO(name,executables_path,plotting_path)
except Exception as e:
    # print(e)
    Input.NEO = False

# Run SPEC
try:
    if Input.SPEC:
        print('Running SPEC...')
        runSPEC(name,executables_path,plotting_path,stel,r_edge)
except Exception as e:
    # print(e)
    Input.SPEC = False

# Run REGCOIL
# Check if user specified coil parameters
try:
    if Input.REGCOIL:
        try:
            coilSeparation = Input.coilSeparation
        except Exception as e:
            coilSeparation = 0.1
        try:
            targetValue = Input.targetValue
        except Exception as e:
            targetValue = 0.08
        try:
            nCoilsPerNFP = Input.nCoilsPerNFP
        except Exception as e:
            nCoilsPerNFP = 6
        print('Running REGCOIL...')
        runREGCOIL(name,stel,r_edge,executables_path,plotting_path,coilSeparation = coilSeparation,targetValue = targetValue,nCoilsPerNFP = nCoilsPerNFP)
except Exception as e:
    # print(e)
    Input.REGCOIL = False

# Run STAGE2
try:
    if Input.STAGE2:
        print('Running STAGE2...')
        runSTAGE2(name, plotting_path, stel, r_edge, run_get_coils=True)
except Exception as e:
    # print(e)
    Input.STAGE2 = False

# Run VMEC free boundary
try:
    if Input.VMECfree:
        print('Running VMEC free boundary...')
        runVMECfree(name, stel, executables_path, plotting_path)
except Exception as e:
    # print(e)
    Input.VMECfree = False

# Run BOOZ_XFORM
try:
    if Input.BOOZ_XFORM_free:
        print('Running BOOZ_XFORM free boundary...')
        runBOOZXFORM(name+"_free")
        stel_from_boozxform = stel.from_boozxform('boozmn_'+name+'_free.nc', max_s_for_fit = 0.4, N_phi = stel.nphi,
                        max_n_to_plot = 2, show=False, vmec_file='wout_'+name+'_free.nc', input_stel=stel, nNormal=stel.helicity)
except Exception as e:
    # print(e)
    Input.BOOZ_XFORM_free = False

# Run NEO
try:
    if Input.NEO_free:
        print('Running NEO free boundary...')
        runNEO(name+"_free",executables_path,plotting_path)
except Exception as e:
    # print(e)
    Input.NEO_free = False

# Check if user specified VMEC rescaling parameters
try:
    B_scale = Input.B_scale
    R_scale = Input.R_scale
except Exception as e:
    B_scale = 6
    R_scale = 15
# Run VMEC rescaled
try:
    if Input.VMECrescale:
        print('Running VMEC rescaled...')
        runVMECrescale(name, stel, executables_path, plotting_path, B_scale=B_scale, R_scale=R_scale)
except Exception as e:
    # print(e)
    Input.VMECrescale = False

# Check if user specified BEAMS3D parameters
try:
    runBEAMS = Input.runBEAMS
    nparticles = Input.nparticles
    s0 = Input.s0
    T_END_IN = Input.T_END_IN
except Exception as e:
    runBEAMS = True
    nparticles = 100000
    s0 = 1e-07
    T_END_IN = 1e-3
# Run BEAMS3D
try:
    if Input.BEAMS3D:
        print('Running BEAMS3D...')
        runBEAMS3D(name,executables_path,plotting_path,runBEAMS,
            nparticles=nparticles,s0=s0,T_END_IN=T_END_IN)
except Exception as e:
    # print(e)
    Input.BEAMS3D = False

# Go back to main
os.chdir(main_path)

print('Done')