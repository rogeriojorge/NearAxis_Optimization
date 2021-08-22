#!/usr/bin/env python3

from stel_repo import get_stel
from util import runqsc, runVMEC, runBOOZXFORM, runNEO, runSPEC, runREGCOIL
from simsopt_driver import optimize
from pathlib import Path
import os

results_folder = 'Results'
executables_folder = 'Executables'
plotting_folder = 'Plotting'

stel, name, r_edge = get_stel(0, nphi=251)

iota_target = 0.41
nIterations = 20
abs_step_array = [1e-2]
rel_step_array = [1e-1]

Path(results_folder+'/'+name).mkdir(parents=True, exist_ok=True)
main_path = str(Path(__file__).parent.resolve())
results_path = str(Path(results_folder+'/'+name).resolve())
executables_path = str(Path(executables_folder).resolve())
plotting_path = str(Path(plotting_folder).resolve())

os.chdir(results_path)

# optimize(stel,iota_target,rel_step_array,abs_step_array,nIterations,grad=False)
# runqsc(stel,name,r_edge,executables_path,plotting_path)
# runVMEC(name,executables_path,plotting_path)
# runBOOZXFORM(name,executables_path,plotting_path)
# runNEO(name,executables_path,plotting_path)
# runSPEC(name,executables_path,plotting_path)
# runREGCOIL(name,executables_path,plotting_path,coilSeparation = 0.2,targetValue = 0.005,nCoilsPerNFP = 6)
# stel.plot(r=0.2,fieldlines=True)

os.chdir(main_path)