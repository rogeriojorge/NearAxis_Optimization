from simsopt import LeastSquaresProblem
from simsopt import least_squares_serial_solve
from glob import glob
from os import remove
import numpy as np

def optimize(stel,iota_target=0.41,rel_step_array=[1e-2],abs_step_array=[1e-2],nIterations=20,grad=False):
    ## Print Initial conditions
    print('Before optimization:')
    print('DMerc mean  = ',np.mean(stel.DMerc_times_r2))
    print('Max elongation  = ',stel.max_elongation)
    print('gradgradB inverse length: ', stel.grad_grad_B_inverse_scale_length)
    print('B20 variation =',stel.B20_variation)
    print('Max X20 =',max(stel.X20))
    print('rc      = [',','.join([str(elem) for elem in stel.rc]),']')
    print('zs      = [',','.join([str(elem) for elem in stel.zs]),']')
    print('etabar  = ',stel.etabar)
    print('nfp     = ',stel.nfp)
    print('B2c     = ',stel.B2c)
    print('p2      = ',stel.p2)
    print('nphi    = ',stel.nphi)
    print('iota    = ',stel.iota)

    stel.all_fixed()
    stel.set_fixed('rc(1)', False) 
    stel.set_fixed('zs(1)', False)
    stel.set_fixed('rc(2)', False) 
    stel.set_fixed('zs(2)', False)
    stel.set_fixed('rc(3)', False) 
    stel.set_fixed('zs(3)', False)
    stel.set_fixed('rc(4)', False)
    stel.set_fixed('zs(4)', False)
    stel.set_fixed('rc(5)', False)
    stel.set_fixed('zs(5)', False)
    # stel.set_fixed('rc(6)', False)
    # stel.set_fixed('zs(6)', False)
    stel.set_fixed('etabar', False)
    stel.set_fixed('B2c', False)
    # stel.set_fixed('p2',False)

    term = [
            (stel, 'iota', iota_target, 1e5),
            #(stel, 'p2', 0.0, 1e-1),
            (stel, 'max_elongation', 0.0, 1e+1),
            (stel, 'elongation', 0.0, 1e0),
            (stel, 'B20_anomaly', 0.0, 2e2),
            (stel, 'B20_variation', 0.0, 1e3),
            (stel, 'X20', 0.0, 5e-1),
            (stel, 'X2c', 0.0, 5e-1),
            (stel, 'X2s', 0.0, 5e-1),
            (stel, 'Y20', 0.0, 5e-1),
            (stel, 'Y2c', 0.0, 5e-1),
            (stel, 'Y2s', 0.0, 5e-1),
            (stel, 'Z20', 0.0, 5e-1),
            (stel, 'Z2c', 0.0, 5e-1),
            (stel, 'Z2s', 0.0, 5e-1),
            # (stel, 'DMerc_times_r2', 0.01,1e6),
            (stel, 'grad_grad_B_inverse_scale_length', 0.0,1e-3)
            ]

    for rel_step in rel_step_array:
        for abs_step in abs_step_array:
            print()
            print(' abs_step =',abs_step)
            print(' rel_step=',rel_step)
            prob = LeastSquaresProblem(term,abs_step=abs_step,rel_step=rel_step,diff_method='centered')

            ## Solve the minimization problem:
            try:
                stelold = stel
                probold = prob
                least_squares_serial_solve(prob, grad=grad, max_nfev=nIterations)#, method='lm')
            except KeyboardInterrupt:
                print("Terminated optimization - no change")
                stel = stelold
                prob = probold
            for f in glob("residuals_2021*.dat"):
                remove(f)
            for f in glob("simsopt_2021*.dat"):
                remove(f)

    ## Print final conditions
    print('After optimization:')
    print('DMerc mean  = ',np.mean(stel.DMerc_times_r2))
    print('Max elongation  = ',stel.max_elongation)
    print('B20 variation =',stel.B20_variation)
    print('Max X20 =',max(stel.X20))
    print('gradgradB inverse length: ', stel.grad_grad_B_inverse_scale_length)
    print('objective function: ', prob.objective())
    nN=stel.iota-stel.iotaN
    if nN==0:
        print('Quasi-axisymmetric solution')
    else:
        print('Quasi-helically symmetric solution with N =',nN)
    print('        rc     = [',','.join([str(elem) for elem in stel.rc]),']')
    print('        zs     = [',','.join([str(elem) for elem in stel.zs]),']')
    print('        etabar = ',stel.etabar)
    print('        nfp    = ',stel.nfp)
    print('        B2c    = ',stel.B2c)
    print('        p2     = ',stel.p2)
    print('        iota   = ',stel.iota)