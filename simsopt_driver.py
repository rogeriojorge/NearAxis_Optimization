from simsopt import LeastSquaresProblem
from simsopt import least_squares_serial_solve
from glob import glob
from os import remove
import numpy as np

def optimize(stel,iota_target=0.41,rel_step_array=[1e-2],abs_step_array=[1e-2],nIterations=20,grad=False):
    ## Print Initial conditions
    print('Before optimization:')
    print('DMerc mean  = ',np.mean(stel.DMerc_times_r2))
    print('DWell mean  = ',np.mean(stel.DWell_times_r2))
    print('Max elongation  = ',stel.max_elongation)
    print('gradgradB inverse length: ', stel.grad_grad_B_inverse_scale_length)
    print('B20 variation =',stel.B20_variation)
    print('Max |X20| =',max(abs(stel.X20)))
    print('Max |X3c1| =',max(abs(stel.X3c1)))
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
    stel.set_fixed('rc(6)', False)
    stel.set_fixed('zs(6)', False)
    stel.set_fixed('rc(7)', False)
    stel.set_fixed('zs(7)', False)
    # stel.set_fixed('zs(8)', False)
    stel.set_fixed('etabar', False)
    stel.set_fixed('B2c', False)
    # stel.set_fixed('p2',False)
    # stel.set_fixed('ds(1)', False)
    # stel.set_fixed('ds(2)', False)
    # stel.set_fixed('ds(3)', False)
    # stel.set_fixed('ds(4)', False)
    # stel.set_fixed('ds(5)', False)
    # stel.set_fixed('ds(6)', False)
    # stel.set_fixed('ds(7)', False)
    # stel.set_fixed('B0(1)', False)
    # stel.set_fixed('alpha0', False)
    # stel.set_fixed('delta', False)

    ## Add integral of J_invariant to optimize for maximum-J
    ## Add NEO to optimize for eps_eff of calculate it analytically
    term = [
            (stel, 'iota', iota_target, 1e6),
            #(stel, 'p2', 0.0, 1e-1),
            (stel, 'max_elongation', 0.0, 3e+0),
            (stel, 'elongation', 0.0, 4e-1),
            (stel, 'B20_anomaly', 0.0, 1e1),
            (stel, 'B20_variation', 0.0, 3e2),
            (stel, 'X20', 0.0, 2e-1),
            (stel, 'X2c', 0.0, 2e-1),
            (stel, 'X2s', 0.0, 2e-1),
            (stel, 'Y20', 0.0, 2e-1),
            (stel, 'Y2c', 0.0, 2e-1),
            (stel, 'Y2s', 0.0, 2e-1),
            (stel, 'Z20', 0.0, 2e-1),
            (stel, 'Z2c', 0.0, 2e-1),
            (stel, 'Z2s', 0.0, 2e-1),
            (stel, 'X3c1', 0.0, 5e-1),
            (stel, 'Y3c1', 0.0, 5e-1),
            (stel, 'Y3s1', 0.0, 5e-1),
            # (stel, 'DMerc_times_r2', 0.5, 1e4),
            # (stel, 'DWell_times_r2', 0.5, 2e3),
            (stel, 'grad_grad_B_inverse_scale_length', 0.0,5e+0),
            (stel.min_R0_penalty, 0.0,1e6)
            # # (stel, 'nlflux_GX', 0.0, 1e1)
            # # (stel, 'd_svals', 0.0, 5e4),
            # # (stel, 'X1c', 0.0, 2e0),
            # # (stel, 'X1s', 0.0, 2e0),
            # # (stel, 'Y1c', 0.0, 2e0),
            # # (stel, 'Y1s', 0.0, 2e0),
            # # # (stel, 'delta', 0.0, 1e4)
            # # (stel, 'curvature', 0.0, 3e2)
            ]

    if grad==False:
        prob = LeastSquaresProblem(term)
        least_squares_serial_solve(prob, max_nfev=nIterations, ftol=1e-10, xtol=1e-10, gtol=1e-10)#, method='lm')
    else:
        for rel_step in rel_step_array:
            for abs_step in abs_step_array:
                if stel.iota == 0: break
                print()
                print(' abs_step =',abs_step)
                print(' rel_step=',rel_step)
                prob = LeastSquaresProblem(term,abs_step=abs_step,rel_step=rel_step,diff_method='centered')

                ## Solve the minimization problem:
                try:
                    stelold = stel
                    probold = prob
                    least_squares_serial_solve(prob, grad=grad, max_nfev=nIterations, ftol=1e-10, xtol=1e-10, gtol=1e-10)#, method='lm')
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
    print('DWell mean  = ',np.mean(stel.DWell_times_r2))
    print('Max elongation  = ',stel.max_elongation)
    print('B20 variation =',stel.B20_variation)
    print('Max |X20| =',max(abs(stel.X20)))
    print('Max |X3c1| =',max(abs(stel.X3c1)))
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
    # print('        B0_vals = [',','.join([str(elem) for elem in stel.B0_vals]),']')
    # print('        d_svals = [',','.join([str(elem) for elem in stel.d_svals]),']')
    # print('        alpha0  =',stel.alpha0)
    # print('        c0      =',stel.c0)
    # print('        delta   =',stel.delta)