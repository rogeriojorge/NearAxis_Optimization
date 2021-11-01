from simsopt import LeastSquaresProblem
from simsopt import least_squares_serial_solve
from simsopt import make_optimizable
from qsc import Qsc
from glob import glob
from os import remove
import numpy as np
from simsopt.solve.mpi import least_squares_mpi_solve
from simsopt.util.mpi import MpiPartition

def optimize(stel,iota_target=0.41,nIterations=20,rel_step_array=[],abs_step_array=[],grad=False,max_fourier_coefficients=5):
    mpi = MpiPartition()
    mpi.write()

    ## Print Initial conditions
    if mpi.proc0_world:
        print('Before optimization:')
        print('Max elongation  = ',stel.max_elongation)
        if stel.order == 'r3':
            print('DMerc mean  = ',np.mean(stel.DMerc_times_r2))
            print('DWell mean  = ',np.mean(stel.DWell_times_r2))
            print('DGeod mean  = ',np.mean(stel.DGeod_times_r2))
            print('B20 variation =',stel.B20_variation)
            print('Max |X20| =',max(abs(stel.X20)))
            print('Max |X3c1| =',max(abs(stel.X3c1)))
            print('gradgradB inverse length: ', stel.grad_grad_B_inverse_scale_length)
            print('B2c     = ',stel.B2c)
            print('p2      = ',stel.p2)
        print('mean gradB inverse length: ', np.mean(stel.inv_L_grad_B))
        print('Max |X1c| =',max(abs(stel.X1c)))
        print('rc      = [',','.join([str(elem) for elem in stel.rc]),']')
        print('zs      = [',','.join([str(elem) for elem in stel.zs]),']')
        print('etabar  = ',stel.etabar)
        print('nfp     = ',stel.nfp)
        print('nphi    = ',stel.nphi)
        print('iota    = ',stel.iota)

    for n_coeffs in range(1,max_fourier_coefficients+1):
        if mpi.proc0_world:
            print()
            print(' number of Fourier coefficients =',n_coeffs)
        if stel.omn == True:
            if n_coeffs < len(stel.d_svals)-1: continue
            if stel.k_second_order_SS == 0:
                if stel.order == 'r1':
                    stel = Qsc(rc=stel.rc,zs=stel.zs, nfp=stel.nfp, B0_vals=stel.B0_vals, d_svals=np.append(stel.d_svals,0), nphi=stel.nphi+20, omn=True, delta=stel.delta)
                else:
                    stel = Qsc(rc=stel.rc,zs=stel.zs, nfp=stel.nfp, B0_vals=stel.B0_vals, d_svals=np.append(stel.d_svals,0), nphi=stel.nphi+20, omn=True, delta=stel.delta, B2s_svals=np.append(stel.B2s_svals,0), p2=stel.p2, order='r2')
            else:
                if stel.order == 'r1':
                    stel = Qsc(rc=stel.rc,zs=stel.zs, nfp=stel.nfp, B0_vals=stel.B0_vals, d_svals=np.append(stel.d_svals,0), nphi=stel.nphi+20, omn=True, delta=stel.delta, k_second_order_SS=stel.k_second_order_SS)
                else:
                    stel = Qsc(rc=stel.rc,zs=stel.zs, nfp=stel.nfp, B0_vals=stel.B0_vals, d_svals=np.append(stel.d_svals,0), nphi=stel.nphi+20, omn=True, delta=stel.delta, B2s_svals=np.append(stel.B2s_svals,0), p2=stel.p2, order='r2', k_second_order_SS=stel.k_second_order_SS)
            stel.change_nfourier(2*n_coeffs+1)
        else:
            if n_coeffs < len(stel.rc): continue
            if stel.order =='r1':
                stel = Qsc(rc=stel.rc, zs=stel.zs, etabar=stel.etabar, nfp=stel.nfp, nphi=stel.nphi+20, order='r1')
            else:
                stel = Qsc(rc=stel.rc, zs=stel.zs, etabar=stel.etabar, nfp=stel.nfp, nphi=stel.nphi+20, B2c=stel.B2c, order='r3', p2=stel.p2)
            stel.change_nfourier(n_coeffs+1)
        try:
            stel.omn
        except:
            stel.omn = False
        stel.min_R0_threshold = 0.6
        stel=make_optimizable(stel)
        stel.all_fixed()
        if stel.omn == False:
            for i in range(1,n_coeffs+1):
                stel.set_fixed('rc('+str(i)+')', False)
                stel.set_fixed('zs('+str(i)+')', False)
            stel.set_fixed('etabar', False)
            if stel.order != 'r1':
                stel.set_fixed('B2c',False)
        else:
            for i in range(1,n_coeffs+1):
                stel.set_fixed('zs('+str(2*i)+')', False)
                stel.set_fixed('B0(1)', False)
                stel.set_fixed('delta', False)
                if stel.k_second_order_SS != 0:
                    stel.set_fixed('k_second_order_SS', False)
                else:
                    stel.set_fixed('ds('+str(i)+')', False)
                if stel.order != 'r1':
                    stel.set_fixed('B2ss('+str(i)+')', False)
        if stel.order=='r1':
            if stel.omn == False:
                term = [
                        # (stel, 'iota', iota_target, 1e5),
                        # (stel, 'max_elongation', 0.0, 3e+0),
                        # (stel, 'elongation', 0.0, 4e-1),
                        (stel, 'inv_L_grad_B', 0.0, 1e-1),
                        # (stel, 'd_X1c_d_varphi', 0.0, 1e-3),
                        # (stel, 'd_Y1c_d_varphi', 0.0, 1e-3),
                        # (stel, 'd_Y1s_d_varphi', 0.0, 1e-3),
                        # (stel, 'sigma', 0.0, 1e-2),
                        (stel.min_R0_penalty, 0.0, 1e9),
                ]
            else:
                term = [
                        # (stel, 'iota', -6.95, 1e4),
                        (stel, 'max_elongation', 0.0, 5e-1),
                        (stel, 'elongation', 0.0, 5e-2),
                        # (stel, 'sigma', 0.0, 1e-1),
                        # (stel, 'torsion', 0.0, 3e-2),
                        # (stel, 'curvature', 1/stel.rc[0], 1e-2),
                        # (stel, 'd', 0.0, 1e-1),
                        (stel, 'd_svals', 0.0, 1e2),
                        # (stel, 'd_X1c_d_varphi', 0.0, 2e-2),
                        # (stel, 'd_Y1c_d_varphi', 0.0, 2e-2),
                        # (stel, 'd_Y1s_d_varphi', 0.0, 2e-2),
                        # (stel.min_R0_penalty, 0.0, 1e9),
                        (stel, 'delta', 0.0, 7e1),
                        (stel, 'B0_well_depth', 0.15, 2e3),
                        (stel, 'inv_L_grad_B', 0.0, 4e-2),
                ]
        else:
            if stel.omn == False:
                term = [
                        (stel, 'iota', iota_target, 1e5),
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
                        # (stel, 'DMerc_times_r2', 0.3, 3e5),
                        (stel, 'd2_volume_d_psi2', -50, 1e-1),
                        # (stel, 'DWell_times_r2', 0.1, 1e3),
                        # (stel, 'DGeod_times_r2', 0.1, 1e3),
                        (stel, 'grad_grad_B_inverse_scale_length', 0.0, 1e+1),
                        (stel.min_R0_penalty, 0.0, 1e9),
                        (stel, 'inv_L_grad_B', 0.0, 1e0)
                ]
            else:
                term = [
                        # (stel, 'iota', iota_target, 1e5),
                        (stel, 'max_elongation', 0.0, 3e+0),
                        (stel, 'elongation', 0.0, 4e-1),
                        (stel, 'B20', 0.0, 1e4),
                        (stel, 'B20_variation', 0.0, 1e3),
                        (stel, 'X20', 0.0, 1e-0),
                        (stel, 'X2c', 0.0, 1e-0),
                        (stel, 'X2s', 0.0, 1e-0),
                        (stel, 'Y20', 0.0, 1e-0),
                        (stel, 'Y2c', 0.0, 1e-0),
                        (stel, 'Y2s', 0.0, 1e-0),
                        (stel, 'Z20', 0.0, 1e-0),
                        (stel, 'Z2c', 0.0, 1e-0),
                        (stel, 'Z2s', 0.0, 1e-0),
                        (stel, 'delta', 0.0, 5e3),
                        (stel, 'B0_well_depth', 0.15, 2e3),
                        # (stel, 'DMerc_times_r2', 0.3, 3e5),
                        # (stel, 'd2_volume_d_psi2', -50, 1e-1),
                        # (stel, 'DWell_times_r2', 0.1, 1e3),
                        # (stel, 'DGeod_times_r2', 0.1, 1e3),
                        # (stel, 'grad_grad_B_inverse_scale_length', 0.0, 1e+1),
                        # (stel.min_R0_penalty, 0.0, 1e9),
                        # (stel, 'inv_L_grad_B', 0.0, 1e0)
                ]

        if grad==False:
            prob = LeastSquaresProblem(term)
            least_squares_serial_solve(prob, max_nfev=nIterations)#, ftol=1e-10, xtol=1e-10, gtol=1e-10)#, method='lm')

        else:
            for rel_step in rel_step_array:
                for abs_step in abs_step_array:
                    if stel.iota == 0: break
                    if mpi.proc0_world:
                        print(' abs_step =',abs_step)
                        print(' rel_step=',rel_step)
                    prob = LeastSquaresProblem(term,abs_step=abs_step,rel_step=rel_step,diff_method='centered')

                    ## Solve the minimization problem:
                    try:
                        stelold = stel
                        probold = prob
                        # least_squares_serial_solve(prob, grad=grad, max_nfev=nIterations, ftol=1e-10, xtol=1e-10, gtol=1e-10)#, method='lm')
                        least_squares_mpi_solve(prob, mpi, grad=grad, max_nfev=nIterations)#, ftol=1e-10, xtol=1e-10, gtol=1e-10)#, method='lm')
                    except KeyboardInterrupt:
                        print("Terminated optimization - no change")
                        stel = stelold
                        prob = probold

    if mpi.proc0_world:
        for f in glob("residuals_2021*.dat"):
            remove(f)
        for f in glob("simsopt_2021*.dat"):
            remove(f)

    ## Print final conditions
    if mpi.proc0_world:
        print('After optimization:')
        nN=stel.iota-stel.iotaN
        if stel.omn:
            print('Quasi-isodynamic solution')
        else:
            if nN==0:
                print('Quasi-axisymmetric solution')
            else:
                print('Quasi-helically symmetric solution with N =',nN)
        print('        rc      = [',','.join([str(elem) for elem in stel.rc]),']')
        print('        zs      = [',','.join([str(elem) for elem in stel.zs]),']')
        if stel.omn:
            if stel.d_svals[-1]==0:
                stel.d_svals = stel.d_svals[0:-1]
            print('        B0_vals = [',','.join([str(elem) for elem in stel.B0_vals]),']')
            print('        d_svals = [',','.join([str(elem) for elem in stel.d_svals]),']')
            if stel.k_second_order_SS != 0:
                print('        k_second_order_SS   =',stel.k_second_order_SS)
            print('        delta   =',stel.delta)
            print('        nfp     =',stel.nfp)
            if stel.order == 'r1':
                if stel.k_second_order_SS == 0:
                    print("        stel    =  make_optimizable(Qsc(rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta))")
                else:
                    print("        stel    =  make_optimizable(Qsc(rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, k_second_order_SS=k_second_order_SS))")
            else:
                if stel.B2s_svals[-1]==0:
                    stel.B2s_svals = stel.B2s_svals[0:-1]
                print('        B2s_svals = [',','.join([str(elem) for elem in stel.B2s_svals]),']')
                print('        p2      = ',stel.p2)
                if stel.k_second_order_SS == 0:
                    print("        stel    =  make_optimizable(Qsc(rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, B2s_svals=B2s_svals, p2=p2, order='r2'))")
                else:
                    print("        stel    =  make_optimizable(Qsc(rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, B2s_svals=B2s_svals, p2=p2, order='r2', k_second_order_SS=k_second_order_SS))")
        else:
            print('        etabar = ',stel.etabar)
            print('        nfp    = ',stel.nfp)
            if stel.order == 'r1':
                print("        stel   = make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, order='r1'))")
            else:
                print('        B2c    = ',stel.B2c)
                print('        p2     = ',stel.p2)
                print("        stel   = make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))")
        print('        iota    = ',stel.iota)
        if stel.order != 'r1':
            print('        # DMerc mean  =',np.mean(stel.DMerc_times_r2))
            print('        # DWell mean  =',np.mean(stel.DWell_times_r2))
            print('        # DGeod mean  =',np.mean(stel.DGeod_times_r2))
            print('        # B20 variation =',stel.B20_variation)
            print('        # B20 mean =',np.mean(stel.B20))
            print('        # Max |X20| =',max(abs(stel.X20)))
            if stel.order == 'r3':
                print('        # Max |X3c1| =',max(abs(stel.X3c1)))
            print('        # gradgradB inverse length:', stel.grad_grad_B_inverse_scale_length)
            print('        # d2_volume_d_psi2 mean =,',np.mean(stel.d2_volume_d_psi2))
        print('        # mean gradB inverse length:', np.mean(stel.inv_L_grad_B))
        print('        # Max elongation =',stel.max_elongation)
        print('        # objective function: ', prob.objective())
