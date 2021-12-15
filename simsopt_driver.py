from simsopt import LeastSquaresProblem, least_squares_serial_solve
from glob import glob
from os import remove
import numpy as np
from simsopt.solve.graph_mpi import least_squares_mpi_solve
from simsopt.util.mpi import MpiPartition
from QSCwrapper import QSCWrapper

def optimize(stel,iota_target=0.41,nIterations=20,rel_step_array=[],abs_step_array=[],grad=False,max_fourier_coefficients=5,ftol=1e-3):
    mpi = MpiPartition()
    mpi.write()

    ## Print Initial conditions
    if mpi.proc0_world:
        print('Before optimization:')
        print('Max elongation  = ',stel.max_elongation)
        if stel.order != 'r1':
            print('DMerc mean  = ',np.mean(stel.DMerc_times_r2))
            print('DWell mean  = ',np.mean(stel.DWell_times_r2))
            print('DGeod mean  = ',np.mean(stel.DGeod_times_r2))
            print('B20 variation =',stel.B20_variation)
            print('Max |X20| =',max(abs(stel.X20)))
            print('Max |X3c1| =',max(abs(stel.X3c1)))
            print('gradgradB inverse length: ', stel.grad_grad_B_inverse_scale_length)
            print('B2c     = ',stel.B2c)
            print('B2s     = ',stel.B2s)
            print('p2      = ',stel.p2)
        print('mean gradB inverse length: ', np.mean(stel.inv_L_grad_B))
        print('Max |X1c| =',max(abs(stel.X1c)))
        print('rc      = [',','.join([str(elem) for elem in stel.rc]),']')
        print('zs      = [',','.join([str(elem) for elem in stel.zs]),']')
        print('rs      = [',','.join([str(elem) for elem in stel.rs]),']')
        print('zc      = [',','.join([str(elem) for elem in stel.zc]),']')
        print('etabar  = ',stel.etabar)
        print('nfp     = ',stel.nfp)
        print('nphi    = ',stel.nphi)
        print('iota    = ',stel.iota)

    try:
        stel.omn
    except:
        stel.omn = False
    for n_coeffs in range(1,max_fourier_coefficients+1):
        if mpi.proc0_world:
            print()
            print(' number of Fourier coefficients =',n_coeffs)
        if stel.omn == True:
            if n_coeffs < (len(stel.rc)-1)/2: continue
            if stel.order == 'r1':
                stel = QSCWrapper(sigma0 = stel.sigma0, p_buffer = stel.p_buffer, k_buffer = stel.k_buffer, omn_method = stel.omn_method, rs=stel.rs,zc=stel.zc, rc=stel.rc,zs=stel.zs, nfp=stel.nfp, B0_vals=stel.B0_vals, d_svals=np.append(stel.d_svals,0), nphi=stel.nphi+20, omn=True, delta=stel.delta, d_over_curvature=stel.d_over_curvature, k_second_order_SS=stel.k_second_order_SS)
            else:
                stel = QSCWrapper(sigma0 = stel.sigma0, p_buffer = stel.p_buffer, k_buffer = stel.k_buffer, omn_method = stel.omn_method, rs=stel.rs,zc=stel.zc, rc=stel.rc,zs=stel.zs, nfp=stel.nfp, B0_vals=stel.B0_vals, d_svals=np.append(stel.d_svals,0), nphi=stel.nphi+20, omn=True, delta=stel.delta, B2c_cvals=np.append(stel.B2c_cvals,0), B2s_svals=np.append(stel.B2s_svals,0), p2=stel.p2, order=stel.order, d_over_curvature=stel.d_over_curvature, k_second_order_SS=stel.k_second_order_SS, B2s_cvals=np.append(stel.B2s_cvals,0), B2c_svals=np.append(stel.B2c_svals,0))
            stel.change_qsc_nfourier(2*n_coeffs+1)
        else:
            if n_coeffs < len(stel.rc): continue
            if stel.order =='r1':
                stel = QSCWrapper(sigma0 = stel.sigma0, rc=stel.rc, zs=stel.zs, rs=stel.rs,zc=stel.zc, etabar=stel.etabar, nfp=stel.nfp, nphi=stel.nphi+20, order='r1')
            else:
                stel = QSCWrapper(sigma0 = stel.sigma0, rc=stel.rc, zs=stel.zs, rs=stel.rs,zc=stel.zc, etabar=stel.etabar, nfp=stel.nfp, nphi=stel.nphi+20, B2s=stel.B2s, B2c=stel.B2c, order='r3', p2=stel.p2)
            stel.change_qsc_nfourier(n_coeffs+1)
        try:
            stel.omn
        except:
            stel.omn = False
        stel.min_R0_threshold = 0.4
        stel.fix_all()
        if stel.omn == False:
            if stel.lasym:
                stel.unfix('sigma0')
            for i in range(1,n_coeffs+1):
                stel.unfix('rc('+str(i)+')')
                stel.unfix('zs('+str(i)+')')
                if stel.lasym:
                    stel.unfix('rs('+str(i)+')')
                    stel.unfix('zc('+str(i)+')')
            stel.unfix('etabar')
            if stel.order != 'r1':
                stel.unfix('B2c')
                if stel.lasym:
                    stel.unfix('B2s')
        else:
            stel.unfix('B0(1)')
            if stel.omn_method == 'buffer':
                stel.unfix('delta')
            stel.unfix('zs(2)')
            stel.unfix('rc(2)')
            if stel.order != 'r1':
                stel.unfix('B2cc(0)')
                stel.unfix('B2sc(0)')
            for i in range(1,n_coeffs+1):
                stel.unfix('zs('+str(2*i)+')')
                if i==2: stel.unfix('rc('+str(2*i)+')')
                if stel.d_svals[1] != 0:
                    stel.unfix('ds('+str(i)+')')
                if stel.order != 'r1':
                    stel.unfix('B2ss('+str(i)+')')
                    stel.unfix('B2cc('+str(i)+')')
                    stel.unfix('B2sc('+str(i)+')')
                    stel.unfix('B2cs('+str(i)+')')
            if stel.k_second_order_SS != 0:
                stel.unfix('k_second_order_SS')
            if stel.d_over_curvature != 0:
                stel.unfix('d_over_curvature')
        if stel.order=='r1':
            if stel.omn == False:
                term = [
                        # (stel, 'iota', iota_target, 1e5),
                        # (stel, 'max_elongation', 0.0, 3e+0),
                        (stel.get_elongation, 0.0, 5e-0/stel.nphi),
                        (stel.get_inv_L_grad_B, 0.0, 3e-2),
                        # (stel, 'd_X1c_d_varphi', 0.0, 1e-3),
                        # (stel, 'd_Y1c_d_varphi', 0.0, 1e-3),
                        # (stel, 'd_Y1s_d_varphi', 0.0, 1e-3),
                        # (stel, 'sigma', 0.0, 1e-2),
                        # (stel.min_R0_penalty, 0.0, 1e9),
                ]
            else:
                term = [
                        (stel.get_elongation, 0.0, 5e-1/stel.nphi),
                        # # (stel.get_sigma, 0.0, 1e+1/stel.nphi),
                        # # (stel.get_torsion, 0.0, 1e+1/stel.nphi),
                        # # (stel.get_curvature, 1.0, 5e-1/stel.nphi),
                        (stel.get_d, 0.0, 2e+1/stel.nphi),
                        (stel.get_d_svals, 0.0, 1e2),
                        # # (stel, 'd_X1c_d_varphi', 0.0, 2e-2),
                        # # (stel, 'd_Y1c_d_varphi', 0.0, 2e-2),
                        # # (stel, 'd_Y1s_d_varphi', 0.0, 2-2),
                        (stel.get_min_R0_penalty, 0.0, 3e1),
                        (stel.get_min_Z0_penalty, 0.0, 3e1),
                        # # (stel.get_delta, 0.0, 1e1),
                        (stel.get_B0_well_depth,0.16, 2e2),
                        (stel.get_inv_L_grad_B, 0.0, 5e+0/stel.nphi), # /stel.nphi
                        (stel.get_d_d_d_varphi_at_0,0.0,2e0), # pseudosymmetry?
                        # (stel.get_d_curvature_d_varphi_at_0,0.0,4e-1),
                        # # (stel.get_d_over_curvature,1.0,1e1),
                        # # (stel.get_k_second_order_SS,0.0,5e0),
                        (stel.get_alpha_deviation,0.0,6e+1/stel.nphi)
                ]
        else:
            if stel.omn == False:
                term = [
                        # (stel.get_iota, iota_target, 1e5),
                        # (stel.get_max_elongation, 0.0, 4e-1),
                        (stel.get_elongation, 0.0, 4e-1/stel.nphi),
                        (stel.get_B20_anomaly, 0.0, 5e0/stel.nphi),
                        # (stel.get_B20_variation, 0.0, 4e-1),
                        (stel.get_X20, 0.0, 5e-1/stel.nphi),
                        (stel.get_X2c, 0.0, 5e-1/stel.nphi),
                        (stel.get_X2s, 0.0, 5e-1/stel.nphi),
                        (stel.get_Y20, 0.0, 5e-1/stel.nphi),
                        (stel.get_Y2c, 0.0, 5e-1/stel.nphi),
                        (stel.get_Y2s, 0.0, 5e-1/stel.nphi),
                        # (stel, 'Z20', 0.0, 2e-1),
                        # (stel, 'Z2c', 0.0, 2e-1),
                        # (stel, 'Z2s', 0.0, 2e-1),
                        (stel.get_X3c1, 0.0, 5e-1/stel.nphi),
                        (stel.get_Y3c1, 0.0, 5e-1/stel.nphi),
                        (stel.get_Y3s1, 0.0, 5e-1/stel.nphi),
                        # (stel, 'DMerc_times_r2', 0.3, 3e5),
                        # (stel.get_d2_volume_d_psi2, -50, 1e-1),
                        # (stel, 'DWell_times_r2', 0.1, 1e3),
                        # (stel, 'DGeod_times_r2', 0.1, 1e3),
                        (stel.get_grad_grad_B_inverse_scale_length, 0.0, 5e-2),
                        (stel.get_min_R0_penalty, 0.0, 2e1),
                        (stel.get_min_Z0_penalty, 0.0, 2e1),
                        (stel.get_inv_L_grad_B, 0.0, 5e-1/stel.nphi)
                ]
            else:
                term = [
                        # (stel, 'iota', iota_target, 1e5),
                            # (stel.get_max_elongation, 0.0, 3e-1),
                            # (stel.get_elongation, 0.0, 5e-1/stel.nphi),
                        # (stel.get_B20QI_deviation, 0.0, 3e-2/stel.nphi),
                        # (stel.get_B2cQI_deviation, 0.0, 3e-2/stel.nphi),
                        # (stel.get_B2sQI_deviation, 0.0, 3e-0/stel.nphi),
                        (stel.get_B20QI_deviation_max, 0.0, 1e-0),
                        (stel.get_B2cQI_deviation_max, 0.0, 1e-0),
                        (stel.get_B2sQI_deviation_max, 0.0, 1e-0),
                        # (stel.get_X20, 0.0, 8e-2),
                        # (stel.get_X2c, 0.0, 8e-2),
                        # (stel.get_X2s, 0.0, 8e-2),
                        # (stel.get_Y20, 0.0, 8e-2),
                        # (stel.get_Y2c, 0.0, 8e-2),
                        # (stel.get_Y2s, 0.0, 8e-2),
                        # (stel.get_Z20, 0.0, 8e-2),
                        # (stel.get_Z2c, 0.0, 8e-2),
                        # (stel.get_Z2s, 0.0, 8e-2),
                        # (stel.get_X3c1, 0.0, 8e-2),
                        # (stel.get_X3s1, 0.0, 8e-2),
                        # (stel.get_Y3c1, 0.0, 8e-2),
                        # (stel.get_Y3s1, 0.0, 8e-2),
                        # (stel.get_delta, 0.0, 2e3),
                            # (stel.get_B0_well_depth, 0.15, 1e4),
                        # (stel.get_k_second_order_SS, 0.0, 1e2),
                        # (stel.get_d_svals, 0.0, 1e2),
                        # (stel, 'DMerc_times_r2', 0.3, 3e5),
                        # (stel, 'd2_volume_d_psi2', -50, 1e-1),
                        # (stel, 'DWell_times_r2', 0.1, 1e3),
                        # (stel, 'DGeod_times_r2', 0.1, 1e3),
                            # (stel.get_min_R0_penalty, 0.0, 1e4),
                            # (stel.get_min_Z0_penalty, 0.0, 1e4),
                            # (stel.get_d_curvature_d_varphi_at_0,0.0,1e+0),
                            # (stel.get_inv_L_grad_B, 0.0, 3e+0/stel.nphi),
                            # (stel.get_grad_grad_B_inverse_scale_length_vs_varphi, 0.0, 1e-0/stel.nphi)
                ]

        if grad==False:
            prob = LeastSquaresProblem.from_tuples(term)
            least_squares_serial_solve(prob, max_nfev=nIterations, ftol=ftol)#, ftol=1e-10, xtol=1e-10, gtol=1e-10)#, method='lm')

        else:
            for rel_step in rel_step_array:
                for abs_step in abs_step_array:
                    if stel.iota == 0: break
                    if mpi.proc0_world:
                        print(' abs_step =',abs_step)
                        print(' rel_step=',rel_step)
                    prob = LeastSquaresProblem.from_tuples(term)

                    ## Solve the minimization problem:
                    try:
                        stelold = stel
                        probold = prob
                        # least_squares_serial_solve(prob, grad=grad, max_nfev=nIterations, ftol=1e-10, xtol=1e-10, gtol=1e-10)#, method='lm')
                        least_squares_mpi_solve(prob, mpi, grad=grad, max_nfev=nIterations, ftol=ftol ,abs_step=abs_step,rel_step=rel_step,diff_method='centered')#, ftol=1e-10, xtol=1e-10, gtol=1e-10)#, method='lm')
                    except KeyboardInterrupt:
                        print("Terminated optimization - no change")
                        stel = stelold
                        prob = probold

    if mpi.proc0_world:
        for f in glob("jac_log*.dat"):
            remove(f)
        for f in glob("objective_*.dat"):
            remove(f)
        for f in glob("residuals_*.dat"):
            remove(f)
        for f in glob("simsopt_*.dat"):
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
        print('        rs      = [',','.join([str(elem) for elem in stel.rs]),']')
        print('        zc      = [',','.join([str(elem) for elem in stel.zc]),']')
        print('        sigma0  = ',stel.sigma0)
        if stel.omn:
            if stel.d_svals[-1]==0:
                stel.d_svals = stel.d_svals[0:-1]
            print('        B0_vals = [',','.join([str(elem) for elem in stel.B0_vals]),']')
            print("        omn_method ='"+stel.omn_method+"'")
            print("        k_buffer =",stel.k_buffer)
            print("        p_buffer =",stel.p_buffer)
            print('        k_second_order_SS   =',stel.k_second_order_SS)
            print('        d_over_curvature   =',stel.d_over_curvature)
            print('        d_svals = [',','.join([str(elem) for elem in stel.d_svals]),']')
            print('        delta   =',stel.delta)
            print('        nfp     =',stel.nfp)
            if stel.order == 'r1':
                    print("        stel    =  QSCWrapper(sigma0 = sigma0, omn_method = omn_method, p_buffer = p_buffer, k_buffer=k_buffer, rs=rs,zc=zc, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)")
            else:
                if stel.B2s_svals[-1]==0:
                    stel.B2s_svals = stel.B2s_svals[0:-1]
                print('        B2s_svals = [',','.join([str(elem) for elem in stel.B2s_svals]),']')
                print('        B2c_cvals = [',','.join([str(elem) for elem in stel.B2c_cvals]),']')
                print('        B2s_cvals = [',','.join([str(elem) for elem in stel.B2s_cvals]),']')
                print('        B2c_svals = [',','.join([str(elem) for elem in stel.B2c_svals]),']')
                print('        p2      = ',stel.p2)
                print("        stel    =  QSCWrapper(sigma0 = sigma0, omn_method = omn_method, p_buffer = p_buffer, k_buffer=k_buffer, rs=rs,zc=zc, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, B2c_cvals=B2c_cvals, B2s_svals=B2s_svals, p2=p2, order='r3', k_second_order_SS=k_second_order_SS, d_over_curvature=d_over_curvature, B2s_cvals=B2s_cvals, B2c_svals=B2c_svals)")
        else:
            print('        etabar = ',stel.etabar)
            print('        nfp    = ',stel.nfp)
            if stel.order == 'r1':
                print("        stel   = QSCWrapper(sigma0 = sigma0, rc=rc, zs=zs, rs=rs,zc=zc, etabar=etabar, nfp=nfp, nphi=nphi, order='r1')")
            else:
                print('        B2c    = ',stel.B2c)
                print('        B2s    = ',stel.B2s)
                print('        p2     = ',stel.p2)
                print("        stel   = QSCWrapper(sigma0 = sigma0, rc=rc, zs=zs, rs=rs,zc=zc, etabar=etabar, nfp=nfp, nphi=nphi, B2s=B2s, B2c=B2c, order='r3', p2=p2)")
        print('        iota    = ',stel.iota)
        if stel.order != 'r1':
            print('        # DMerc mean  =',np.mean(stel.DMerc_times_r2))
            print('        # DWell mean  =',np.mean(stel.DWell_times_r2))
            print('        # DGeod mean  =',np.mean(stel.DGeod_times_r2))
            print('        # B20 mean =',np.mean(stel.B20))
            if stel.omn:
                print('        # B20QI_deviation_max =',stel.B20QI_deviation_max)
                print('        # B2cQI_deviation_max =',stel.B2cQI_deviation_max)
                print('        # B2sQI_deviation_max =',stel.B2sQI_deviation_max)
            else:
                print('        # B20 variation =',stel.B20_variation)
            print('        # Max |X20| =',max(abs(stel.X20)))
            if stel.order == 'r3':
                print('        # Max |X3c1| =',max(abs(stel.X3c1)))
            print('        # gradgradB inverse length:', stel.grad_grad_B_inverse_scale_length)
            print('        # d2_volume_d_psi2 mean =',np.mean(stel.d2_volume_d_psi2))
        if stel.omn:
            print("        # max curvature'(0):", stel.d_curvature_d_varphi_at_0)
            print("        # max d'(0):", stel.d_d_d_varphi_at_0)
        print('        # max gradB inverse length:', np.max(stel.inv_L_grad_B))
        print('        # Max elongation =',stel.max_elongation)
        print('        # objective function: ', prob.objective())
