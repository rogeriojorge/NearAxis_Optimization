from numpy import savetxt
from subprocess import run
import sys
import shutil
import matplotlib.pyplot as plt

## Function to replace text in files
def replace(file_path, pattern, subst):
    #Create temp file
    from tempfile import mkstemp
    from os import fdopen, remove
    from shutil import copymode, move
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    #Copy the file permissions from the old file to the new file
    copymode(file_path, abs_path)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)
    return 0

## Output and run quasisymmetry Fortran code
# Create input file
def output2qsc(stel,name,rr,executables_path):
    vmecExample='"'+executables_path+'/quasisymmetry.VMECexample"'
    filename = 'quasisymmetry_in.'+name
    with open(filename, "w") as txt_file:
        txt_file.write('&quasisymmetry\n')
        txt_file.write(' general_option="single"\n')
        txt_file.write(' finite_r_option="nonlinear"\n')
        if stel.order=='r1':
            txt_file.write(' order_r_option="r1"\n')
        else:
            txt_file.write(' order_r_option="r3_flux_constraint"\n')
        txt_file.write(' N_phi=451\n')
        txt_file.write(' nfp = '+str(stel.nfp)+'\n')
        txt_file.write(' eta_bar = '+str(stel.etabar)+'\n')
        txt_file.write(' R0c = ')
        savetxt(txt_file, stel.rc, newline=',')
        txt_file.write('\n')
        txt_file.write(' Z0s = ')
        savetxt(txt_file, stel.zs, newline=',')
        txt_file.write('\n')
        txt_file.write(' B2c = '+str(stel.B2c)+'\n')
        txt_file.write(' B2s = '+str(stel.B2s)+'\n')
        txt_file.write(' sigma_initial = '+str(stel.sigma0)+'\n')
        txt_file.write(' I2_over_B0 = '+str(stel.I2)+'\n')
        txt_file.write(' p2 = '+str(stel.p2)+'\n')
        txt_file.write(' r = '+str(rr)+'\n')
        txt_file.write(' vmec_template_filename = '+vmecExample+'\n')
        txt_file.write(' r_singularity_high_order = T'+'\n')
        txt_file.write(' r_singularity_Newton_iterations = 20'+'\n')
        txt_file.write(' r_singularity_line_search = 10'+'\n')
        txt_file.write('/\n')
# Run
def runqsc(stel,name,rr,executables_path,plotting_path):
    print("Output to quasisymmetry code")
    output2qsc(stel,name,rr,executables_path)
    print("Run quasisymmetry code")
    run([executables_path+"/./quasisymmetry", 'quasisymmetry_in.'+name])
    replace("input."+name,'NTOR = 0101','NTOR = 0014')
    print("Plot quasisymmetry result")
    sys.path.insert(1, plotting_path)
    import quasisymmetryPlotSingle
    quasisymmetryPlotSingle.main("quasisymmetry_out."+name+".nc")

# Run VMEC code
def runVMEC(name,executables_path,plotting_path):
    print("Run VMEC")
    # import vmec
    # import numpy as np
    # from mpi4py import MPI
    # ictrl = np.zeros(5, dtype=np.int32)
    # verbose = True
    # reset_file = ''

    # # Flags used by runvmec():
    # restart_flag = 1
    # readin_flag = 2
    # timestep_flag = 4
    # output_flag = 8
    # cleanup_flag = 16
    # reset_jacdt_flag = 32

    # fcomm = MPI.COMM_WORLD.py2f()
    # ictrl[:] = 0
    # ictrl[0] = restart_flag + readin_flag + timestep_flag + output_flag + cleanup_flag
    # print("Calling runvmec. ictrl={} comm={}".format(ictrl, fcomm))
    # vmec.runvmec(ictrl, "input."+name, verbose, fcomm, reset_file)

    # np.testing.assert_equal(ictrl[1], 0)

    bashCommand = executables_path+"/./xvmec2000 input."+name
    # run(bashCommand.split())
    print("Plot VMEC result")
    sys.path.insert(1, plotting_path)
    import vmecPlot2
    vmecPlot2.main("wout_"+name+".nc")

# Run booz_xform
def runBOOZXFORM(name):
    import booz_xform as bx
    print("Run BOOZ_XFORM")
    b1 = bx.Booz_xform()
    b1.read_wout("wout_"+name+".nc")
    b1.compute_surfs = [0,2,5,8,16,25,32,40,50,60,70,80,90,100,110,125,135,145]
    b1.mboz = 150
    b1.nboz = 60
    b1.run()
    b1.write_boozmn("boozmn_"+name+".nc")
    # b1.read_boozmn("boozmn_"+name+".nc")
    print("Plot BOOZ_XFORM")
    fig = plt.figure(); bx.surfplot(b1, js=5,  fill=False, ncontours=35)
    plt.savefig("Boozxform_surfplot_1_"+name+'.pdf', bbox_inches = 'tight', pad_inches = 0); plt.close()
    fig = plt.figure(); bx.surfplot(b1, js=11, fill=False, ncontours=35)
    plt.savefig("Boozxform_surfplot_2_"+name+'.pdf', bbox_inches = 'tight', pad_inches = 0); plt.close()
    fig = plt.figure(); bx.surfplot(b1, js=16, fill=False, ncontours=35)
    plt.savefig("Boozxform_surfplot_3_"+name+'.pdf', bbox_inches = 'tight', pad_inches = 0); plt.close()
    if name[0:2] == 'QH':
        helical_detail = True
    else:
        helical_detail = False
    fig = plt.figure(); bx.symplot(b1, helical_detail = helical_detail, sqrts=True)
    plt.savefig("Boozxform_symplot_"+name+'.pdf', bbox_inches = 'tight', pad_inches = 0); plt.close()
    fig = plt.figure(); bx.modeplot(b1, sqrts=True)
    plt.savefig("Boozxform_modeplot_"+name+'.pdf', bbox_inches = 'tight', pad_inches = 0); plt.close()
    fig = bx.wireplot(b1, orig = False)
    fig.write_image("Boozxform_wireplot_"+name+'.pdf')

# Run NEO
def runNEO(name,executables_path,plotting_path):
    print("Run NEO")
    shutil.copy(executables_path+'/neo_in.example', 'neo_in.'+name)
    replace("neo_in."+name,'example',name)
    bashCommand = executables_path+"/./xneo "+name
    run(bashCommand.split())
    print("Plot NEOresult")
    sys.path.insert(1, plotting_path)
    import neoPlot
    neoPlot.main("neo_out."+name,name)

## Run SPEC
def output2spec(qvfilename,qscfile,executables_path,nmodes,stel,r_edge):
    stel.to_vmec('Input.'+qvfilename,r=r_edge,
            params={"ns_array": [16, 49, 101, 151],
                    "ftol_array": [1e-17,1e-16,1e-15,1e-14],
                    "niter_array": [3000,3000,4000,5000]})
    RBC=stel.RBC
    ZBS=stel.ZBS
    rc=stel.rc
    zs=stel.zs
    nfp=stel.nfp
    with open("input."+qvfilename, 'r') as read_obj:
        for line in read_obj:
            if "PHIEDGE" in line:
                PHIEDGE=line.split()[2]
    rbc=[0 for i in range(len(RBC))]
    zbs=[0 for i in range(len(RBC))]
    nmodesTot=len(RBC[1])
    for count in range(len(RBC)):
        if count==0:
            val = next((index for index,value in enumerate(RBC[0]) if value != 0), None)
            rbc[count]=RBC[count][val:val+2*nmodes]
            val = next((index for index,value in enumerate(ZBS[0]) if value != 0), None)
            zbs[count]=ZBS[count][val-1:val-1+2*nmodes]
        else:
            rbc[count]=RBC[count][int((nmodesTot-1)/2-nmodes):int((nmodesTot-1)/2+nmodes)]
            zbs[count]=ZBS[count][int((nmodesTot-1)/2-nmodes):int((nmodesTot-1)/2+nmodes)]
    text="! Magnetic axis shape:\n"
    text=text+"rac(0:"+str(len(rc)-1)+")="+", ".join([str(elem) for elem in rc])+"\n\n"
    text=text+"zas(0:"+str(len(zs)-1)+")="+", ".join([str(-elem) for elem in zs])+"\n"
    shutil.copyfile(executables_path+"/input.SPECexample.sp", qvfilename+".sp")
    replace(qvfilename+".sp","! Magnetic axis shape:",text)
    text="!----- Boundary Parameters -----\n"
    for countn in range(len(rbc)-1):
        if countn==0:
            for countm in range(nmodes):
                text=text+"rbc("+str(countm)+","+str(countn)+")= "+str(rbc[countn][countm])+", zbs("+str(countm)+","+str(countn)+")= "+str(zbs[countn][countm])+",\n"
        elif countn==len(rbc)-2:
            for countm in range(2*nmodes):
                if countm==2*nmodes-1:
                    text=text+"rbc("+str(countm-nmodes)+","+str(countn)+")= "+str(rbc[countn][countm])+", zbs("+str(countm-nmodes)+","+str(countn)+")= "+str(zbs[countn][countm])+"\n"
                else:
                    text=text+"rbc("+str(countm-nmodes)+","+str(countn)+")= "+str(rbc[countn][countm])+", zbs("+str(countm-nmodes)+","+str(countn)+")= "+str(zbs[countn][countm])+",\n"
        else:
            for countm in range(2*nmodes):
                text=text+"rbc("+str(countm-nmodes)+","+str(countn)+")= "+str(rbc[countn][countm])+", zbs("+str(countm-nmodes)+","+str(countn)+")= "+str(zbs[countn][countm])+",\n"
    replace(qvfilename+".sp","!----- Boundary Parameters -----",text)
    replace(qvfilename+".sp"," Nfp         =         4"," Nfp         =         "+str(nfp))
    replace(qvfilename+".sp"," phiedge     =   2.000000000000000E+00"," phiedge     =   "+str(PHIEDGE))
    #replace(qvfilename+".sp","pressure    =   0.000000000000000E+00","pressure    =   "+p2)

def runSPEC(name,executables_path,plotting_path,stel,r_edge, nmodes=25):
    print("Output to SPEC")
    qscfile="quasisymmetry_out."+name+".nc"
    output2spec(name,qscfile,executables_path,nmodes,stel,r_edge)
    print("Run SPEC")
    bashCommand = executables_path+"/./xspec "+name+".sp"
    run(bashCommand.split())
    print("Plot SPEC result")
    sys.path.insert(1, plotting_path)
    import specPlot
    specPlot.main(name+'.sp.h5',name)

## Run REGCOIL
def output2regcoil(regcoilFile,vmecFile,nescinfilename,coilSeparation,targetValue):
	with open(regcoilFile, "w") as txt_file:
		txt_file.write('&regcoil_nml\n')
		txt_file.write(' general_option=5\n')
		txt_file.write(' nlambda=100\n')
		txt_file.write(' lambda_min = 1e-19\n')
		txt_file.write(' lambda_max = 1e-12\n')
		txt_file.write(' target_option = "max_Bnormal"\n')
		txt_file.write(' target_value = '+str(targetValue)+'\n')
		txt_file.write(' ntheta_plasma=64\n')
		txt_file.write(' ntheta_coil  =64\n')
		txt_file.write(' nzeta_plasma =64\n')
		txt_file.write(' nzeta_coil   =64\n')
		txt_file.write(' mpol_potential = 8\n')
		txt_file.write(' ntor_potential = 8\n')
		txt_file.write(' geometry_option_plasma = 2\n')
		txt_file.write(" wout_filename = '"+vmecFile+"'\n")
		txt_file.write(' geometry_option_coil=2\n')
		txt_file.write(' separation = '+str(coilSeparation)+'\n')
		txt_file.write(" nescin_filename = '"+nescinfilename+"'\n")
		txt_file.write("/\n")

def runREGCOIL(name,stel,r_edge,executables_path,plotting_path,coilSeparation,targetValue,nCoilsPerNFP):
    print("Output to REGCOIL")
    nescinfilename = name+"_nescin.out"
    output2regcoil("regcoil_in."+name,"wout_"+name+".nc",nescinfilename,coilSeparation,targetValue)
    print("Run REGCOIL")
    bashCommand = executables_path+"/./regcoil regcoil_in."+name
    run(bashCommand.split())
    print("Cut coils from regcoil")
    bashCommand = plotting_path+"/./cutCoilsFromRegcoil regcoil_out."+name+".nc "+nescinfilename+" "+str(nCoilsPerNFP)+" 0 -1"
    run(bashCommand.split())
    print("Plot REGCOIL result")
    sys.path.insert(1, plotting_path)
    import REGCOILplot
    REGCOILplot.main(name, stel, r_edge)

def runSTAGE2(name, plotting_path, stel, r_edge, ncoils=16, R0=1.0, R1=0.6, order=5, ALPHA=1e-4, MIN_DIST=0, BETA=0, MAXITER=500, run_get_coils=True):
    if run_get_coils:
        import numpy as np
        from scipy.optimize import minimize
        from simsopt.geo.surfacerzfourier import SurfaceRZFourier
        from simsopt.objectives.fluxobjective import SquaredFlux, CoilOptObjective
        from simsopt.geo.curve import RotatedCurve, curves_to_vtk, create_equally_spaced_curves
        from simsopt.field.biotsavart import BiotSavart
        from simsopt.field.coil import Current, coils_via_symmetries, ScaledCurrent
        from simsopt.geo.curveobjectives import CurveLength, MinimumDistance
        filename = "input."+name
        # Initialize the boundary magnetic surface:
        nphi = 32
        ntheta = 32
        s = SurfaceRZFourier.from_vmec_input(filename, range="half period", nphi=nphi, ntheta=ntheta)
        # Create the initial coils:
        base_curves = create_equally_spaced_curves(ncoils, s.nfp, stellsym=True, R0=R0, R1=R1, order=order)
        # Since the target field is zero, one possible solution is just to set all
        # currents to 0. To avoid the minimizer finding that solution, we fix one
        # of the currents:
        # base_currents[0].fix_all()
        # base_currents = [Current(3.6e5) for i in range(ncoils)]
        base_currents = []
        for i in range(ncoils):
            curr = Current(1.)
            # since the target field is zero, one possible solution is just to set all
            # currents to 0. to avoid the minimizer finding that solution, we fix one
            # of the currents
            if i == 0:
                curr.fix_all()
            base_currents.append(ScaledCurrent(curr, 3.6e5))

        coils = coils_via_symmetries(base_curves, base_currents, s.nfp, True)
        bs = BiotSavart(coils)
        bs.set_points(s.gamma().reshape((-1, 3)))

        curves = [c.curve for c in coils]
        curves_to_vtk(curves, name+"_curves_init")
        pointData = {"B_N": np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)[:, :, None]}
        s.to_vtk(name+"_surf_init", extra_data=pointData)

        # Define the objective function:
        Jf = SquaredFlux(s, bs)
        Jls = [CurveLength(c) for c in base_curves]
        Jdist = MinimumDistance(curves, MIN_DIST)

        JF = CoilOptObjective(Jf, Jls, ALPHA, Jdist, BETA)

        # We don't have a general interface in SIMSOPT for optimisation problems that
        # are not in least-squares form, so we write a little wrapper function that we
        # pass directly to scipy.optimize.minimize
        def fun(dofs):
            JF.x = dofs
            J = JF.J()
            grad = JF.dJ()
            cl_string = ", ".join([f"{J.J():.3f}" for J in Jls])
            mean_AbsB = np.mean(bs.AbsB())
            jf = Jf.J()
            print(f"J={J:.3e}, Jflux={jf:.3e}, sqrt(Jflux)/Mean(|B|)={np.sqrt(jf)/mean_AbsB:.3e}, CoilLengths=[{cl_string}], ||∇J||={np.linalg.norm(grad):.3e}")
            return J, grad

        # print("""
        # ################################################################################
        # ### Perform a Taylor test ######################################################
        # ################################################################################
        # """)
        # f = fun
        dofs = JF.x
        # np.random.seed(1)
        # h = np.random.uniform(size=dofs.shape)
        # J0, dJ0 = f(dofs)
        # dJh = sum(dJ0 * h)
        # for eps in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
        #     J1, _ = f(dofs + eps*h)
        #     J2, _ = f(dofs - eps*h)
        #     print("err", (J1-J2)/(2*eps) - dJh)

        print("""
        ################################################################################
        ### Run the optimisation #######################################################
        ################################################################################
        """)
        res = minimize(fun, dofs, jac=True, method='L-BFGS-B', options={'maxiter': MAXITER, 'maxcor': 400}, tol=1e-15)
        curves_to_vtk(curves, name+"_curves_opt")
        pointData = {"B_N": np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)[:, :, None]}
        s.to_vtk(name+"_surf_opt", extra_data=pointData)

        # Output the coil curves to a txt file
        coilsFourier = np.array([coils[i].curve.get_dofs() for i in range(ncoils)])
        coilsOutput = []
        for j in range(ncoils):
            coilsOutput.append(np.insert([coilsFourier[j][2*(i+0*order)+1] for i in range(order)],0,0.0))
            coilsOutput.append([coilsFourier[j][2*(i+0*order)+0] for i in range(order+1)])
            coilsOutput.append(np.insert([coilsFourier[j][2*(i+1*order)+2] for i in range(order)],0,0.0))
            coilsOutput.append([coilsFourier[j][2*(i+1*order)+1] for i in range(order+1)])
            coilsOutput.append(np.insert([coilsFourier[j][2*(i+2*order)+3] for i in range(order)],0,0.0))
            coilsOutput.append([coilsFourier[j][2*(i+2*order)+2] for i in range(order+1)])
        coilsOutput = np.array(coilsOutput).transpose()
        np.savetxt(name+"_coil_curve_data.txt", coilsOutput, delimiter=',')

        # Check that it was outputted correctly
        coil_data = np.loadtxt(name+"_coil_curve_data.txt", delimiter=',')
        np.testing.assert_allclose(coilsOutput,coil_data)

        # Check that loading the data from the txt file yields the same coils
        from simsopt.geo.curvexyzfourier import CurveXYZFourier
        coil_data = CurveXYZFourier.load_curves_from_file(name+"_coil_curve_data.txt", order=order)
        coils_data = coils_via_symmetries(coil_data, base_currents, s.nfp, True)
        coilsFourier_data = np.array([coils_data[i].curve.get_dofs() for i in range(ncoils)])
        np.testing.assert_allclose(coilsFourier,coilsFourier_data)

        # Output the current data
        base_currents_Output = np.array([coils[i].current.get_value() for i in range(ncoils)])
        np.savetxt(name+"_coil_current_data.txt", base_currents_Output, delimiter=',')

        # Check that it was outputted correctly
        current_data = np.loadtxt(name+"_coil_current_data.txt", delimiter=',')
        np.testing.assert_allclose(base_currents_Output,current_data)

    print("Plot STAGE2 result")
    sys.path.insert(1, plotting_path)
    import STAGE2plot
    STAGE2plot.main(name, order, stel, r_edge)