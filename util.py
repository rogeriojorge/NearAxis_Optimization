from numpy import savetxt
from subprocess import run
import sys
import shutil
from scipy.io import netcdf
from shutil import copyfile
import booz_xform as bx
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
    # ictrl[0] = restart_flag + readin_flag# + timestep_flag + output_flag + cleanup_flag
    # print("Calling runvmec. ictrl={} comm={}".format(ictrl, fcomm))
    # vmec.runvmec(ictrl, "input."+name, verbose, fcomm, reset_file)

    # np.testing.assert_equal(ictrl[1], 0)

    bashCommand = executables_path+"/./xvmec2000 input."+name
    run(bashCommand.split())
    print("Plot VMEC result")
    sys.path.insert(1, plotting_path)
    import vmecPlot2
    vmecPlot2.main("wout_"+name+".nc")

# Run booz_xform
def runBOOZXFORM(name,executables_path,plotting_path):
    print("Run BOOZ_XFORM")
    b1 = bx.Booz_xform()
    b1.read_wout("wout_"+name+".nc")
    b1.compute_surfs = [1,2,5,10,18,25,32,40,50,60,70,80,90,100,110,125,150,170,190,210,230,249]
    b1.run()
    print("Plot BOOZ_XFORM")
    bx.surfplot(b1, js=5,  fill=False, ncontours=35)
    plt.savefig("Boozxform_surfplot_1_"+name+'.pdf', bbox_inches = 'tight', pad_inches = 0);plt.close()
    bx.surfplot(b1, js=15, fill=False, ncontours=35)
    plt.savefig("Boozxform_surfplot_2_"+name+'.pdf', bbox_inches = 'tight', pad_inches = 0);plt.close()
    bx.surfplot(b1, js=21, fill=False, ncontours=35)
    plt.savefig("Boozxform_surfplot_3_"+name+'.pdf', bbox_inches = 'tight', pad_inches = 0);plt.close()
    if name[0:2] == 'QH':
        helical_detail = True
    else:
        helical_detail = False
    bx.symplot(b1, helical_detail = helical_detail, sqrts=True)
    plt.savefig("Boozxform_symplot_"+name+'.pdf', bbox_inches = 'tight', pad_inches = 0); plt.close()
    bx.modeplot(b1, sqrts=True)
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
def output2spec(qvfilename,qscfile,executables_path,nmodes):
    f = netcdf.netcdf_file(qscfile,mode='r',mmap=False)
    RBC = f.variables['RBC'][()]
    ZBS = f.variables['ZBS'][()]
    rc = f.variables['R0c'][()]
    zs = f.variables['Z0s'][()]
    nfp = f.variables['nfp'][()]
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
    text="! Axis shape\n"
    # text=text+"rac(0:"+str(len(rc)-1)+")="+", ".join([str(elem) for elem in rc])+"\n"
    # text=text+"zas(0:"+str(len(zs)-1)+")="+", ".join([str(elem) for elem in zs])+"\n"
    text = text+"Rac   =   "+" ".join([str(elem) for elem in rc])+"\n"
    text = text+"Zas   =   "+" ".join([str(elem) for elem in zs])+"\n"
    copyfile(executables_path+"/input.SPECexample.sp", qvfilename+".sp")
    replace(qvfilename+".sp","! Axis shape",text)
    text="!----- Boundary Parameters -----\n"
    for countn in range(len(rbc)-1):
        if countn==0:
            for countm in range(nmodes):
                text=text+"Rbc("+str(countm)+","+str(countn)+")= "+str(rbc[countn][countm])+" Zbs("+str(countm)+","+str(countn)+")= "+str(zbs[countn][countm])+",\n"
        elif countn==len(rbc)-2:
            for countm in range(2*nmodes):
                if countm==2*nmodes-1:
                    text=text+"Rbc("+str(countm-nmodes)+","+str(countn)+")= "+str(rbc[countn][countm])+" Zbs("+str(countm-nmodes)+","+str(countn)+")= "+str(zbs[countn][countm])+"\n"
                else:
                    text=text+"Rbc("+str(countm-nmodes)+","+str(countn)+")= "+str(rbc[countn][countm])+" Zbs("+str(countm-nmodes)+","+str(countn)+")= "+str(zbs[countn][countm])+",\n"
        else:
            for countm in range(2*nmodes):
                text=text+"Rbc("+str(countm-nmodes)+","+str(countn)+")= "+str(rbc[countn][countm])+" Zbs("+str(countm-nmodes)+","+str(countn)+")= "+str(zbs[countn][countm])+",\n"
    replace(qvfilename+".sp","!----- Boundary Parameters -----",text)
    replace(qvfilename+".sp"," Nfp         =         4"," Nfp         =         "+str(nfp))
    replace(qvfilename+".sp"," phiedge     =   2.000000000000000E+00"," phiedge     =   "+str(PHIEDGE))
    #replace(qvfilename+".sp","pressure    =   0.000000000000000E+00","pressure    =   "+p2)

def runSPEC(name,executables_path,plotting_path, nmodes=25):
    print("Output to SPEC")
    qscfile="quasisymmetry_out."+name+".nc"
    output2spec(name,qscfile,executables_path,nmodes)
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

def runREGCOIL(name,executables_path,plotting_path,coilSeparation,targetValue,nCoilsPerNFP):
	print("Output to REGCOIL")
	nescinfilename = name+"_nescin.out"
	output2regcoil("regcoil_in."+name,"wout_"+name+".nc",nescinfilename,coilSeparation,targetValue)
	print("Run REGCOIL")
	bashCommand = executables_path+"/./regcoil regcoil_in."+name
	run(bashCommand.split())
	print("Cut coils from regcoil")
	bashCommand = plotting_path+"/./cutCoilsFromRegcoil regcoil_out."+name+".nc "+nescinfilename+" "+str(nCoilsPerNFP)+" 0 -1"
	run(bashCommand.split())