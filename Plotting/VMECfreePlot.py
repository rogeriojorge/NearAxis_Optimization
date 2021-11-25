#!/usr/bin/env python3

import numpy as np
from scipy.io import netcdf
from pandas import read_csv

def main(name, savefig=True):

    # Import coils
    allCoils = read_csv("coils."+name)
    allCoilsValues = allCoils.values
    coilN=0
    coilPosN=0
    coilPos=[[],[]]
    for nVals in range(len(allCoilsValues)-2):
        listVals=allCoilsValues[nVals+2][0]
        vals=listVals.split()
        try:
            floatVals = [float(nVals) for nVals in vals][0:3]
            coilPos[coilN].append(floatVals)
            coilPosN=coilPosN+1
        except:
            try:
                floatVals = [float(nVals) for nVals in vals[0:3]][0:3]
                coilPos[coilN].append(floatVals)
                coilN=coilN+1
                coilPos.append([])
            except:
                break
    current=allCoilsValues[6][0].split()
    current=float(current[3])
    coilPos=np.array(coilPos[:-2])

    # Import VMEC
    filename = "wout_"+name+"_free.nc"
    f = netcdf.netcdf_file(filename,'r',mmap=False)
    ns = f.variables['ns'][()]
    nfp = f.variables['nfp'][()]
    xn = f.variables['xn'][()]
    xm = f.variables['xm'][()]
    xn_nyq = f.variables['xn_nyq'][()]
    xm_nyq = f.variables['xm_nyq'][()]
    rmnc = f.variables['rmnc'][()]
    zmns = f.variables['zmns'][()]
    bmnc = f.variables['bmnc'][()]
    lasym = f.variables['lasym__logical__'][()]
    if lasym==1:
        rmns = f.variables['rmns'][()]
        zmnc = f.variables['zmnc'][()]
        bmns = f.variables['bmns'][()]
    else:
        rmns = 0*rmnc
        zmnc = 0*rmnc
        bmns = 0*bmnc
    f.close()
    nmodes = len(xn)

    ntheta = 80
    nzeta = int(150*nfp)
    theta1D = np.linspace(0,2*np.pi,num=ntheta)
    zeta1D = np.linspace(0,2*np.pi,num=nzeta)
    zeta2D, theta2D = np.meshgrid(zeta1D,theta1D)
    iradius = ns-1
    R = np.zeros((ntheta,nzeta))
    Z = np.zeros((ntheta,nzeta))
    B = np.zeros((ntheta,nzeta))
    for imode in range(nmodes):
        angle = xm[imode]*theta2D - xn[imode]*zeta2D
        R = R + rmnc[iradius,imode]*np.cos(angle) + rmns[iradius,imode]*np.sin(angle)
        Z = Z + zmns[iradius,imode]*np.sin(angle) + zmnc[iradius,imode]*np.cos(angle)

    for imode in range(len(xn_nyq)):
        angle = xm_nyq[imode]*theta2D - xn_nyq[imode]*zeta2D
        B = B + bmnc[iradius,imode]*np.cos(angle) + bmns[iradius,imode]*np.sin(angle)

    X = R * np.cos(zeta2D)
    Y = R * np.sin(zeta2D)
    # Rescale to lie in [0,1]:
    # B_rescaled = (B - B.min()) / (B.max() - B.min())

    import mayavi.mlab as mlab
    from matplotlib.pyplot import cm
    fig = mlab.figure(bgcolor=(1,1,1), size=(430,720))
    # Plot coils
    for count, coil in enumerate(coilPos):
        mlab.plot3d(np.array(coil).transpose()[0], np.array(coil).transpose()[1], np.array(coil).transpose()[2], tube_radius=0.02, color=(1,0,0))

    # Plot VMEC
    mlab.mesh(X, Y, Z, scalars=B, colormap='viridis')
    mlab.view(azimuth=0, elevation=0, distance=8.5, focalpoint=(-0.15,0,0), figure=fig)

    cb = mlab.colorbar(orientation='horizontal', title='|B| [T]', nb_labels=7)
    cb.scalar_bar_representation.position = [0.1, 0.9]
    cb.scalar_bar_representation.position2 = [0.8, 0.05]
    cb.scalar_bar.unconstrained_font_size = True
    cb.label_text_property.font_family = 'times'
    cb.label_text_property.bold = 0
    cb.label_text_property.font_size=18
    cb.label_text_property.color=(0,0,0)
    cb.title_text_property.font_family = 'times'
    cb.title_text_property.font_size=30
    cb.title_text_property.color=(0,0,0)
    cb.title_text_property.bold = 1

    if savefig != None:
        mlab.savefig(filename=name+'_3D_VMEC_free_boundary.png', figure=fig)

    mlab.show()