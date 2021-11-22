#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv

def main(name, stel, r_edge, savefig=True):
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

    import mayavi.mlab as mlab
    from matplotlib.pyplot import cm
    fig = mlab.figure(bgcolor=(1,1,1), size=(430,720))
    for count, coil in enumerate(coilPos):
        mlab.plot3d(np.array(coil).transpose()[0], np.array(coil).transpose()[1], np.array(coil).transpose()[2], tube_radius=0.02, color=(1,0,0))

    ntheta=40
    nphi=130
    X_qsc, Y_qsc, Z_qsc, R_qsc = stel.get_boundary(r=r_edge, ntheta=ntheta, nphi=nphi)

    theta1D = np.linspace(0, 2 * np.pi, ntheta)
    phi1D = np.linspace(0, 2 * np.pi, nphi)
    phi2D, theta2D = np.meshgrid(phi1D, theta1D)
    Bmag = stel.B_mag(r_edge, theta2D, phi2D)

    mlab.mesh(X_qsc, Y_qsc, Z_qsc, scalars=Bmag, colormap='viridis')
    mlab.view(azimuth=0, elevation=0, distance=8.5, focalpoint=(-0.15,0,0), figure=fig)
    # Create the colorbar and change its properties
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
        mlab.savefig(filename=name+'_3D_REGCOIL_coils.png', figure=fig)

    mlab.show()
