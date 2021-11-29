#!/usr/bin/env python3

import numpy as np
import mayavi.mlab as mlab

def main(name, stel, r_edge, show=True, savefig=True):
    fig = mlab.figure(bgcolor=(1,1,1), size=(430,720))

    ntheta=80
    nphi=200
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
    cb.label_text_property.font_size=20
    cb.label_text_property.color=(0,0,0)
    cb.title_text_property.font_family = 'times'
    cb.title_text_property.font_size=20
    cb.title_text_property.color=(0,0,0)
    cb.title_text_property.bold = 1

    if savefig != None:
        mlab.savefig(filename=name+'_simple_3Dplot.png', figure=fig)

    if show:
        mlab.show()
