#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main(name, order, stel, r_edge, savefig=True):
    from simsopt.field.coil import Current, coils_via_symmetries
    from simsopt.geo.curvexyzfourier import CurveXYZFourier
    from simsopt.geo.surfacerzfourier import SurfaceRZFourier

    filename = "input."+name
    # Initialize the boundary magnetic surface:
    nphi = 32
    ntheta = 32
    s = SurfaceRZFourier.from_vmec_input(filename, range="full torus", nphi=nphi, ntheta=ntheta)

    current_data = np.loadtxt(name+"_coil_current_data.txt", delimiter=',')
    base_currents = [Current(currents) for currents in current_data]
    coil_data = CurveXYZFourier.load_curves_from_file(name+"_coil_curve_data.txt", order=order)
    coils = coils_via_symmetries(coil_data, base_currents, s.nfp, True)
    ncoils = len(coils)
    nbasecoils = len(base_currents)

    import mayavi.mlab as mlab
    from matplotlib.pyplot import cm
    # mlab.options.offscreen = True

    fig = mlab.figure(bgcolor=(1,1,1), size=(430,720))
    colors = tuple(map(tuple,np.delete(cm.rainbow(np.linspace(0, 1, int(ncoils/nbasecoils))),3,1)))
    i=-1
    for count, coil in enumerate(coils):
        if count % nbasecoils == 0:
            i += 1
        mlab.plot3d(coil.curve.gamma()[:, 0], coil.curve.gamma()[:, 1], coil.curve.gamma()[:, 2], color=colors[i], tube_radius=0.02)

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
        mlab.savefig(filename=name+'_3D_coils.png', figure=fig)

    mlab.show()