#!/usr/bin/env python

from glob import glob
import numpy as np
from scipy.io import netcdf
from simsopt.geo.curve import curves_to_vtk
from simsopt.geo.curvexyzfourier import CurveXYZFourier
from simsopt.field.coil import Current, coils_via_symmetries
from simsopt.geo.surfacerzfourier import SurfaceRZFourier

#well = False
well = True

if well:
    prefix = 'output_well_True_lengthbound_24.0_kap_5.0_msc_5.0_dist_0.1_fil_0_ig_5_order_16_expquad_qfm_None_curve_'
    wout_filename = '/Users/mattland/Box Sync/work21/wout_20210728-01-010_QA_nfp2_A6_magwell_weight_1.00e+01_rel_step_3.00e-06_centered.nc'
else:
    prefix = 'output_well_False_lengthbound_24.0_kap_5.0_msc_5.0_dist_0.1_fil_0_ig_2_order_16_expquad_qfm_None_curve_'
    wout_filename = '/Users/mattland/Box Sync/work21/wout_20210704-01-050_weight_3.00e+01_rel_step_1.00e-05_forward_nfp2_QA.nc'
    
f = netcdf.netcdf_file(wout_filename, 'r', mmap=False)
phi = f.variables['phi'][()]
iotaf = f.variables['iotaf'][()]
presf = f.variables['presf'][()]
iotas = f.variables['iotas'][()]
pres = f.variables['pres'][()]
ns = f.variables['ns'][()]
nfp = f.variables['nfp'][()]
xn = f.variables['xn'][()]
xm = f.variables['xm'][()]
xn_nyq = f.variables['xn_nyq'][()]
xm_nyq = f.variables['xm_nyq'][()]
rmnc = f.variables['rmnc'][()]
zmns = f.variables['zmns'][()]
lmns = f.variables['lmns'][()]
bmnc = f.variables['bmnc'][()]

mpol = f.variables['mpol'][()]
ntor = f.variables['ntor'][()]
data = f.variables['aspect'][()]
print("aspect:            ",data)

f.close()
nmodes = len(xn)
print("bmnc.shape:",bmnc.shape)
iota = iotaf[-1]


######################################
# Initialize figure
######################################

ntheta = 60
nphi = 250
theta1D = np.linspace(0, 2 * np.pi, ntheta) + np.pi / 2
phi1D = np.linspace(0, 2 * np.pi, nphi)
phi2D, theta2D = np.meshgrid(phi1D, theta1D)
iradius = ns - 1
R = np.zeros((ntheta, nphi))
Z = np.zeros((ntheta, nphi))
B = np.zeros((ntheta, nphi))
for imode in range(nmodes):
    angle = xm[imode] * theta2D - xn[imode] * phi2D
    R = R + rmnc[iradius, imode] * np.cos(angle)
    Z = Z + zmns[iradius, imode] * np.sin(angle)

for imode in range(len(xn_nyq)):
    angle = xm_nyq[imode] * theta2D - xn_nyq[imode] * phi2D
    B = B + bmnc[iradius, imode] * np.cos(angle)

#shift = np.pi / 2 # qsc and simsopt configs are shifted by 1/2 a period
shift = 0
X = R * np.cos(phi2D + shift)
Y = R * np.sin(phi2D + shift)
# Rescale to lie in [0,1]:
B_rescaled = (B - B.min()) / (B.max() - B.min())

from mayavi import mlab
#fig = mlab.figure(bgcolor=(1,1,1), size=(300, 240))
#fig = mlab.figure(bgcolor=(1,1,1), size=(150, 130))
#fig = mlab.figure(bgcolor=(1,1,1), size=(225, 176))
#fig = mlab.figure(bgcolor=(1,1,1), size=(1200, 800))
fig = mlab.figure(bgcolor=(1,1,1), size=(225, 188))

#mlab.mesh(x_2D_plot_rotated, y_2D_plot_rotated-shift_array[i], z_2D_plot_rotated, scalars=Bmag, colormap='viridis')
# Flip sign of Z to match helicity of the qsc config.
mlab.mesh(X, Y, Z, scalars=B, colormap='viridis')


currents = [Current(1.0) for j in range(4)]

print("Handling configuration", prefix)

curves = []
for j in range(4):
    c = CurveXYZFourier(quadpoints=100, order=16)
    data = np.loadtxt(prefix + f"{j}.txt")
    c.x = data
    curves.append(c)

coils = coils_via_symmetries(curves, currents, nfp=2, stellsym=True)

colors = [(0.4, 0.4, 0.4),
        (0.6, 0.6, 0.6),
        (0.8, 0.8, 0.8),
        (1.0, 1.0, 1.0)]
for j in range(16):
    color = colors[j % 4]
    coils[j].curve.plot(engine="mayavi", close=True, show=False, color=color)

#[coil.curve.plot(engine="mayavi", close=True, show=False) for coil in coils]
#surf.plot(engine="mayavi", close=True)
#mlab.show()

nframes = 120
j = 0
#mlab.view(azimuth=360.0 * j / nframes, elevation=-70, distance=7, focalpoint=[0, 0, 0])
filenames = []
# For some reason, the first frame is darker, so repeat that frame at the end.
for j in range(nframes + 1):
    mlab.view(azimuth=180.0 * j / nframes, elevation=-70, distance=6., focalpoint=[0, 0, -0.30])
    filename = f'temp{j:04}.png'
    mlab.savefig(filename=filename)
    filenames.append(filename)

import imageio
# Duration is the number of seconds each frame is shown
with imageio.get_writer(__file__ + f'_well{int(well)}.gif', mode='I', duration=0.03) as writer:
    for filename in filenames[1:]:
        image = imageio.imread(filename)
        writer.append_data(image)