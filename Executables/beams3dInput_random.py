#!/usr/bin/env python3

"""
This version creates initial conditions that are randomly
distributed over the first field period of a surface and over
velocity-space.
"""

import sys
import os
import numpy as np
from scipy.io import netcdf
from scipy.interpolate import interp1d, RectBivariateSpline
import matplotlib.pyplot as plt

def main(name, s0 = 0.3, nparticles = 1000, ntheta = 60, nphi = 100):

    filename = "wout_"+name+".nc"
    justwout = os.path.basename(filename)

    outfile = 'particles.{}_volumeJacobian_s{}_n{}'.format(justwout[5:-3], s0, nparticles)
    beams3dfile = 'beams3d_in.{}_volumeJacobian_s{}_n{}'.format(justwout[5:-3], s0, nparticles)

    print('filename: ', filename)
    print('outfile: ', outfile)
    print('beams3dfile: ', beams3dfile)
    print('Requested normalized flux: ', s0)
    print('nparticles: ', nparticles)

    # Data for alpha particles:
    energy_eV = 3.5e6 # Alpha birth energy, in eV
    e_C =1.602176634e-19 # Charge of the electron, NOT of the alpha!
    m_kg = 6.644657230e-27 # Mass of the alpha, NOT of a proton!
    energy_J = energy_eV * e_C # Alpha birth energy, in Joules
    v = np.sqrt(2 * energy_J / m_kg) # Alpha birth speed [meter/second]
    print("v [m/s]:", v)

    f = netcdf.netcdf_file(filename,'r',mmap=False)

    nfp = f.variables['nfp'][()]
    ns = f.variables['ns'][()]
    xm = f.variables['xm'][()]
    xn = f.variables['xn'][()]
    xm_nyq = f.variables['xm_nyq'][()]
    xn_nyq = f.variables['xn_nyq'][()]
    rmnc = f.variables['rmnc'][()]
    zmns = f.variables['zmns'][()]
    gmnc = f.variables['gmnc'][()]
    bmnc = f.variables['bmnc'][()]

    # bmnc and gmnc are on the half grid and use the Nyquist modes.
    # rmnc and zmns are on the full grid and use the non-Nyquist modes.

    f.close()

    s_full = np.linspace(0, 1, ns)
    ds = s_full[1] - s_full[0]
    s_half = s_full[1:] - ds / 2
    print('rmnc.shape: ', rmnc.shape)

    interp_method='linear'

    index = ns - 1
    print('s_full[index]: ', s_full[index])
    old_rmnc = rmnc[index,:]
    print('rmnc[index,:]: ', old_rmnc)
    rc = interp1d(s_full, rmnc, axis=0, kind=interp_method)(s0)
    zs = interp1d(s_full, zmns, axis=0, kind=interp_method)(s0)
    print('New rmnc:', rmnc)

    print('difference:', np.max(np.abs(old_rmnc - rc)))

    # If s0 is close to 1, we may need to extrapolate off the end of the half grid.
    bc = interp1d(s_half, bmnc[1:, :], axis=0, kind=interp_method, fill_value='extrapolate')(s0)
    gc = interp1d(s_half, gmnc[1:, :], axis=0, kind=interp_method, fill_value='extrapolate')(s0)

    #theta1d = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    #phi1d = np.linspace(0, 2 * np.pi / nfp, nphi, endpoint=False)
    #dtheta = theta1d[1] - theta1d[0]
    #dphi = phi1d[1] - phi1d[0]

    theta1d = np.linspace(0, 2 * np.pi, ntheta, endpoint=True)
    phi1d = np.linspace(0, 2 * np.pi / nfp, nphi, endpoint=True)

    dtheta = theta1d[1] - theta1d[0]
    dphi = phi1d[1] - phi1d[0]

    phi, theta = np.meshgrid(phi1d, theta1d)
    modB = np.zeros((ntheta, nphi))
    sqrtg = np.zeros((ntheta, nphi))
    r = np.zeros((ntheta, nphi))
    x = np.zeros((ntheta, nphi))
    y = np.zeros((ntheta, nphi))
    z = np.zeros((ntheta, nphi))
    dxdtheta = np.zeros((ntheta, nphi))
    dydtheta = np.zeros((ntheta, nphi))
    dzdtheta = np.zeros((ntheta, nphi))
    dxdphi = np.zeros((ntheta, nphi))
    dydphi = np.zeros((ntheta, nphi))
    dzdphi = np.zeros((ntheta, nphi))
    sinphi = np.sin(phi)
    cosphi = np.cos(phi)
    for imn in range(len(xm_nyq)):
        m = xm_nyq[imn]
        n = xn_nyq[imn]
        angle = m * theta - n * phi
        cosangle = np.cos(angle)
        modB += bc[imn] * cosangle
        sqrtg += gc[imn] * cosangle

    for imn in range(len(xm)):
        m = xm[imn]
        n = xn[imn]
        angle = m * theta - n * phi
        sinangle = np.sin(angle)
        cosangle = np.cos(angle)
        rmnc = rc[imn]
        zmns = zs[imn]
        r += rmnc * cosangle
        x += rmnc * cosangle * cosphi
        y += rmnc * cosangle * sinphi
        z += zmns * sinangle

        dxdtheta += rmnc * (-m * sinangle) * cosphi
        dydtheta += rmnc * (-m * sinangle) * sinphi
        dzdtheta += zmns * m * cosangle

        dxdphi += rmnc * (n * sinangle * cosphi \
                        + cosangle * (-sinphi))
        dydphi += rmnc * (n * sinangle * sinphi \
                        + cosangle * cosphi)
        dzdphi += zmns * (-n * cosangle)
        if False:
            # Eventually we might include non-stellarator-symmetric cases:
            rmns = self.get_rs(m, n_without_nfp)
            zmnc = self.get_zc(m, n_without_nfp)
            r += rmns * sinangle
            x += rmns * sinangle * cosphi
            y += rmns * sinangle * sinphi
            z += zmnc * cosangle

            dxdtheta += rmns * (m * cosangle) * cosphi
            dydtheta += rmns * (m * cosangle) * sinphi
            dzdtheta += zmnc * (-m * sinangle)
                        
            dxdphi += rmns * (-n * cosangle * cosphi \
                            + sinangle * (-sinphi))
            dydphi += rmns * (-n * cosangle * sinphi \
                            + sinangle * cosphi)
            dzdphi += zmnc * (n * sinangle)

    normalx = dydphi * dzdtheta - dzdphi * dydtheta
    normaly = dzdphi * dxdtheta - dxdphi * dzdtheta
    normalz = dxdphi * dydtheta - dydphi * dxdtheta
    norm_normal = np.sqrt(normalx * normalx + normaly * normaly + normalz * normalz)
    area = np.sum(norm_normal[:-1, :-1]) * dtheta * dphi * nfp # Note we must drop repeated grid points in norm_normal
    print("Computed area of the surface: ", area)

    # Use |sqrtg| instead of sqrtg from here onward:
    sqrtg = np.abs(sqrtg)
    print('Before ** 4, min(sqrtg)={}, max(sqrtg)={}'.format(np.min(sqrtg), np.max(sqrtg)))
    #sqrtg = sqrtg ** 4
    print('After ** 4, min(sqrtg)={}, max(sqrtg)={}'.format(np.min(sqrtg), np.max(sqrtg)))

    # Initialize 2D splines
    r_spl = RectBivariateSpline(theta1d, phi1d, r)
    z_spl = RectBivariateSpline(theta1d, phi1d, z)
    modB_spl = RectBivariateSpline(theta1d, phi1d, modB)
    sqrtg_spl = RectBivariateSpline(theta1d, phi1d, sqrtg)
    norm_normal_spl = RectBivariateSpline(theta1d, phi1d, norm_normal)

    thetas = np.zeros(nparticles)
    phis = np.zeros(nparticles)
    #rs = np.zeros(nparticles)
    #zs = np.zeros(nparticles)
    #vlls = np.zeros(nparticles)
    #mus = np.zeros(nparticles)

    #rng = np.random.default_rng()
    rng = np.random
    max_norm_normal = np.max(norm_normal)
    max_sqrtg = np.max(sqrtg)
    for j in range(nparticles):
        print('j={}: '.format(j), end='')
        while True:
            print('.', end='')
            theta0 = rng.random() * 2 * np.pi
            phi0 = rng.random() * 2 * np.pi / nfp
            #f = rng.random() * max_norm_normal
            #if f <= norm_normal_spl(theta0, phi0):
            #    break
            f = rng.random() * max_sqrtg
            if f <= sqrtg_spl(theta0, phi0):
                break
        thetas[j] = theta0
        phis[j] = phi0
        print()

    rs = r_spl.ev(thetas, phis)
    zs = z_spl.ev(thetas, phis)

    # Random numbers on [-1, 1]:
    xlls = rng.random((nparticles,)) * 2 - 1
    vlls = xlls * v
    # In BEAMS3D, "MU" is defined as 0.5 * m * vperp^2 / B
    mus = 0.5 * m_kg * (1 - xlls * xlls) * v * v / modB_spl.ev(thetas, phis)

    f = open(outfile, 'w')
    f.write('# nparticles\n')
    f.write('{}\n'.format(nparticles))
    f.write('# theta, phi, v||, R, Z, mu, weight\n')
    fmt_master = '{:24.15e}'
    fmt_left = '{:<24.15e}'
    fmt = 7 * fmt_master
    fmt = fmt + '\n'
    for j in range(nparticles):
        f.write(fmt.format(thetas[j], phis[j], vlls[j], rs[j], zs[j], mus[j], 1.0 / nparticles))
    f.close()

    # Now write the info needed by BEAMS3D:
    f = open(beams3dfile, 'w')

    numstr = str(nparticles) + '*'

    f.write('T_END_IN = ' + numstr + '1.0d-5\n')
    f.write('ZATOM_IN = ' + numstr + '2.0d+0\n')
    f.write('CHARGE_IN = ' + numstr + fmt_left.format(2 * e_C) + '\n')
    f.write('MASS_IN = ' + numstr + fmt_left.format(m_kg) + '\n')
    f.write('! Using volume Jacobian\n')
    mystr = 'R_START_IN = '
    for j in range(nparticles):
        mystr += fmt_master.format(rs[j])
        if np.mod(j, 5) == 4:
            mystr += '\n'
        f.write(mystr)
        mystr = ' '
    f.write('\n')

    mystr = 'PHI_START_IN = '
    for j in range(nparticles):
        mystr += fmt_master.format(phis[j])
        if np.mod(j, 5) == 4:
            mystr += '\n'
        f.write(mystr)
        mystr = ' '
    f.write('\n')

    mystr = 'Z_START_IN = '
    for j in range(nparticles):
        mystr += fmt_master.format(zs[j])
        if np.mod(j, 5) == 4:
            mystr += '\n'
        f.write(mystr)
        mystr = ' '
    f.write('\n')

    mystr = 'VLL_START_IN = '
    for j in range(nparticles):
        mystr += fmt_master.format(vlls[j])
        if np.mod(j, 5) == 4:
            mystr += '\n'
        f.write(mystr)
        mystr = ' '
    f.write('\n')

    mystr = 'MU_START_IN = '
    for j in range(nparticles):
        mystr += fmt_master.format(mus[j])
        if np.mod(j, 5) == 4:
            mystr += '\n'
        f.write(mystr)
        mystr = ' '
    f.write('\n')

    f.close()

    print('min(sqrtg)={}, max(sqrtg)={}'.format(np.min(sqrtg), np.max(sqrtg)))

    #######################################
    # Make figure of the points
    #######################################

    exit(0)

    fig = plt.figure()
    plt.contourf(phi, theta, sqrtg, 25)
    plt.colorbar()
    plt.xlabel('phi')
    plt.ylabel('theta')
    plt.plot(phis, thetas, '.k', ms=2)
    plt.title('Color = sqrtg')

    plt.show()
