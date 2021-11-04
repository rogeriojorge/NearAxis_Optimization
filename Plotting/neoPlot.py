#!/usr/bin/env python3

def main(file,qvfilename):

    import matplotlib.pyplot as plt
    import numpy as np

    filename = file

    token = open('neo_out.'+qvfilename,'r')
    linestoken=token.readlines()
    column_number = 1
    eps_eff=[]
    s_radial=[]
    for x in linestoken:
        s_radial.append(float(x.split()[0])/150)
        eps_eff.append(float(x.split()[1]))
    token.close()
    s_radial = np.array(s_radial)
    eps_eff = np.array(eps_eff)
    s_radial = s_radial[np.argwhere(~np.isnan(eps_eff))[:,0]]
    eps_eff = eps_eff[np.argwhere(~np.isnan(eps_eff))[:,0]]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(s_radial,eps_eff, label='eps eff '+qvfilename)
    # ax.set_yscale('log')
    plt.xlabel(r'$s=\psi/\psi_a$', fontsize=12)
    plt.ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=12)

    from scipy.optimize import curve_fit
    def objective(x, a, b):
        return a * x + b
    popt, _ = curve_fit(objective, s_radial, eps_eff)
    a, b = popt
    y_new = objective(s_radial, a, b)
    plt.plot(s_radial,y_new, label='linear fit'+str(popt))

    plt.legend()
    
    fig.savefig('neo_out_'+qvfilename+'.pdf', dpi=fig.dpi)