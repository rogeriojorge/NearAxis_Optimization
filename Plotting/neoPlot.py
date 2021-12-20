#!/usr/bin/env python3
import sys

def main(file,qvfilename):

    import matplotlib.pyplot as plt
    import numpy as np

    token = open('neo_out.'+qvfilename,'r')
    linestoken=token.readlines()
    eps_eff=[]
    s_radial=[]
    for x in linestoken:
        s_radial.append(float(x.split()[0])/150)
        eps_eff.append(float(x.split()[1])**(2/3))
    token.close()
    s_radial = np.array(s_radial)
    eps_eff = np.array(eps_eff)
    s_radial = s_radial[np.argwhere(~np.isnan(eps_eff))[:,0]]
    eps_eff = eps_eff[np.argwhere(~np.isnan(eps_eff))[:,0]]
    fig = plt.figure(figsize=(7, 3), dpi=200)
    ax = fig.add_subplot(111)
    plt.plot(s_radial,eps_eff, label='eps eff '+qvfilename)
    # ax.set_yscale('log')
    plt.xlabel(r'$s=\psi/\psi_b$', fontsize=12)
    plt.ylabel(r'$\epsilon_{eff}$', fontsize=14)

    from scipy.optimize import curve_fit
    def objective_quadratic(x, a, b, c):
        return a * x*x + b * x + c
    popt_quadratic, _ = curve_fit(objective_quadratic, s_radial, eps_eff)
    a_quadratic, b_quadratic, c_quadratic = popt_quadratic
    y_new_quadratic = objective_quadratic(s_radial, a_quadratic, b_quadratic, c_quadratic)
    # plt.plot(s_radial,y_new_quadratic, label='quadratic fit'+str(popt_quadratic))

    # plt.legend()
    plt.tight_layout()
    fig.savefig('neo_out_'+qvfilename+'.pdf', dpi=fig.dpi)#, bbox_inches = 'tight', pad_inches = 0)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[1])