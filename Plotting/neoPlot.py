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

    # def objective_linear(x, a, b):
    #     return a * x + b
    # popt_linear, _ = curve_fit(objective_linear, s_radial, eps_eff)
    # a_linear, b_linear = popt_linear
    # y_new_linear = objective_linear(s_radial, a_linear, b_linear)
    # plt.plot(s_radial,y_new_linear, label='linear fit'+str(popt_linear))

    def objective_quadratic(x, a, b, c):
        return a * x*x + b * x + c
    popt_quadratic, _ = curve_fit(objective_quadratic, s_radial, eps_eff)
    a_quadratic, b_quadratic, c_quadratic = popt_quadratic
    y_new_quadratic = objective_quadratic(s_radial, a_quadratic, b_quadratic, c_quadratic)
    plt.plot(s_radial,y_new_quadratic, label='quadratic fit'+str(popt_quadratic))

    plt.legend()
    
    fig.savefig('neo_out_'+qvfilename+'.pdf', dpi=fig.dpi)