#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import py_spec

def main(file,qvfilename):
    s = py_spec.SPECout(file)
    s.plot_iota()
    plt.savefig('SPECiota'+qvfilename+'.pdf', bbox_inches = 'tight', pad_inches = 0)
    plt.close()
    fig = plt.figure(figsize=(14,7))
    # indPoincare=[0,12,24,36]
    # indPoincare=[0,24,48,72]
    indPoincare_array=[[0,12,24,36],[48,60,72,84]]
    for count, indPoincare in enumerate(indPoincare_array):
        fig = plt.figure(figsize=(14,7))
        fig.patch.set_facecolor('white')
        for toroidalIdx in indPoincare:
            rr = s.poincare.R[:, :, toroidalIdx]
            zz = s.poincare.Z[:, :, toroidalIdx]
            nptrj = rr.shape[0]
            for ii in range(nptrj):
                plt.plot(rr[ii, :],zz[ii, :],'o',ms=1.0,markeredgecolor=None,mew=0)
        print("Save PDF "+str(count))
        plt.xlabel('R (meters)', fontsize=14)
        plt.ylabel('Z (meters)', fontsize=14)
        plt.gca().set_aspect('equal',adjustable='box')
        plt.savefig('SPECplot_'+str(count)+'_'+qvfilename+'.pdf', bbox_inches = 'tight', pad_inches = 0.1)
        plt.savefig('SPECplot_'+str(count)+'_'+qvfilename+'.png', bbox_inches = 'tight', pad_inches = 0.1)
