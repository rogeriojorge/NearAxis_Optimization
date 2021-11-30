#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import h5py   

def main(name):
    beams3d_file = "beams3d_"+name+".h5"
    with h5py.File(beams3d_file, "r") as f:
        ## List all groups
        # print("Keys: %s" % f.keys())
        # a_group_key = list(f.keys())[0]
        # print(a_group_key)
        ## Get the data
        R_lines   = np.array(f['R_lines'])
        Z_lines   = np.array(f['Z_lines'])
        PHI_lines = np.array(f['PHI_lines'])
        NPOINC    = f['npoinc'][()][0]
        t_end     = f['t_end'][()][0]
        nsteps    = f['nsteps'][()][0]
        vll_lines = f['vll_lines'][()][0]
        
    X_lines = R_lines * np.cos(PHI_lines)
    Y_lines = R_lines * np.sin(PHI_lines)

    time    = np.linspace(0,t_end,NPOINC+1,endpoint=True)

    nParticles_plot = 1000

    print(R_lines[1])
    print(Z_lines[1])
    print(PHI_lines[1])
    # print(NPOINC)
    # print(t_end)
    # print(nsteps)
    # print(vll_lines)

    # import mayavi.mlab as mlab
    # fig = mlab.figure(bgcolor=(1,1,1), size=(430,720))
    # # for i in range(len(X_lines)):
    # for i in range(nParticles_plot):
    #     mlab.plot3d(X_lines[i][X_lines[i] != 0], Y_lines[i][Y_lines[i] != 0], Z_lines[i][Z_lines[i] != 0], color=(0.5,0.5,0.5), tube_radius=0.005)
    # mlab.show()
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # for i in range(nParticles_plot):
    #     plt.plot(X_lines[i][X_lines[i] != 0], Y_lines[i][Y_lines[i] != 0], Z_lines[i][Z_lines[i] != 0])

    # loss_over_time = np.array([len(R_lines[:,i][R_lines[:,i] == 0])/len(R_lines[:,i]) for i in range(len(R_lines[0,:]))])
    # plt.plot(time,loss_over_time)
    # # plt.xscale('log') 
    # plt.show()
