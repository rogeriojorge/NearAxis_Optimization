from qsc import Qsc
from simsopt import make_optimizable
import numpy as np

def get_stel(ind,nphi=251):
    if ind==0:
        name   = 'QA_NFP2_vac'
        rc     = [ 1.0,-0.16269301331725952,0.018893800182906907,-0.0020843525869404015,0.0002032623426799173,-1.859988194152255e-05,1.7428172152882315e-06 ]
        zs     = [ 0.0,0.17341226970640833,-0.019937729110505058,0.0021985139740552787,-0.00022321977387611668,2.2615738120388246e-05,-2.0970137122833628e-06 ]
        etabar =  0.7699121990257152
        nfp    =  2
        B2c    =  0.4551751876618826
        p2     =  0.0
        iota   =  0.40376036676890553
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2))
        r_edge =  0.1
        coilSeparation = 0.2
        targetValue = 0.005
        nCoilsPerNFP = 6
    if ind==1:
        name   = 'QA_NFP2_plasma'
        rc     = [ 1.0,-0.15178890701086753,0.011580426258003833,-0.00012387677229852621,-2.9080534353763382e-05,7.70146723929258e-08 ]
        zs     = [ 0.0,0.1544147689230838,-0.011000741042765286,-0.00013865026415631695,8.990615381841263e-05,1.3076196767557744e-07 ]
        etabar =  0.6501173908230352
        nfp    =  2
        B2c    =  -0.915195507869042
        p2     =  -500000.0
        iota   =  0.39160446975249485
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2))
        r_edge =  0.1
        coilSeparation = 0.18
        targetValue = 0.03
        nCoilsPerNFP = 6
    if ind==2:
        name   = 'QA_NFP3_vac'
        rc     = [ 1.0,0.0455013255637958,0.0021091197706754297,-9.489832223804648e-05,-9.601800003988155e-06,2.806985889865061e-06 ]
        zs     = [ 0.0,-0.055767062771199534,-0.002235860163818757,8.771253863082454e-05,1.2823168297031853e-05,-2.43079218142775e-06 ]
        etabar =  0.7811463903776809
        nfp    =  3
        B2c    =  0.6726004478276207
        p2     =  0.0
        iota   =  0.3999857272555072
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2))
        r_edge =  0.1
        coilSeparation = 0.18
        targetValue = 0.03
        nCoilsPerNFP = 5
    if ind==3:
        name   = 'QA_NFP3_plasma'
        rc     = [ 1.0,0.049617467646414413,0.0020501552448295222,-7.952960068630817e-05,-4.1525851418898675e-06 ]
        zs     = [ 0.0,-0.05618080515173391,-0.0022200323317885145,8.348675846023883e-05,1.3207647066584537e-05 ]
        etabar =  0.6741477835089387
        nfp    =  3
        B2c    =  -1.1823423764823269
        p2     =  -500000.0
        iota   =  0.4183742218322489
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2))
        r_edge =  0.1
        coilSeparation = 0.17
        targetValue = 0.05
        nCoilsPerNFP = 5
    if ind==4:
        name   = 'QA_NFP4_vac'
        rc     = [ 1.0,0.01958274954753399,0.0018243506698573202,0.00010965039865412234,4.73894722623086e-06,-1.5461083677051984e-06 ]
        zs     = [ 0.0,-0.024358455876563764,-0.0006678125089242426,-0.00013757861553785294,-1.6184908843018587e-05,-3.126955879553316e-06 ]
        etabar =  0.9390078410893093
        nfp    =  4
        B2c    =  0.055309698986764286
        p2     =  0.0
        iota   =  0.3369558032942888
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2))
        r_edge =  0.05
        coilSeparation = 0.16
        targetValue = 0.05
        nCoilsPerNFP = 4
    if ind==5:
        name   = 'QH_NFP4_vac'
        rc     = [ 1.0,0.3567120764364836,0.09309473742665696,0.02068580554650727,0.003352763602115162,0.00029590481856858514 ]
        zs     = [ 0.0,0.30605344736274837,0.08535505841371785,0.020300831381059183,0.003604828690431057,0.0003636106532187047 ]
        etabar =  1.7131771845862305
        nfp    =  4
        B2c    =  0.5528452930602289
        p2     =  0.0
        iota   =  -1.638024385998008
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2))
        r_edge =  0.1
        coilSeparation = 0.1
        targetValue = 0.25
        nCoilsPerNFP = 4
    if ind==6:
        name    = 'QI_NFP1_vac'
        rc     = [ 1.0,0.0,-0.2,0.0,0.0,0.0,0.0,0.0,0.0 ]
        zs     = [ 0.0,0.0,0.3538620388942289,0.0,-0.03190697863469498,0.0,-0.019442780647161585,0.0,0.005700818220902032 ]
        etabar =  0.0
        nfp    =  1
        B2c    =  0.0
        p2     =  0.0
        iota   =  0.9866126874975886
        B0_vals = [ 1.0,-0.003759072330629602 ]
        d_svals = [ 0.0,1.0594127385850904,3.65993529337379e-05,0.3892573325353277,0.00029829921421492366,-0.026720277634932094,0.00026679451432573466,-0.25437872550520174 ]
        alpha0  = -4.712534887915769
        c0      = -1.5707963267948966
        delta   = 0.2
        stel   =  make_optimizable(Qsc(rc=rc,zs=zs, nfp=1, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi,phi_shift=1/3, omn=True, delta=delta, alpha0=alpha0, c0=c0))
        r_edge =  0.1
        coilSeparation = 0.05
        targetValue = 0.25
        nCoilsPerNFP = 4
    if ind==7:
        name   = 'QA_NFP2_plasma_2'
        rc     = [ 1.0,-0.1944554440057859,0.024896649310590737,-0.00324347056650685,0.0003957790550103649,-3.5124028349806346e-05 ]
        zs     = [ 0.0,0.16947886640795443,-0.024367760721762312,0.003252403175044678,-0.00039710651838295813,3.605237502341612e-05 ]
        etabar =  0.5910044289478692
        nfp    =  2
        B2c    =  -1.3148495687068937
        p2     =  -500000.0
        iota   =  0.4377598860159437
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2))
        r_edge =  0.08
        coilSeparation = 0.2
        targetValue = 0.14
        nCoilsPerNFP = 6
    return stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP