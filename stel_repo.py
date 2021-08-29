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
        r_edge =  0.1
        coilSeparation = 0.31
        targetValue = 0.02
        nCoilsPerNFP = 6
    if ind==8:
        name   = 'QH_NFP4_vac_2_wo_gradgradB'
        rc     = [ 1.0,0.24166255996559133,0.03507927831848278,0.003072101422988407,-0.00021640475768656763,-0.0001363835654519441,-2.223736161913154e-05,-1.1792425120563113e-06 ]
        zs     = [ 0.0,0.21525825259832043,0.03221825938601267,0.0029388783262989686,-0.00018327610234595618,-0.000131613504217177,-2.3862394771575614e-05,-1.4509761624798847e-06 ]
        etabar =  1.9697675913521753
        nfp    =  4
        B2c    =  0.09655841297334804
        p2     =  0.0
        iota   =  -1.2865657607802559
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  0.15
        coilSeparation = 0.15
        targetValue = 0.08
        nCoilsPerNFP = 4
        # DMerc mean  =  0.0
        # Max elongation  =  2.503410982828086
        # B20 variation = 0.13561097312595394
        # Max |X20| = 0.6195365938105948
        # Max |X3c1| = 0.337460934254934
        # gradgradB inverse length:  2.5361897898745758
        # objective function:  829.9393490682932
        # Quasi-helically symmetric solution with N = -4.0
    if ind==8:
        name   = 'QH_NFP4_vac_2'
        rc     = [ 1.0,0.24166101615219052,0.03507911723897173,0.0030714758562874213,-0.00021710659761328953,-0.00013811263507303878,-2.2973633833468426e-05,-1.1109924853845887e-06 ]
        zs     = [ 0.0,0.21477390023482557,0.03218163581234407,0.002938884475059996,-0.0001826969770848549,-0.0001307464910730861,-2.2424676912296043e-05,-7.394997689498148e-07 ]
        etabar =  1.973712439428752
        nfp    =  4
        B2c    =  0.09140051741116116
        p2     =  0.0
        iota   =  -1.28572842947075
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  0.15
        coilSeparation = 0.15
        targetValue = 0.08
        nCoilsPerNFP = 4
        # DMerc mean  =  0.0
        # Max elongation  =  2.5037679503738364
        # B20 variation = 0.1497674070149666
        # Max |X20| = 0.6236107132927807
        # Max |X3c1| = 0.3400976067483434
        # gradgradB inverse length:  2.4961342898583
        # objective function:  3946.4353590826954
        # Quasi-helically symmetric solution with N = -4.0
    if ind==9:
        name   = 'QH_NFP4_plasma_2'
        rc     = [ 1.0,0.5178784404326631,-0.04424782585887367,0.009423311753794298,0.0006855570108155253,0.002208553146897593,-0.0013409321817794506,0.00029218322820064133 ]
        zs     = [ 0.0,0.5160627570612976,-0.01970682168578155,0.020006251376125224,-0.0008786505301953536,0.002388275820849969,-0.0005682651952807946,0.00045725410482900576 ]
        etabar =  1.096487001424194
        nfp    =  4
        B2c    =  -1.8171228141489864
        p2     =  -500000.0
        iota   =  -2.903776561280769
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  0.05
        coilSeparation = 0.15
        targetValue = 0.08
        nCoilsPerNFP = 4
        # DMerc mean  =  0.07021989443710663
        # DWell mean  =  1.6178637456192302
        # Max elongation  =  3.523778403753809
        # B20 variation = 0.6876978773560349
        # Max |X20| = 4.589959585465102
        # Max |X3c1| = 0.9484335362286295
        # gradgradB inverse length:  5.042967658568867
        # objective function:  2950.951650176265
    if ind==10:
        name   = 'QH_NFP5_vac'
        rc     = [ 1.0,0.1523094867470661,0.013157767646675112,0.0005005599925381512,-8.210272043711552e-05,-1.9025105246839816e-05,-1.5831051579313817e-06,3.0694707972137405e-09 ]
        zs     = [ 0.0,0.14162848726009797,0.012463696312271524,0.00046997304207690044,-8.194324747551428e-05,-1.9702979966612966e-05,-2.103228060854488e-06,-2.755334578689467e-08 ]
        etabar =  2.447727452515191
        nfp    =  5
        B2c    =  -0.22565885565267826
        p2     =  0.0
        iota   =  -1.302123928937576
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/6
        coilSeparation = 0.15
        targetValue = 0.08
        nCoilsPerNFP = 4
        # DMerc mean  =  0.0
        # DWell mean  =  0.0
        # Max elongation  =  2.2818027473251723
        # B20 variation = 0.06509325814329947
        # Max |X20| = 0.6689709806269559
        # Max |X3c1| = 0.32553930452700025
        # gradgradB inverse length:  3.1075312043749643
        # objective function:  1160.5572185715787
    if ind==11:
        name   = 'QH_NFP5_plasma'
        rc     = [ 1.0,0.15222503745286767,0.008643014779434653,-0.001528645358505054,-0.0003459957233227212,1.1655721389926272e-05,1.1307568806314385e-05,9.118071743105142e-07 ]
        zs     = [ 0.0,0.12529932261137572,0.007678578044264318,-0.0013127428337132072,-0.0003164361498580318,9.631344765435701e-06,1.0769685715917072e-05,3.849811197964361e-07 ]
        etabar =  2.788143445397788
        nfp    =  5
        B2c    =  0.34632953958054613
        p2     =  -1000000.0
        iota   =  -1.467194623787214
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/12
        coilSeparation = 0.15
        targetValue = 0.08
        nCoilsPerNFP = 4
        # DMerc mean  =  0.011457934470162856
        # DWell mean  =  0.55319027048864
        # Max elongation  =  2.4990345338701117
        # B20 variation = 0.07299741678870753
        # Max |X20| = 3.3566189051288644
        # Max |X3c1| = 0.9297476417013978
        # gradgradB inverse length:  5.076142325951124
        # objective function:  3353.7058762894367
    return stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP
