from qsc import Qsc
from simsopt import make_optimizable
import numpy as np

def get_stel(ind,nphi=251,r_edge=0.06,coilSeparation = 0.1,targetValue = 0.08,nCoilsPerNFP = 6):
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
        delta   = 0.2
        stel   =  make_optimizable(Qsc(rc=rc,zs=zs, nfp=1, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta))
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
        rc     = [ 1.0,0.15192176276052208,0.008513568102712886,-0.0015165238838419468,-0.0003309423982280693,1.4058643082522878e-05,1.0787776194275778e-05,7.234161726187371e-07 ]
        zs     = [ 0.0,0.12539540811386873,0.00763044187826761,-0.0012985697094884045,-0.0003048185734085636,1.1728461451924552e-05,1.0356948718321514e-05,1.63620547094035e-07 ]
        etabar =  2.799836576915385
        nfp    =  5
        B2c    =  0.3406199032028608
        p2     =  -1000000.0
        iota   =  -1.4767262784829156
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/14
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 4
        # DMerc mean  =  0.06856092856963991
        # DWell mean  =  0.6186922013208368
        # Max elongation  =  2.513512334926185
        # B20 variation = 0.08175755660651873
        # Max |X20| = 3.4921400643679172
        # Max |X3c1| = 0.8888517092742774
        # gradgradB inverse length:  5.132758420843008
        # objective function:  3487.706869820473
    if ind==12:
        name   = 'QA_NFP2_vac_3'
        rc     = [ 1.0,-0.17012557797445627,0.02424514888601057,-0.0038379761742939208,0.0006167734737180552,-9.241415312428205e-05,1.1223574404549308e-05,-7.895289752143264e-07 ]
        zs     = [ 0.0,0.1856522434048165,-0.02490314192591721,0.0038872855790383835,-0.000607802031532201,8.622822576822008e-05,-8.823962948687036e-06,1.9218858876105486e-07 ]
        etabar =  0.8402138648190306
        nfp    =  2
        B2c    =  0.40056018426804035
        p2     =  0.0
        iota   =  0.4075404037099618
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/6
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # DMerc mean  =  0.0
        # DWell mean  =  0.0
        # Max elongation  =  3.0004650857301463
        # B20 variation = 0.31232073723082865
        # Max |X20| = 0.7715497266667112
        # Max |X3c1| = 0.32307112724204834
        # gradgradB inverse length:  1.8789608084254732
        # objective function:  1915.0558249149124
    if ind==13:
        name   = 'QA_NFP2_plasma_3'
        rc     = [ 1.0,-0.1581077057789473,0.009353750115945785,-0.0004333034563751681,0.00014991469720175874,-4.608369248451757e-05,1.5891197499457243e-05,-5.687144986446836e-06,1.4426861835371922e-06 ]
        zs     = [ 0.0,0.13974198471550384,-0.00936051798763442,-0.0006362118365517266,0.00023793862122763505,-4.3004251806679115e-05,1.4379505424900648e-05,-6.443506237125145e-06,1.3338555611606268e-06 ]
        etabar =  0.5638571249135147
        nfp    =  2
        B2c    =  -1.998280780569248
        p2     =  -500000.0
        iota   =  0.40267356191820125
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/200
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # DMerc mean  =  0.1100560754666067
        # DWell mean  =  0.23601756461607715
        # Max elongation  =  5.799767819549366
        # B20 variation = 0.575425433011612
        # Max |X20| = 3.783584114559825
        # Max |X3c1| = 1.709466004997212
        # gradgradB inverse length:  2.6302743285009753
        # objective function:  18910.68644864485
    if ind==14:
        name   = 'QH_NFP2_vac'
        rc     = [ 1.0,-0.4920919671769627,-0.1515285587774181,-0.03291373781171326,-0.00869589281363886,-0.0028873798690479944,-0.001211205091384429,-0.0004183690660154582,-0.00014192564344839165,-1.7056197840101293e-05 ]
        zs     = [ 0.0,0.4529808323399858,0.1540716760408071,0.03506279317499001,0.009381964494053907,0.00321833101837483,0.0012654212518138836,0.0003911925161434127,0.00011280462111120577,5.542913205568535e-06 ]
        etabar =  0.9166239780365532
        nfp    =  2
        B2c    =  0.2923484029933731
        p2     =  0.0
        iota   =  1.9132977908943651
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/100
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # DMerc mean  =  0.0
        # DWell mean  =  0.0
        # Max elongation  =  6.168100986462805
        # B20 variation = 1.3875369286719117
        # Max |X20| = 3.5483983426192385
        # Max |X3c1| = 0.5005519164293439
        # gradgradB inverse length:  3.663906263015699
        # objective function:  9025.891130916021
    if ind==15:
        name   = 'QH_NFP2_plasma'
        rc     = [ 1.0,-0.728522911259118,0.023009266116872287,0.021266172975764885,0.002809711754763311,0.0031348954054080565,-0.00022602852110138063,-0.0009201669022045536,8.020840841274659e-05 ]
        zs     = [ 0.0,0.6043120487431142,-0.08528928868975037,0.012932584671270031,-0.018044335147057437,0.0029504773541039584,-0.0008483912934742645,0.0010017207023001915,-5.689524784332993e-05 ]
        etabar =  0.5979274719360609
        nfp    =  2
        B2c    =  -1.56289502358437
        p2     =  -500000.0
        iota   =  1.3827723353676435
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/100
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # DMerc mean  =  0.44644186064692326
        # DWell mean  =  0.7567999119923016
        # Max elongation  =  6.067616547715366
        # B20 variation = 1.9673140621869392
        # Max |X20| = 10.070618262934925
        # Max |X3c1| = 1.5093242120957913
        # gradgradB inverse length:  5.805802707802815
        # objective function:  7536.864772783507
    if ind==16:
        name   = 'QA_NFP3_vac'
        rc     = [ 1.0,0.050224142623339005,0.0024648591290499142,0.0001119266286864647,4.182627242567558e-06,2.0737623981586257e-07 ]
        zs     = [ 0.0,-0.051219095337737844,-0.0024220896528694405,-0.00011375374153045005,-4.344635778503563e-06,-3.0080665491370894e-07 ]
        etabar =  0.9256884588145107
        nfp    =  3
        B2c    =  2.1121434634154563
        p2     =  0.0
        iota   =  0.4076794866225493
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  0.1
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # DMerc mean  =  0.0
        # DWell mean  =  0.0
        # Max elongation  =  2.3201965391823505
        # B20 variation = 0.1926020792957739
        # Max |X20| = 5.270117027956928
        # Max |X3c1| = 1.016402947243469
        # gradgradB inverse length:  2.50749246531078
        # objective function:  1646.6747729026563
    if ind==17:
        name   = 'QA_NFP3_plasma'
        rc     = [ 1.0,0.05089081600407128,0.0017011902646358252,0.0001105470109461782,1.8120572839703084e-05,-2.069567917922587e-06 ]
        zs     = [ 0.0,-0.0596243608216842,-0.0031234854259240916,-5.769159501725202e-05,1.2708716854131408e-05,1.6106901459541407e-06 ]
        etabar =  0.5891397568351997
        nfp    =  3
        B2c    =  -3.212844879693586
        p2     =  -500000.0
        iota   =  0.4074195627908438
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/12
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # DMerc mean  =  0.28079321053974654
        # DWell mean  =  0.38968505483903443
        # Max elongation  =  5.484476223101584
        # B20 variation = 0.7147628259667416
        # Max |X20| = 5.164336109862278
        # Max |X3c1| = 1.7507156446389611
        # gradgradB inverse length:  3.4190791073094893
        # objective function:  5651.888129153469
    if ind==18:
        name   = 'QH_NFP3_vac'
        rc     = [ 1.0,0.6017188333339238,-0.05013160291495114,-0.0059509156974514804,0.0046395002057484855,0.0005181849976044792,-0.0009021840397330221,3.2091159445160866e-05,0.00027903445923394487,-9.405544457236126e-05 ]
        zs     = [ 0.0,-0.5670143663357554,0.04828670600798458,0.0012712834089255308,-0.0019043062624112754,0.00048209215902020384,5.270817704207781e-05,5.599588943640311e-05,-0.0001442108409911305,0.0001292277548494176 ]
        etabar =  1.1468556706440256
        nfp    =  3
        B2c    =  -0.25130771525326506
        p2     =  0.0
        iota   =  2.2850300916357633
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/50
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # DMerc mean  =  0.0
        # DWell mean  =  0.0
        # Max elongation  =  3.7123657848388616
        # B20 variation = 0.9305790803911311
        # Max |X20| = 2.6911230299614943
        # Max |X3c1| = 1.1917258552291847
        # gradgradB inverse length:  3.5160403007603254
        # objective function:  2825.8995800383113
    if ind==19:
        name   = 'QH_NFP4_vac_3'
        rc     = [ 1.0,0.23534671348568964,0.032370741935558243,0.0028930588676929696,-9.979729683116764e-06,-4.947079146209497e-05,-6.7276829969791605e-06,-2.1686735790258402e-07,-2.1995177184739665e-08 ]
        zs     = [ 1.0,0.21350397651117003,0.030823452391198828,0.0029006895978590185,1.2731238475299125e-05,-4.385318742563811e-05,-5.201229091310175e-06,1.9385975965326248e-07,3.3406806412235074e-08 ]
        etabar =  2.009861401881025
        nfp    =  4
        B2c    =  0.014978229496948308
        p2     =  0.0
        iota   =  -1.2850584906730413
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/6
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # DMerc mean  =  0.0
        # DWell mean  =  0.0
        # Max elongation  =  2.5840289881380625
        # B20 variation = 0.09588435631261749
        # Max |X20| = 0.6648220678804774
        # Max |X3c1| = 0.3579745843831077
        # gradgradB inverse length:  2.4650583266372568
        # objective function:  543.7323144010619
    if ind==20:
        name   = 'QH_NFP4_plasma_3'
        rc     = [ 1.0,0.48901681580866646,0.07866833907333268,-0.003988981256181756,-0.0038518075650708157,0.0009832042786224607,0.0014020574385464684,0.0005469202893270081,8.469474285514098e-05 ]
        zs     = [ 1.0,-0.40443216690777006,-0.06860068919954514,-0.0017381573489314553,-4.663777216347447e-05,-0.001882013016647241,-0.0014023624219547796,-0.0005515360428790581,-0.00010482758630573898 ]
        etabar =  1.1774213377118163
        nfp    =  4
        B2c    =  -0.9336381190413221
        p2     =  -500000.0
        iota   =  2.44870850972142
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2))
        r_edge =  1/22
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # DMerc mean  =  0.47429340726938696
        # DWell mean  =  1.1135905853158308
        # Max elongation  =  2.791909049149931
        # B20 variation = 0.8168770201186457
        # Max |X20| = 3.4000759759030963
        # Max |X3c1| = 2.7302771176927543
        # gradgradB inverse length:  4.696956784587547
        # objective function:  2394.4241486287015
    if ind==21:
        name   = 'QA_NFP2_r1'
        rc     = [ 1.0,-0.18196013413025666,0.020999614567506025,-0.0021866610299392767,0.0001738008300926121,-3.4546322306952025e-06,-1.9298141483514003e-06 ]
        zs     = [ 0.0,0.15577654319784595,-0.01839877656873839,0.001868901595355642,-0.0001391416351411979,-1.360580336280459e-07,1.706575772828251e-06 ]
        etabar =  0.8008233540885676
        nfp    =  2
        B2c    =  0.0
        p2     =  0.0
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, order='r1'))
        iota   =  0.4198790908227896
        r_edge =  0.1
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # Max |X1c| = 1.6191472702006777
        # Max elongation  =  3.3445701395505636
        # objective function:  23.501823574911157
    if ind==22:
        name   = 'QH_NFP4_r1'
        rc     = [ 1.0,0.5244368868248144,0.1699930341221515,0.04549135595299858,0.0091392319599508,0.0011702750200607542,5.622120041171967e-05 ]
        zs     = [ 0.0,0.49268004097859763,0.16924278110761584,0.047119948399724026,0.00970139003772989,0.0011972374113402445,2.8842577610049084e-05 ]
        etabar =  1.5049682592480425
        nfp    =  4
        B2c    =  0.0
        p2     =  0.0
        stel   =  make_optimizable(Qsc(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, order='r1'))
        iota   =  -2.2663872307770703
        r_edge =  0.06
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # Max |X1c| = 1.4298303819278673
        # Max elongation  =  2.881315728145213
        # objective function:  52.66809016147327
    if ind==23:
        name   = 'QI_NFP1_r1'
        rc     = [ 1.0,0.0,-0.2,0.0,0.0,0.0,0.0 ]
        zs     = [ 0.0,0.0,0.06150199148349551,-0.024525761300959653,0.0035708217656314364,0.0036905543378600875,0.0009216112298971144 ]
        B0_vals = [ 1.0,0.1 ]
        d_svals = [ 0.0,1.4248780525666942,-0.02873217816997125,0.01960207141637209,0.003156746990570063,0.005069833884773372,0.0 ]
        delta   = 0.6283185307179586
        stel   =  make_optimizable(Qsc(rc=rc,zs=zs, nfp=1, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi,phi_shift=1/3, omn=True, delta=delta))
        iota   =  0.4198331396220715
        r_edge =  0.03
        coilSeparation = 0.1
        targetValue = 0.08
        nCoilsPerNFP = 6
        # Max |curvature| = 1.3967408412720728
        # Max |X1c| = 0.8427827161080429
        # Max elongation  =  4.2440217493138075
        # objective function:  240.9861459479161
    if ind==24:
        name   = 'QI_NFP2_r1'
        rc     = [ 1.0,0.0,-0.058823529411764705,0.0,0.0,0.0,0.0,0.0,0.0 ]
        zs     = [ 0.0,0.0,0.01209128995076437,0.0,0.003042581799538996,0.0,7.945611170593634e-06,0.0,-1.9331543658597153e-05 ]
        B0_vals = [ 1.0,0.17696326718965785 ]
        d_svals = [ 0.0,0.8725104433010822,-0.007125723995617353,-0.007859854466243633,0.025890615049836384 ]
        delta   = 0.49344993841566404
        nfp     = 2
        stel    =  make_optimizable(Qsc(rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta))
        iota   =  3.0792207535685585
        r_edge = 0.001
        # mean gradB inverse length: 7.886964648107345
        # Max elongation = 6.351384904251
        # objective function:  662.2277709720398
    try:
        stel.omn
    except:
        stel.omn = False
    return stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP
