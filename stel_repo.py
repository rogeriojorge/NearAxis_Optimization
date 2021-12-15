from QSCwrapper import QSCWrapper
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2)
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
        stel   =  QSCWrapper(rc=rc,zs=zs, nfp=1, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r2', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, B2c=B2c, order='r3', p2=p2)
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, order='r1')
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
        stel   =  QSCWrapper(rc=rc, zs=zs, etabar=etabar, nfp=nfp, nphi=nphi, order='r1')
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
        rc      = [ 1.0,0.0,-0.2,0.0,0.0,0.0,0.0,0.0,0.0 ]
        zs      = [ 0.0,0.0,0.06299249203124066,0.0,-0.003505048607115181,0.0,0.002332062917754121,0.0,0.00017792827916776854 ]
        B0_vals = [ 1.0,0.1362166563376213 ]
        d_svals = [ 0.0,1.0406687319250656,-0.008459506894986953,0.11710250798554353,0.0607991113831632 ]
        delta   = 0.4854813680461505
        nfp     = 1
        stel    =  QSCWrapper(rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta)
        iota    =  0.5317023102263984
        # mean gradB inverse length: 1.770579342267488
        # Max elongation = 5.367464071620833
        # objective function:  365.69346679821996
    if ind==24:
        name   = 'QI_NFP2_r1'
        rc      = [ 1.0,0.0,-0.058823529411764705,0.0,0.0,0.0,0.0,0.0,0.0 ]
        zs      = [ 0.0,0.0,-0.018643361988803665,0.0,-0.0013564016185653712,0.0,-0.00046521641214137557,0.0,6.561558212708543e-05 ]
        B0_vals = [ 1.0,0.12150771699343063 ]
        omn_method ='buffer'
        k_buffer = 3
        d_svals = [ 0.0,1.2066424359614483,0.0839996920762284,-0.0481193237846653,0.009748360571682517 ]
        delta   = 0.8577502034294207
        nfp     = 2
        stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta)
        iota    =  -2.7967628055761162
        # mean gradB inverse length: 2.165558583516918
        # Max elongation = 4.81497617522344
        # objective function:  2587.4620377027527
    if ind==25:
        name   = 'QI_NFP1_r2'
        rc      = [ 1.0,0.0,-0.3128567868276656]#,0.0,0.033421754548621244,0.0,-0.00010502414022251839,0.0,0.0,0.0,0.0 ]
        zs      = [ 0.0,0.0,-0.27212625571920984]#,0.0,0.044167669563470506,0.0,-6.574629652081523e-05,0.0,-9.893447860006291e-06,0.0,9.277608267145067e-06 ]
        B0_vals = [ 1.0,0.15810011629203774 ]
        omn_method ='non-zone'
        k_buffer = 3
        k_second_order_SS   = 0.1
        d_over_curvature   = 0.1
        d_svals = [ 0.0,0.1]
        delta   = 0.8
        nfp     = 1
        B2s_svals = [ 0.0,0.1]#,0.03283462033677047,0.004837653606915959,0.03466436789051,0.012588451775448475,0.0034798016087720404,0.012368193922088775,0.03844118710392006,0.002403469150739275,0.01163670504914525 ]
        B2c_cvals = [ 0.1,0.1]#,0.06344485624821924,0.026187767224228442,0.05839086377958894,0.007915937165665798,0.036743548309213274,0.036842540401104525,0.0142687978041074,-0.013198287181027281,-0.0072079191108149875,0.0 ]
        B2s_cvals = [ -0.0,0.0]#,-0.022473047295254987,-0.08265049514305536,0.41342725761575694,-0.2931657967710277,0.31542024309141736,-0.21005748857764805,0.6787636103131851,-0.1240125056565712,0.1242848766908088,0.0 ]
        B2c_svals = [ 0.0,0.0]#,-0.57516292896785,0.8948531647751348,-3.0875825530983967,-0.06329756906736621,-1.7890702692964318,1.039759175233141,-0.1690212827616995,-0.16762613422413772,0.015322897225364228,0.0,0.0 ]
        p2      =  0.0
        stel    =  QSCWrapper(d_over_curvature=d_over_curvature, d_svals = d_svals, omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, nphi=nphi, omn=True, delta=delta, B2c_cvals=B2c_cvals, B2s_svals=B2s_svals, p2=p2, order='r3', k_second_order_SS=k_second_order_SS, B2s_cvals=B2s_cvals, B2c_svals=B2c_svals)
        iota    =  0.5327817951931755
        # DMerc mean  = 0.0
        # DWell mean  = 0.0
        # DGeod mean  = 0.0
        # B20 mean = 0.02294560802894031
        # B20QI_deviation_max = 5.453883074067173
        # B2cQI_deviation_max = 3.0062562564563393
        # B2sQI_deviation_max = 0.19536838040285726
        # Max |X20| = 7.198294663064016
        # Max |X3c1| = 1.5390774968163563
        # gradgradB inverse length: 4.2411938841600865
        # d2_volume_d_psi2 mean =, -15.495596160684793
        # mean gradB inverse length: 1.6170809917518607
        # Max elongation = 9.601832283747992
        # objective function:  18559.816948294818
    if ind==26:
        name   = 'QI_NFP1_r1_k_second_order_SS'
        rc      = [ 1.0,0.0,-0.3215933546627156,0.0,0.04374904855292157,0.0,-0.0036693797861105035,0.0,0.0 ]
        zs      = [ 0.0,0.0,0.1761894885594696,0.0,-0.03212484216647393,0.0,0.004288782017242512,0.0,-0.0002197285688705138 ]
        B0_vals = [ 1.0,0.1494580658305839 ]
        omn_method ='non-zone'
        k_buffer = 4
        k_second_order_SS   = -0.2945414999865881
        d_over_curvature   = 0.0
        d_svals = [ 0.0,-1.2734735798960166,-0.011728655941655858,0.06995110936196922,0.02501702564275091 ]
        delta   = 0.1
        nfp     = 1
        stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        iota    =  0.5171335081684099
        # max curvature'(0): 1.3882301507705863
        # max d'(0): -0.9487763169897289
        # max gradB inverse length: 2.405695177762281
        # Max elongation = 3.8972846007080695
        # objective function:  29.070861104483427
    if ind==27:
        name   = 'QI_NFP1_r1_d_over_curvature'
        rc      = [ 1.0,0.0,-0.37510983559429445,0.0,0.0637209503308546,0.0,-0.005613702098731244,0.0,0.0 ]
        zs      = [ 0.0,0.0,0.26383174888812866,0.0,-0.06099954215159932,0.0,0.0061269955241693205,0.0,0.00010252920801029004 ]
        B0_vals = [ 1.0,0.162749191944802 ]
        omn_method ='non-zone'
        k_buffer = 3
        k_second_order_SS   = 0.0
        d_over_curvature   = 0.6481405174969549
        d_svals = [ 0.0,0.05367173687497751,0.05302220930256652,-0.03477317136890787,-0.02028376459985586 ]
        delta   = 0.1
        nfp     = 1
        stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        iota    =  0.5998169708699606
        # max curvature'(0): 0.4844311966624933
        # max d'(0): 0.2882883585584308
        # max gradB inverse length: 1.5653267149603343
        # Max elongation = 4.161353494303795
        # objective function:  26.78191718287084
    if ind==28:
        name   = 'QI_NFP1_r1_best'
        rc      = [ 1.0,0.0,-0.31156698348096523,0.0,0.03281381867087213,0.0,0.0,0.0,0.0,0.0,0.0 ]
        zs      = [ 0.0,0.0,0.32596628468713607,0.0,-0.05148075653531195,0.0,0.00035721010864138127,0.0,-0.00025120004844148816,0.0,-9.681058605149824e-06 ]
        B0_vals = [ 1.0,0.1611858881581418 ]
        omn_method ='non-zone'
        k_buffer = 3
        d_over_curvature   = 0.4882661440608365
        d_svals = [ 0.0,0.005743059145943575,0.0010017023515223697,-0.0030768469053832734,-0.0011614585411610532,-5.826410883440976e-06 ]
        delta   = 0.1
        nfp     = 1
        stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature)
        iota    =  0.7056219688553746
        # max curvature'(0): 0.2068350328820353
        # max d'(0): 0.09483603670101105
        # mean gradB inverse length: 1.2397283031308683
        # Max elongation = 5.402736072814011
        # objective function:  25243.16894816218
    if ind==29:
        name   = 'QI_NFP1_r1_test'
        rc      = [ 1.0,0.0,-0.2]
        zs      = [ 0.0,0.0,-0.3]
        B0_vals = [ 1.0,0.16]
        omn_method ='non-zone'
        k_buffer = 3
        k_second_order_SS   = 0.0
        d_over_curvature   = 0.5
        d_svals = [ 0.0,0.01]
        delta   = 0.1
        nfp     = 1
        stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature)
        # (stel.get_elongation, 0.0, 5e-1/stel.nphi),
        # (stel.get_d, 0.0, 2e+1/stel.nphi),
        # (stel.get_d_svals, 0.0, 1e2),
        # (stel.get_min_R0_penalty, 0.0, 2e1),
        # (stel.get_min_Z0_penalty, 0.0, 2e1),
        # (stel.get_B0_well_depth,0.15, 3e2),
        # (stel.get_inv_L_grad_B, 0.0, 3e-2),
        # (stel.get_d_d_d_varphi_at_0,0.0,2e0),
        # (stel.get_alpha_deviation,0.0,6e+1/stel.nphi)
        # nIterations = 150
        # abs_step_array = [1e-1,1e-2,1e-3,1e-4,1e-6]
        # rel_step_array = [1e-1,1e-2]
        # max_fourier_coefficients = 5
        # r_edge = 1/8
        # ftol = 1e-3
        # nCoilsPerNFP   = 20
        # targetValue    = 0.042
        # coilSeparation = 0.375
        # stel.min_R0_threshold = 0.4
        # rc      = [ 1.0,0.0,-0.3128567868276656,0.0,0.033421754548621244,0.0,-0.00010502414022251839,0.0,0.0,0.0,0.0 ]
        # zs      = [ 0.0,0.0,-0.27212625571920984,0.0,0.044167669563470506,0.0,-6.574629652081523e-05,0.0,-9.893447860006291e-06,0.0,9.277608267145067e-06 ]
        # B0_vals = [ 1.0,0.15810011629203774 ]
        # omn_method ='non-zone'
        # k_buffer = 3
        # k_second_order_SS   = 0.0
        # d_over_curvature   = 0.5356625931840477
        # d_svals = [ 0.0,-0.00216107652370974,0.0004036139788440245,1.2115865686228548e-06,-1.5895958585498093e-06,-1.268277691787932e-05 ]
        # delta   = 0.1
        # nfp     = 1
        # stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        # iota    =  -0.6545438578820277
        # # max curvature'(0): 0.16346295965506322
        # # max d'(0): 0.08614102277953053
        # # max gradB inverse length: 1.4851665891326913
        # # Max elongation = 5.02330125051389
        # # objective function:  35.67338668936966
        #
        #
        # rc      = [ 1.0,0.0,-0.40324108353845844,0.0,0.07641000989682407,0.0,-0.007642290555505867,0.0,0.0,0.0,0.0,0.0,0.0 ]
        # zs      = [ 0.0,0.0,-0.27043849940003517,0.0,0.07222260076564689,0.0,-0.007717369810843748,0.0,-0.00016773293342324324,0.0,-0.0002892845740449629,0.0,4.807952480165844e-05 ]
        # B0_vals = [ 1.0,0.16658668391385073 ]
        # omn_method ='non-zone'
        # k_buffer = 3
        # k_second_order_SS   = 0.0
        # d_over_curvature   = 0.5271692188343875
        # d_svals = [ 0.0,0.008160054531665732,0.0021653150897661755,-0.002708012747006448,-0.0010492663711443879,4.962266658327616e-05,-0.0002642297314102418 ]
        # delta   = 0.1
        # nfp     = 1
        # stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        # iota    =  -0.6797193989743312
        # # max curvature'(0): 0.25151830693877797
        # # max d'(0): 0.1314278745988563
        # # max gradB inverse length: 1.5683895752828814
        # # Max elongation = 4.928278600925702
        # # objective function:  34.652798293914636
        #
        #
        #
        #
        # (stel.get_elongation, 0.0, 4e-1/stel.nphi),
        # (stel.get_d, 0.0, 2e+1/stel.nphi),
        # (stel.get_d_svals, 0.0, 1e2),
        # (stel.get_min_R0_penalty, 0.0, 3e1),
        # (stel.get_min_Z0_penalty, 0.0, 3e1),
        # (stel.get_B0_well_depth,0.16, 2e2),
        # (stel.get_inv_L_grad_B, 0.0, 3e-2),
        # (stel.get_d_d_d_varphi_at_0,0.0,2e0),
        # (stel.get_alpha_deviation,0.0,6e+1/stel.nphi)
        # rc      = [ 1.0,0.0,-0.4056622889934463,0.0,0.07747378220100756,0.0,-0.007803860877024245,0.0,0.0,0.0,0.0,0.0,0.0 ]
        # zs      = [ 0.0,0.0,-0.24769666390049602,0.0,0.06767352436978152,0.0,-0.006980621303449165,0.0,-0.0006816270917189934,0.0,-1.4512784317099981e-05,0.0,-2.839050532138523e-06 ]
        # B0_vals = [ 1.0,0.16915531046156507 ]
        # omn_method ='non-zone'
        # k_buffer = 3
        # k_second_order_SS   = 0.0
        # d_over_curvature   = 0.5183783762725197
        # d_svals = [ 0.0,0.003563114185517955,0.0002015921485566435,-0.0012178616509882368,-0.00011629450296628697,-8.255825435616736e-07,3.2011540526397e-06 ]
        # delta   = 0.1
        # nfp     = 1
        # stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        # iota    =  -0.6802282588093194
        # # max curvature'(0): 0.18193259676864582
        # # max d'(0): 0.0941731237760574
        # # max gradB inverse length: 1.6322554885663985
        # # Max elongation = 5.1003624635531954
        # # objective function:  32.61831563657005
        #
        #
        # (stel.get_elongation, 0.0, 5e-1/stel.nphi),
        # # (stel.get_sigma, 0.0, 1e+1/stel.nphi),
        # # (stel.get_torsion, 0.0, 1e+1/stel.nphi),
        # # (stel.get_curvature, 1.0, 5e-1/stel.nphi),
        # (stel.get_d, 0.0, 2e+1/stel.nphi),
        # (stel.get_d_svals, 0.0, 1e2),
        # # (stel, 'd_X1c_d_varphi', 0.0, 2e-2),
        # # (stel, 'd_Y1c_d_varphi', 0.0, 2e-2),
        # # (stel, 'd_Y1s_d_varphi', 0.0, 2-2),
        # (stel.get_min_R0_penalty, 0.0, 3e1),
        # (stel.get_min_Z0_penalty, 0.0, 3e1),
        # # (stel.get_delta, 0.0, 1e1),
        # (stel.get_B0_well_depth,0.16, 3e2),
        # (stel.get_inv_L_grad_B, 0.0, 3e-2),
        # (stel.get_d_d_d_varphi_at_0,0.0,1e0),
        # # (stel.get_d_curvature_d_varphi_at_0,0.0,5e-1),
        # # (stel.get_d_over_curvature,1.0,1e1),
        # # (stel.get_k_second_order_SS,0.0,5e0),
        # (stel.get_alpha_deviation,0.0,6e+1/stel.nphi)
        # rc      = [ 1.0,0.0,-0.4083538847571259,0.0,0.07891332771166197,0.0,-0.008101544521962805,0.0,0.0,0.0,0.0,0.0,0.0 ]
        # zs      = [ 0.0,0.0,-0.2464387627592985,0.0,0.06740672382378578,0.0,-0.006936056024592203,0.0,-0.0006797374565423773,0.0,-8.44760125464751e-05,0.0,6.87969853773234e-06 ]
        # B0_vals = [ 1.0,0.16721710271418216 ]
        # omn_method ='non-zone'
        # k_buffer = 3
        # k_second_order_SS   = 0.0
        # d_over_curvature   = 0.543005297653247
        # d_svals = [ 0.0,0.0010704905600812038,0.003444928960092528,-0.00040665107254926445,-0.0018086255197985862,1.925958756619766e-05,-1.5495741126771823e-05 ]
        # delta   = 0.1
        # nfp     = 1
        # stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        # iota    =  -0.663708939244891
        # # max curvature'(0): 0.3067294946728777
        # # max d'(0): 0.16606661468118356
        # # max gradB inverse length: 1.6447422928562212
        # # Max elongation = 4.928383967140109
        # # objective function:  34.49056743092893
        #
        #
        #
        #
        # rc      = [ 1.0,0.0,-0.2]
        # zs      = [ 0.0,0.0,-0.3]
        # B0_vals = [ 1.0,0.16]
        # omn_method ='non-zone'
        # k_buffer = 3
        # k_second_order_SS   = 0.0
        # d_over_curvature   = 0.5
        # d_svals = [ 0.0,0.01]
        # delta   = 0.1
        # nfp     = 1
        # (stel.get_elongation, 0.0, 4e-1/stel.nphi),
        # (stel.get_d, 0.0, 2e+1/stel.nphi),
        # (stel.get_d_svals, 0.0, 1e2),
        # (stel.get_min_R0_penalty, 0.0, 3e1),
        # (stel.get_min_Z0_penalty, 0.0, 3e1),
        # (stel.get_B0_well_depth,0.16, 2e2),
        # (stel.get_inv_L_grad_B, 0.0, 3e-2),
        # (stel.get_d_d_d_varphi_at_0,0.0,2e0),
        # (stel.get_alpha_deviation,0.0,6e+1/stel.nphi)
        rc      = [ 1.0,0.0,-0.4056622889934463,0.0,0.07747378220100756,0.0,-0.007803860877024245,0.0,0.0,0.0,0.0,0.0,0.0 ]
        zs      = [ 0.0,0.0,-0.24769666390049602,0.0,0.06767352436978152,0.0,-0.006980621303449165,0.0,-0.0006816270917189934,0.0,-1.4512784317099981e-05,0.0,-2.839050532138523e-06 ]
        B0_vals = [ 1.0,0.16915531046156507 ]
        omn_method ='non-zone-fourier'
        k_buffer = 3
        k_second_order_SS   = 0.0
        d_over_curvature   = 0.5183783762725197
        d_svals = [ 0.0,0.003563114185517955,0.0002015921485566435,-0.0012178616509882368,-0.00011629450296628697,-8.255825435616736e-07,3.2011540526397e-06 ]
        delta   = 0.1
        nfp     = 1
        stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        iota    =  -0.6802282588093194
        # max curvature'(0): 0.18193259676864582
        # max d'(0): 0.0941731237760574
        # max gradB inverse length: 1.6322554885663985
        # Max elongation = 5.1003624635531954
        # objective function:  32.61831563657005
    if ind==30:
        name   = 'QI_NFP2_r1'
        rc      = [ 1.0,0.0,-0.01]
        zs      = [ 0.0,0.0,0.01]
        B0_vals = [ 1.0,0.15]
        omn_method ='non-zone-fourier'
        k_buffer = 1
        k_second_order_SS   = 0.0
        d_over_curvature   = 0.5
        d_svals = [ 0.0,0.01]
        delta   = 0.1
        nfp     = 2
        stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature)
        #
        #
        rc      = [ 1.0,0.0,-0.08409539411347328,0.0,0.0075907297331848506,0.0,-0.000439832639503239,0.0,0.0 ]
        zs      = [ 0.0,0.0,0.054107403635000426,0.0,-0.006677539244689272,0.0,0.0005528215202865916,0.0,-3.6968992508928826e-05 ]
        rs      = [ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ]
        zc      = [ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ]
        sigma0  =  0.0
        B0_vals = [ 1.0,0.17873176948252556 ]
        omn_method ='non-zone-fourier'
        k_buffer = 1
        k_second_order_SS   = 0.0
        d_over_curvature   = 0.6758786723901443
        d_svals = [ 0.0,-0.29492716510830536,0.008093286550201534,0.043669675571235564,-0.0011914162939404638 ]
        delta   = 0.1
        nfp     = 2
        stel    =  QSCWrapper(sigma0 = sigma0, omn_method = omn_method, k_buffer=k_buffer, rs=rs,zc=zc, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        iota    =  1.2544901618034312
        # max curvature'(0): 1.2264129677330398
        # max d'(0): 0.5238114939385041
        # max gradB inverse length: 2.924661456312724
        # Max elongation = 5.361837970777982
        # objective function:  56.245844190585444
    if ind==31:
        name   = 'QI_NFP3_r1'
        rc      = [ 1.0,0.0,-0.01]
        zs      = [ 0.0,0.0,0.01]
        B0_vals = [ 1.0,0.15]
        omn_method ='non-zone-fourier'
        k_buffer = 1
        k_second_order_SS   = 0.0
        d_over_curvature   = 0.5
        d_svals = [ 0.0,0.01]
        delta   = 0.1
        nfp     = 3
        stel    =  QSCWrapper(omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature)
        #
        rc      = [ 1.0,0.0,-0.032774251422958894,0.0,0.0015303115216324316,0.0,-2.8454978422225952e-05,0.0,0.0 ]
        zs      = [ 0.0,0.0,0.024635017625994605,0.0,-0.001507030960075822,0.0,3.816171660582472e-05,0.0,-3.854618005163362e-06 ]
        rs      = [ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ]
        zc      = [ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ]
        sigma0  =  0.0
        B0_vals = [ 1.0,0.1968392611688621 ]
        omn_method ='non-zone-fourier'
        k_buffer = 1
        k_second_order_SS   = 0.0
        d_over_curvature   = 0.8066150911500598
        d_svals = [ 0.0,-0.5906721968730675,0.020087801467747076,0.06119519341933523,-0.00024309272075431754 ]
        delta   = 0.1
        nfp     = 3
        stel    =  QSCWrapper(sigma0 = sigma0, omn_method = omn_method, k_buffer=k_buffer, rs=rs,zc=zc, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        iota    =  1.7499216269893367
        # max curvature'(0): 2.778048368625772
        # max d'(0): 1.1369192206581655
        # max gradB inverse length: 3.422272769562569
        # Max elongation = 5.3348766008669495
        # objective function:  149.20785124911535
    if ind==32:
        name   = 'QI_NFP1_r1_smootherAlpha'
        rc      = [ 1.0,0.0,-0.2]
        zs      = [ 0.0,0.0,-0.3]
        B0_vals = [ 1.0,0.16]
        omn_method ='non-zone-smoother'
        k_buffer = 6
        p_buffer = 5
        k_second_order_SS   = 0.0
        d_over_curvature   = 0.5
        d_svals = [ 0.0,0.01]
        delta   = 0.1
        nfp     = 1
        stel    =  QSCWrapper(p_buffer=p_buffer, omn_method = omn_method, k_buffer=k_buffer, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature)
        #
        #
        rc      = [ 1.0,0.0,-0.4033958781759469,0.0,0.07635203004611708,0.0,-0.007594732970385295,0.0,0.0,0.0,0.0 ]
        zs      = [ 0.0,0.0,-0.24131310065314193,0.0,0.06382450374382645,0.0,-0.006016219817385928,0.0,-0.0006194662608072705,0.0,-6.333753439738566e-05 ]
        rs      = [ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ]
        zc      = [ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ]
        sigma0  =  0.0
        B0_vals = [ 1.0,0.17158156401368893 ]
        omn_method ='non-zone-smoother'
        k_buffer = 6
        p_buffer = 5
        k_second_order_SS   = 0.0
        d_over_curvature   = 0.5441152805054742
        d_svals = [ 0.0,0.0034051190420990563,0.003317625768761041,-0.0012623220508664378,-0.0017117338190760154,2.471035249896242e-05 ]
        delta   = 0.1
        nfp     = 1
        stel    =  QSCWrapper(sigma0 = sigma0, omn_method = omn_method, p_buffer = p_buffer, k_buffer=k_buffer, rs=rs,zc=zc, rc=rc,zs=zs, nfp=nfp, B0_vals=B0_vals, d_svals=d_svals, nphi=nphi, omn=True, delta=delta, d_over_curvature=d_over_curvature, k_second_order_SS=k_second_order_SS)
        iota    =  -0.6659650002271166
        # max curvature'(0): 0.16210885264588087
        # max d'(0): 0.08773816656521899
        # max gradB inverse length: 1.6540394959627556
        # Max elongation = 4.925160413321286
        # objective function:  30.98593629156751
    return stel, name, r_edge, coilSeparation, targetValue, nCoilsPerNFP