&physicslist
 Igeometry   =         2
 Istellsym   =         1
 Lfreebound  =         0
 phiedge     =   2.000000000000000E+00
 curtor      =   0.000000000000000E+00
 curpol      =   0.000000000000000E+00
 gamma       =   0.000000000000000E+00
 Nfp         =         4
 Nvol        =         1
 Mpol        =         10
 Ntor        =         12
 Lrad        =         6
 tflux       =   2.764900000000000E+00 
 pflux       =   1.475300000000000E-01                  
 helicity    =  -1.000000000000000E-01 
 pscale      =   1.000000000000000E+00
 Ladiabatic  =         0
 pressure    =   2.500000000000000E-02 
 adiabatic   =   0.000000000000000E+00 
 mu          =         1.0 
 Lconstraint =         0
 pl          =                       0                      0                      0
 ql          =                       0                      0                      0
 pr          =                       0                      0                      0
 qr          =                       0                      0                      0
 iota        =                       0                      0                      0
 lp          =                       0                      0                      0
 lq          =                       0                      0                      0
 rp          =                       0                      0                      0
 rq          =                       0                      0                      0
 oita        =                       0                      0                      0
 mupftol     =   1.000000000000000E-16
 mupfits     =         8

 ! Axis shape
 Ras         =   0.000000000000000E+00
 Zac         =   0.000000000000000E+00

!----- Boundary Parameters -----
Rwc(0,0)    =  0.000000000000000E+00 Zws(0,0)    =  0.000000000000000E+00 Rws(0,0)    =  0.000000000000000E+00 Zwc(0,0)    =  0.000000000000000E+00
Vns(0,0)    =  0.000000000000000E+00 Bns(0,0)    =  0.000000000000000E+00 Vnc(0,0)    =  0.000000000000000E+00 Bnc(0,0)    =  0.000000000000000E+00
/
&numericlist
 Linitialize =         1
 Lzerovac    =         0
 Ndiscrete   =         2
 Nquad       =        -1
 iMpol       =        -4
 iNtor       =        -4
 Lsparse     =         0
 Lsvdiota    =         0
 imethod     =         3
 iorder      =         2
 iprecon     =         1
 iotatol     =  -1.000000000000000E+00
 Lextrap     =         0
 Mregular    =        -1
/
&locallist
 LBeltrami   =         4
 Linitgues   =         1
 Lposdef     =         0
/
&globallist
 Lfindzero   =         2
 escale      =   0.000000000000000E+00
 opsilon     =   1.000000000000000E+00
 pcondense   =   2.000000000000000E+00
 epsilon     =   0.000000000000000E+00
 wpoloidal   =   1.000000000000000E+00
 upsilon     =   1.000000000000000E+00
 forcetol    =   1.000000000000000E-12
 c05xmax     =   1.000000000000000E-06
 c05xtol     =   1.000000000000000E-10
 c05factor   =   1.000000000000000E-02
 LreadGF     =         F
 mfreeits    =         0
 gBntol      =   1.000000000000000E-06
 gBnbld      =   6.660000000000000E-01
 vcasingeps  =   1.000000000000000E-12
 vcasingtol  =   1.000000000000000E-08
 vcasingits  =         8
 vcasingper  =         1
/
&diagnosticslist
 odetol      =   1.000000000000000E-07
 nPpts       =         500
 nPtrj       =         10    
 LHevalues   =         F
 LHevectors  =         F
 LHmatrix    =         F
 Lperturbed  =         0
 dpp         =        -1
 dqq         =        -1
 Lcheck      =         1
 Ltiming     =         F
/
&screenlist
Wjo00aa      = F
Wpp00aa      = F
/
