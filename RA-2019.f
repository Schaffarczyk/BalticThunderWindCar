C
C----------------------------------------------------------------------------------------
C
C     PROGRAMM opt  (C)  A.P. SCHAFFARCZYK UAS Kiel
C     march 1998, nov 2004, oct 2009, for CHINOOK, 23. august 2011
c     
C     basic equtions are from De Vries' AGARDOgraph 243 (1979), chapter 4.4 ff
C     and  
c     Mac Gaunaa et.al. proc. ewec 2009, Equation  (36) 
c
c     thorough revision on Sept 2019
c     note that de Vries heaviliy relies on "higher order" tip correction
c	
c---------------------------------------------------------------------------------------------
c
C----------------------------------------------------------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      CHARACTER*100 NAMEOUT1,NAMEOUT2, NAMEOUT3, NAMEOUT4 
      CHARACTER*100 NAMEIN
      CHARACTER*1   s1, s2 
      CHARACTER*6   steuer
c
      REAL PI, LAMBDA, la, B, GZ, etadr
C
C---- NUMERiC constants
C
      PI       =    4. * ATAN(1.)
      EPS      =    1.E-07
      EPSG     =    1.E-05
      MAXI     =  200
      MAXJ     =  100
      thetamax = pi/2.
C
c*************************************************************************
c********** start parameter *********************************************
c
c
	write(*,*)'#################### INPUT #####################'
c
C----  to change: re-compile
c
c     (1) car
c
c     vratio = design speed ratio 
c     etradr = drive-train efficiency 
c
      vratio   =   0.5
c
c     eta-Mech = 0.8
c     eta-Aero aus tip- und profil loss ca.  0.8
c
      etaM = 0.88
      etaae= 0.8
c
c     because  tip-loss and drag aerodynamic losses are included !
c
      etadr    =  etaM
c
      write(*,'(2(a10,f6.3))')'vra= ',vratio,'eta= ',etadr
c
c     (2) Rotor
c
C     number of B(lades), Lambda = Omega R_tip = TSR , G(leit)Z(ahl) = L2D ratio
C     Rtip = Tip-Radius
C
      B        =    3.0
      LAMBDA   =    5.
      Rtip     =    0.95
      write(*,'(3(a10,f6.3))')'B= ',B,'TSR= ',lambda,'Rtip= ',rtip
c
c     (3) Chosen profile Eppler 387, ReN = 200 k
c
c     design lift
c
      cLdes = 1.0
c
c     aoa-des(CL-des)
c
      aoa      = 6.0
c
c     L2D(at aoa-des)
c
      GZ       =    80
      write(*,'(3(a10,f6.3))')'Cl-des= ',cldes,'aoa-des= ',aoa,
     +                        'L2D= ',GZ
      write(*,*)
c
cc****************************************************************************
cc************** end parameter ***********************************************
cc
c
c !!!!!!!!!!!!!!!!! DON'T CHANGE from here !!!!!!!!!!!!!!!!!!
c     
      drag     = 1./gz
c
      intla    = int(10*lambda)
      i1 = intla/10
      i2 = intla-10*i1
      s1 = ACHAR(48+i1)
      s2 = achar(48+i2)
C
      IOOUT21   = 21
      OPEN(UNIT=IOOUT21,FILE='list',FORM='FORMATTED',STATUS='UNKNOWN')
c
      WRITE(21,'(5a12  )')'B','LAMBDA','GZ','Rtip','Design CL'
      WRITE(21,'(5f12.4)') B, LAMBDA, GZ, Rtip,CLdes
c
      Write(21,*)
C
      NAMEOUT1 = 'BTdes.' // s1 //s2
C
      IOOUT1   = 12
      OPEN(UNIT=IOOUT1,FILE=NAMEOUT1,FORM='FORMATTED',STATUS='UNKNOWN')
      WRITE(IOOUT1,'(7a12)')
     + 'x=r/Rtip','xx=x*la','c Cl/Rtip','r/mm','c/mm','Phi/Deg',' F'
c     
      WRITE(*,*)
      Write(*,*)"############### Blade Desiggn ################"
      WRITE(*,'(9a12)')
     + 'x=r/Rtip','xx=x*la','c Cl/Rtip','r/mm','c/mm','Phi/Deg',' F',
     + 'dcP loc','dcT'
C
      NAMEOUT2 = 'INDUCTION_FACTORS'
C
      IOOUT2   = 13
      OPEN(UNIT=IOOUT2,FILE=NAMEOUT2,FORM='FORMATTED',STATUS='UNKNOWN')
      WRITE(IOOUT2,'(4a12)')'X','A ax','A phi','dcp'
c
      NAMEOUT3 = 'steuer.dat'
C
      IOOUT3   = 14
      OPEN(UNIT=IOOUT3,FILE=NAMEOUT3,FORM='FORMATTED',STATUS='UNKNOWN')
      WRITE(IOOUT3,'(5a10,a5)')'fn','z','pitch','chord','aec','ipr'
c
      NAMEOUT4 = 'wtperf.out'
C
      IOOUT4   = 15
      OPEN(UNIT=IOOUT4,FILE=NAMEOUT4,FORM='FORMATTED',STATUS='UNKNOWN')
c
      aec    = 0.25
      ipr    = 0
      steuer = 'p1.dat'
c
C----------------------------------------------------------------------
c
      cptot = 0.0
      cttot = 0.0
c
      nx = 20
      dx = 1./nx
c
c     blade loop 
c
      DO IX=0,nx-1
c
c     x = r/rtip
c     xx = x * lambda
c
        X   = (0.5+ix)*dx
        XX  = LAMBDA*X
C
c       start from AD
C---    estimate from ACT DISK 
C
c	deVries Eq  4.4.18 on page 4-11
c
        THETANULL  = .5*ATAN(1./XX)
        SIGMANULL  = 4.0*(1.-COS(THETANULL))
C
        CCLACT     = 2.*PI*X*SIGMANULL/B
c        write(*,*)'x sigmanull ',x,sigmanull
        DCPALT     = DELTACP(X,SIGMANULL,THETANULL)
        DCTALT     = DELTACT(X,SIGMANULL,THETANULL)
C
c       WRITE(*,*)'from ACTUATOR DISK-THEORY:'
c       WRITE(*,*)'---------------------------'
c       WRITE(*,*)'THETANULL(GRAD): ',180.*THETANULL/PI
c       WRITE(*,*)'CCL              ',CCLACT
c       WRITE(*,*)'DCPALT:          ',DCPALT   
c       WRITE(*,*)'DCTALT:          ',DCTALT 
c       WRITE(*,*)'---------------------------'   
c       WRITE(*,*)
C
        DELTATHETA = .05*THETANULL
C
c       WRITE(*,*)'DELTATHETA(GRAD):',180.*DELTATHETA/PI
C
C----------------------------------------------------------------------------
C---    outer (theta) ITERATION
C       maximize (local !) NOT cP but from  Eq 36:
c
c       dcP-car = eta*(1 + v-wind/v-car)*dcP - dcT
c
C       WRITE(*,101)'THETA         DCP           DF          '
C101    FORMAT(19X,A42)
C
        J = 0
C
 10     DO WHILE(.TRUE..and.theta.lt.thetamax)
           J = J + 1
C
           IF (J.EQ.1)THEN
              THETA = THETANULL
              FALT  = DELTACPauto(X,SIGMANULL,THETANULL)
           ENDIF
C
C---------------------------------------------------------------------
C          INNer (SIGMA) Iteration to
C          check if constraint G(SIG,THETA) = 0 is fullfilled
c
c         from wilson-lissaman (deVries Eq. 4.4.14) 
c         may be some kind of  orthogonalty condition
c
           I = 0
              DO WHILE(.TRUE.)
                 I=I+1
                 IF (I.EQ.1)THEN
                     SIGMA = SIGMANULL
                 ELSE
 1                   SIGMA = SIGMA 
     +                - GRX(X,SIGMA,THETA)/DGDSIGMA(X,SIGMA,THETA)
                 ENDIF
C
                 GRXD = GRX(X,SIGMA,THETA)
C  
C---             stop INNER ITERATION:
C 
                 IF(ABS(GRXD).LT.EPSG.OR.I.GT.MAXI)EXIT 
              END DO   
              IF(I.GT.MAXI)WRITE(*,*)'WARNing: MAXI reached'
C
C---          END INNER ITERATION
C----------------------------------------------------------------------
C
C             WRITE(*,*)I,'.INN. ITER.: SIG,TH,GRX ',SIGMA,THETA,GRXD
C
              FNEU  = DELTACPauto(X,SIGMA,THETA)
              DF    = FNEU - FALT
              DFEPS = ABS((FNEU - FALT)/FNEU)
C
c             WRITE(*,100)J,'.AEUSS. ITER: THETA,DCP ',180.*THETA/PI,FNEU,DF
c100          FORMAT(I4,A14,3E14.5)
C
              IF (DF.GT.EPS)THEN
                 THETA = THETA + DELTATHETA
                 FALT  = FNEU
              ELSE
                 DELTATHETA = -0.5*DELTATHETA
                 THETA      = THETA + DELTATHETA
                 FALT       = FNEU
              ENDIF   
C
              IF (DFEPS.LT.EPS.OR.J.GT.MAXJ)EXIT
          END DO   
C
C---      END outer iteration
C-----------------------------------------------------------------------
C
c         WRITE(21,*)
c         WRITE(21,*)'-----------------------------------------------'
C         WRITE(21,*)'ENDE DER ITERATION: I,J',I,J
c         WRITE(21,*)'ENDE DER ITERATION: '
c         WRITE(21,*)'-----------------------------------------------'
C
          CCLOUT  =  CCL(X,SIGMA,THETA)
          AXOUT   =  AAX(X,SIGMA,THETA)
          APHIOUT = APHI(X,SIGMA,THETA)
C
c         WRITE(21,*)'CCL,THETA     ',CCLOUT,180.*THETA/PI
c         WRITE(21,*)'AAX,APHI      ',AXOUT,PHIOUT
C
C----------------------------------------------------------------------
C
c         compute increment of  cP and  cT bestimmen
c
          sth  = sin(theta)
          cth  = cos(theta)
          tth  = sth/cth
          la   = LAMBDA
          ftip = F(x,theta)
          AA   = AAX(X,SIGMA,THETA)*ftip
c
c--------------------------------------------------------------------------
c
c         de Vries Eq 4 from page 4.17
c         drag = 1/GZ
c
c         About factor 8 lambda:
c         Eq (1) on page 4.17 and X  = lambda x gives X² dX -> lambda³ dlambda

          dcp  = 8.*lambda*aa*(1.-aa)*(tth - drag)*x*x*dx
c
c         de Vries Eq 4.4.7 + 4.4.10 on page 4-10
c
          dct  = 4.* aa*(1.-aa)*(1. + tth*drag)*x*dx
c---------------------------------------------------------------------------
c
          cptot = cptot + dcp
          cttot = cttot + dct
c
C----    output local values
C
	  thetadeg = 180.*THETA/PI
          thetadeg = thetadeg - aoa
c
          WRITE(IOOUT1,'(7F12.4)')X,xx,CCLOUT,1.e3*X*Rtip,
     +          1.e3*CCLOUT*Rtip/CLdes,
     +          THETAdeg,ftip
c
          WRITE(*      ,'(9f12.4)')X,xx,CCLOUT,1.e3*X*Rtip,
     +          1.e3*CCLOUT*Rtip/CLdes,
     +          THETAdeg,ftip,dcp/dx,dct/dx
c
          WRITE(IOOUT2,'(4F12.6)')X, AXOUT,APHIOUT, dcp
c
c         output for blade desing (CAD) tool 
c
          WRITE(IOOUT3,'(a10,F10.4,f10.1,2F10.4,i5)')
     +      steuer,X*Rtip,thetadeg,cclout*rtip,aec,ipr
c
          WRITE(IOOUT4,'(3F10.4)')
     +      cclout*rtip,0.08*cclout*rtip,thetadeg
c
       END DO	
c
C---  END blade loop =  1, nx
C
c
       cpauto = etadr*(1.+1./vratio)*cptot - cttot
       aa     = 0.5*(1.-sqrt(1.-cttot))
       aaa    = 1. - cptot/cttot
c
       WRITE(*,*)
       Write(*,*)"############### Global values ################"
       WRITE(*,*)'ROTOR '
       WRITE(*, '(4a15  )')'cP','cT','a from cP/cT',' a from cT'
       WRITE(* ,'(4f15.3)') cptot,cttot,aaa,aa
       WRITE(*,*)'CAR '
       WRITE(*, '(a11,f11.3)')'cP-car/eta = ',cpauto
c------------------------------------------------------------------
c       WRITE(ioout1, '(a4,a8,   6a9  )')'B','LAMBDA','GZ',
c     +                               'cP','cT','cP-car/K','a1','a2'
c       WRITE(ioout1 ,'(i4,2f8.1,6f9.3)') int(B) , LAMBDA , GZ, 
c     +                                   cptot,cttot,cpauto,aa,aaa
c------------------------------------------------------------------
       CLOSE(IOOUT1)
       CLOSE(IOOUT2)
C
      END
C 
C------------------------------------------------------------------------
C------------------------------------------------------------------------
c
C----------------------------------------------
      FUNCTION AAX(X,SIG,THETA)
c
c     aax is 6b of eqs 6 a - d) on page 4-17
c    
C----------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      REAL LAMBDA
C
      AAX = SQRT(F(X,THETA)*F(X,THETA)+4.*A(SIG,THETA)*F(X,THETA)*
     1      (1.-F(X,THETA)))
      AAX = (2.*A(SIG,THETA)+F(X,THETA)-AAX)
     1      /(2.*(A(SIG,THETA)+F(X,THETA)**2))
C
c      WRITE(*,*)'debug: AAX= ',AAX
c
      END FUNCTION AAX
C
C------------------------------------------
      FUNCTION APHI(X,SIG,THETA) 
C------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      REAL LAMBDA
C
      APHI  = SIG/(4.*F(X,THETA)*COS(THETA) - SIG)
C
      IF(APHI.LT.0.)WRITE(21,*)'WARNING: APHI:' ,APHI
      END FUNCTION APHI
C
C------------------------------------------
      FUNCTION A(SI,TH) 
c
c     Eq 6d onpage 4.17
c     factor cL missing ?
c
C------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      REAL LAMBDA
C
      A   = SI* COS(TH)/(4.*SIN(TH)*SIN(TH))
C
      IF(A.LT.0)WRITE(21,*)'WARNING: A= ',A
      END FUNCTION A
C
C---------------------------------------------------------------------------
      FUNCTION F(X,THETA) 
C---------------------------------------------------------------------------
C
C     Tip-Loss model from L. Prandtl (1919)
c
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      REAL LAMBDA
C
c------------------------------------------------
      F = .5*B*(1.-X)/(X*SIN(THETA))
c-------------------------------------------------
c
      F = (2./PI) * ACOS(EXP(-F))
c
c      WRITE(*,*)'debug: F= ',F
      END FUNCTION F
c
C
C------------------------------------------------------------------------------
      FUNCTION CCL(X,SIG,THETA) 
C------------------------------------------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      REAL LAMBDA
C
      CCL = 2.*PI*X*SIG/B
C
C      WRITE(*,*)'CCL= ',CCL
      END FUNCTION CCL
C
C------------------------------------------
      FUNCTION GRX(X,SIG,THETA) 
C------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      REAL LAMBDA
C
      GRX=LAMBDA * X -
     +    AAX(X,SIG,THETA)*(1.-AAX(X,SIG,THETA)*F(X,THETA))*TAN(THETA)
     +    /(APHI(X,SIG,THETA)*(1.-AAX(X,SIG,THETA)))
C
C      WRITE(21,*)'GRX= ',GRX
C      WRITE(21,*)'L,X,AAX,TH,APHI',LAMBDA,X,AAX(X,SIG,THETA),THETA,
C     +                            APHI(X,SIG,THETA)
      END FUNCTION GRX
C
C--------------------------------------------------------------------
      FUNCTION DELTACPauto(X,S,T) 
C---------------------------------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      REAL LAMBDA
C
c     From Gaunaa, Oeye and Mikkelsen EWEC 2009 page 6 Eqs 30/31/36
c
       DELTACPauto = 
     + etadr*(1.+1./vratio)*deltacp(x,s,t)-deltact(x,s,t)
C
      END FUNCTION DELTACPauto
C
C---------------------------------------------------------------------
      FUNCTION DELTACT(X,S,T) 
C---------------------------------------------------------------------
c     Eqs 4.4.7 and 4.4.10 on page 4-10 
c
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      REAL LAMBDA
c
C------------------------------------------
      DELTACT = (TAN(T)/GZ + 1.)
     +          * 4.* AAX(X,S,T)*F(X,T)
     +          *(1.- AAX(X,S,T)*F(X,T))
C------------------------------------------
      END FUNCTION DELTACT
C
C------------------------------------------
      FUNCTION DELTACP(X,S,T) 
c------------------------------------------
c     Eq (6) on page 4-17
c
C------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      REAL LAMBDA
C
      DELTACP = 8.*lambda*(-1./GZ + TAN(T))
     +          *    AAX(X,S,T)*F(X,T)
     +          *(1.-AAX(X,S,T)*F(X,T))
C
      END FUNCTION DELTACP
C
C------------------------------------------
      FUNCTION DGDTHETA(X,S,T) 
C------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      DT = T/100.
      DGDTHETA   = (GRX(X,S,T+DT)-GRX(X,S,T))/DT
C
C      WRITE(21,*)'DGDTHETA= ',DGDTHETA
      END FUNCTION DGDTHETA
C
C------------------------------------------
      FUNCTION DGDSIGMA(X,S,T) 
C------------------------------------------
C
      COMMON PI,B,LAMBDA,GZ,vratio,etadr
C
      DS = S/100.
      DGDSIGMA    = (GRX(X,S+DS,T)-GRX(X,S,T))/DS
C
C      WRITE(21,*)'DGDSIGMA= ',DGDSIGMA
      END FUNCTION DGDSIGMA
C
C------------------------------------------------------------------
