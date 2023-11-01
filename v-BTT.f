c
C---------------------------------------------------------------------
C     PROGRAMM  v-auto (C)  P. SCHAFFARCZYK  02/2009 ... 08/2011
C     gibt  v-auto aus
c     aerodyn - weber blaetter (2012 - 2014)
c     BTT V1          22 aug  2015
c         V2 (Anfang) 09 sept 2016  
c                        feb  2018
c	                 aug  2018
c     BTS redux sten     aug  2023
c--------------------------------------------------------------------
C
      CHARACTER*30 NO1, NO2
      real om,  a, lamax, la, mrad, laneu, laalt, laopt
      real mradmax, mass
c
      real ks(0:3), kct(0:3), gear(1:10)
C
      logical gearflag, cpmaxflag, optflag, BTT
c
      NO2 = 'dummy'
C
c     zahlen
c
      PI     = 3.1415926
      rho    = 1.25
      g      = 9.81
c
c     Rotordaten
c     Radius inkl diffusor
c
      R1      = 0.0
      R2      = 0.95
c
c     bezugs radius
      R       = r2
c     
      A      = pi*(r1*r1+r2*r2)
c
c  ************************
c
      btt = .false.
c
      if (btt)then
        write(*,'(2a15,f6.2,a5)') "**** BTT **** ","Rotor-Area = ",
     +  A," m**2"
        a = 4.  
      endif	
c
c    SW (aerodyn) Rotor
c    einfache empirische rotorkennlinie nach manwell
c    BEM Daten in aerodyn-SW1-MW-cP-cT.xls
c
      lamax = 5.24
      cpmax = 0.40
c
      aa = 3.*cpmax/lamax**2
      bb = 2.*cpmax/lamax**3
c
c    koeffizienten schub (N) vs v (m/s) aerodyn wt_perf
c
      ks(0) = -5.386e-3
      ks(1) = -6.454
      ks(2) =  2.468
      ks(3) = -7.007e-2
c
c    ct (tsr) BEM 
c    
      kct(0) = -3.642e-1
      kct(1) =  3.571e-1
      kct(2) = -4.200e-2
      kct(3) =  1.471e-3
c
c     Maxwerte Rotor vergl. data sheet 3WTR WindDynamic 
c
c     1200 zu klein
c
      upmmax  =  1500.
      pomax   = 20000.
c
c     maximales zul. moment am rad
c
      mradmax =  500.
c
c    car konstants
c=================
c    radius antriebsrad 
c
      rrad   = 0.0256*26./2.
c
c   triebstrangeff
c
c      effdr   =   0.86 th
      effdr   =   0.85
c
c   masse BTT (inkl. fahrer)
c
c    wiegen: 23.08.18: 210 + 75 
c
      mass = 285.
c
c     Masse BTS 2015: mass    = 320
c     
c     konst. Rollwiderstand 25 N, aus Messung inkl. Getriebe
c
      croll  =   0.01
c                0.007
c
c   cW wert verkleidung (ink. turm), stirnflaeche
c
c     APS:                       0.62
c     Müller Studienarbeit 2015: 0.25
c
      cw = 0.62
c
c     streamlined ?
c
c      cw = 0.1
c      
      as = 0.4
c        
c    3 FAHR modi moeglich (nur einzeln !)
c
c    konstanter gang   : gearflag
c    cpmax geregelt    : cpmaxflag 
c    opt unklar -> max: optflag 
c
c
       optflag   = .true.
       cpmaxflag = .false.
       gearflag  = .false.
c
c     gears BTS 
c
c     Den Helder 22. 08. 2023
c
       gear(1)  = 9.
       gear(2)  = 8.
       gear(3)  = 7.3
       gear(4)  = 6.
       gear(5)  = 5.
       gear(6)  = 4.
       gear(7)  = 3.
       gear(8)  = 2.
       gear(9)  = 1.
       gear(10) = 1.0
c
c      max anzahl gänge      
       ngmax    = 9
c
c       
c
       gearup   = 1.1
c
c***** Daten fuer opt-mode: P-rest = P-rotor - v-car*thrust) -> max
c      Aerodyn-blaetter-wtperf-pitchen.xls
c
c************************************************************************
c
          laopt = 5.24
          cpopt = 0.4
          ctopt = 0.57
c	  
c        pitch 0°
c        aps 29 08 18 sehr unklar, vergl. xls von oben
c*********************************************************************
c
c    open output file(s) 
c
      NO1 = 'VP-BTS-cp-max.dat'
c
      IO1  = 11 	
      OPEN(UNIT=IO1,FILE=NO1,FORM='FORMATTED',STATUS='UNKNOWN') 
      IO2  = 12	
      OPEN(UNIT=IO2,FILE=NO2,FORM='FORMATTED',STATUS='UNKNOWN') 
c
	vwa =    2.0
        vwe =   12.0
	dvw =    0.25
c
	nv = (vwe-vwa)/dvw
c
c	write(*,*)'debug 1'
        if (gearflag)then
           write(*,109)'driving on gear-mode'
           do i =1,ngmax
               write(*,108)'gear = ',i,'Ueberset.:', gear(i)
           end do 
        endif
        if (cpmaxflag)then
           ngear = 0
           write(*,109)'driving on cp-max'
        endif 
        if (optflag)then
           ngear = 0
           write(*,109)'driving on P-rest max'
        endif   
c	write(*,*)'debug 2'
C
	WRITE(*,*)
        write(*     ,107)'v-wi(m/s)','upm-r','ro-P(W)','car-P(W)',
     +   'schub(N)','Rollw','Drag',
     +   'lambda','cP','cT',
     +   'M-rad(Nm)','v-au(m/s)','ratio ','Uebers.',' Gang'
c
        write(11     ,107)'v-wi(m/s)','upm-r','ro-P(W)','car-P(W)',
     +   'schub(N)','Rollw','Drag',
     +   'lambda','cP','cT',
     +   'M-rad(Nm)','v-au(m/s)','ratio ','Uebers.',' Gang' 
c
c------------------------------------------------------------------------------
c
c    loop windgeschwindigkeiten
c
      do i=0,nv
c
         pautoalt = 0.0
         pautoneu = 0.0
         porotalt = 0.0
         porotneu = 0.0
c
         Vw    = vwa + i*dvw
c
c         write(*,*)'***** v-wind: ',vw
c
c       stepping v-auto 
c
         vautomin =  0.05*vw
         vautomax = 15.0
         vauto    = vautomin
         dvauto   =  0.05*vautomin
c
c       erhoehe autogeschwindigkeit bis leistung aus dem wind erschoepft ist
c
 801     vauto = vauto + dvauto
         if(vauto.ge.vautomax)then
            write(*,*) '***** vc Max'
            goto 802
	 endif
c
c         write(*,*)'vauto ',vauto   
         vwrel = vw + vauto 	
c
c********************************************************************************
c ************** Schaltmodus ****************************************************
c
c       beginne im ersten gang
c
         if(gearflag)then
	    ngear = 1
555         upma = 30.*vauto/(pi*rrad)
            upm  = gear(ngear)*upma
            om   = upm*pi/30.
            la   = om*r/vwrel  
            cp   = aa*la**2-bb*la**3
c
c           hochschalten - falls rotordrehzahl (gearup-1)*100 % zu gross
c***************************************************************************
c
c	    if (la.gt.lamax.or.ngear.gt.ngmax) then
	    if (la.gt.gearup*laopt.and.ngear.le.ngmax) then	 
               ngear = ngear+1
	       goto 555
            endif
c
            ct = 0.
            do j=0,3
               ct = ct+kct(j)*la**j
            end do   
            if (ct.lt.0)ct = 0.
         endif
c
c        end gear modus
c***************************************************************************
c
c***************************************************************************
c*****        cP - max *****************************************************
c
         if(cpmaxflag)then
	    om  = lamax*vwrel/r
	    upm = 30.*om/pi
            la  = om*r/vwrel  
            cp  = aa*la**2-bb*la**3
c
            ct = 0.
            do j=0,3
               ct = ct+kct(j)*la**j
            end do  
            if (ct.lt.0)ct = 0.
c	    write(*,*)'debug ','la cT ',la, ct 
         endif   
c
c*****************************************************************
c*****      P-rest-> max "optimum"  ******************************
c
c      emp korrelation; verlg. xls Auswertung Windkanal 27 07 11
c          
c
         if(optflag)then 
            la    = laopt
	    om    = la*vwrel/r
	    upm   = 30.*om/pi
            cp    = cpopt
            ct    = ctopt
         endif   
c
c	 drehzahl begrenzen wenn zu gross
c
	 if (upm.gt.upmmax)then
            upm = upmmax
            om  = pi*upm/30.
            la  = om*r/vwrel
	    write(*,*) '***** upm Max'
         endif
c
c       leistung und schub absolut bestimmen
c
         porotalt = porotneu
         porotneu = effdr*cp*0.5*rho*a*vwrel**3
c
         th = ct*0.5*rho*a*vwrel**2
c
c****** drehzahl abregeln wenn p > pmax
c
	 if (porotneu.gt.pomax)then
	    write(*,*)'***** Power Max'
	    porotneu  = pomax
	    cp        = porotneu/(0.5*rho*a*vwrel**3)
c
c        iterative suche nach neuem lambda
c
           maxiit = 300
           iit    = 0
           laneu  = la
c
 902       iit    = iit+1
           laalt = laneu
c
c          manwell aufgeloest nach lambda
c
	   laneu = sqrt(cp/(aa-bb*laalt))
c
           if(abs(laalt-laneu).lt.01.or.iit.ge.maxiit)then
             goto 900
           else
             goto 902
           endif
c   
 900       la = laneu
c
	   upm = 30*vwrel*la/(pi*r)
c
           ct = 0.
           do k=0,3
              ct = ct + kct(k)*la**k
              if (ct.lt.0)   ct = 0.
              if (ct.gt.1.5) ct = 1.5
           end do   
c
c           neuen schub und ct bestimmen
c
            th = 0.
            do k=0,3
               th = th + ks(k)*vwrel**k
            enddo
c
            ct = th/(0.5*rho*a*vwrel**2)
c
	 endif
c
c******* drehzahl abregeln wenn p > pmax
c
c       Rollwiderstand  
c
         r0 = mass*g*croll
c
c       bestimme leistungsbedarf bei v-auto
c
         drag  = cw*0.5*rho*vwrel**2*as
	 reibk = r0 + drag
	 fauto = reibk + th
c
         pautoalt = pautoneu
 	 pautoneu = fauto*vauto
c
         rotormom = pautoalt/om
c
c     begin debug
c2345
         write(io2,'(a30,2f8.2,5f8.2)')'wi ca alt PAa PAn PRa PRn UPM ',
     +              vw,vauto,pautoalt,pautoneu,porotalt,porotneu,upm
c     end debug
c
         ueb = 30.*vauto/(pi*rrad)
         ueb = upm/ueb
c
         quo  = vauto/vw
         mrad = rrad*fauto
c
c       wenn  moment am rad zu gross: abbrechen
c
         if(mrad.gt.mradmax)then
	    write(*,*)'***** Mo Max'
            goto 800
         end if
c   
c	auto faehrt bei v-wind mit v-auto wenn erstmals pauto < porot
c
         if (pautoalt.gt.porotalt.and.pauotneu.lt.porotneu)goto 800
c
c       wenn nicht, dann mit erhoehtem v-auto weiter iterieren
c
	 goto 801
c
c	datenausgabe
c
800     ngear = 0	
	if(ueb.gt.0.9*gear(10).and.ueb.lt.1.1*gear(10))ngear = 10
	if(ueb.gt.0.9*gear( 9).and.ueb.lt.1.1*gear( 9))ngear = 9
	if(ueb.gt.0.9*gear( 8).and.ueb.lt.1.1*gear( 8))ngear = 8
	if(ueb.gt.0.9*gear( 7).and.ueb.lt.1.1*gear( 7))ngear = 7
	if(ueb.gt.0.9*gear( 6).and.ueb.lt.1.1*gear( 6))ngear = 6
	if(ueb.gt.0.9*gear( 5).and.ueb.lt.1.1*gear( 5))ngear = 5
	if(ueb.gt.0.9*gear( 4).and.ueb.lt.1.1*gear( 4))ngear = 4
	if(ueb.gt.0.9*gear( 3).and.ueb.lt.1.1*gear( 3))ngear = 3
	if(ueb.gt.0.9*gear( 2).and.ueb.lt.1.1*gear( 2))ngear = 2
	if(ueb.gt.0.9*gear( 1).and.ueb.lt.1.1*gear( 1))ngear = 1
c
	pcar = porotalt-vauto*th
c
        write(IO1,106)vw,upm,porotalt,pcar,th,r0,drag,
     +                la,cp,ct, mrad,vauto,quo,ueb,ngear
         write(* ,106)vw,upm,porotalt,pcar,th,r0,drag,
     +                la,cp,ct,mrad,vauto,quo,ueb,ngear
c
c       naechste windgeschw.
c
	write(io2,*)
 802  continue
c
      enddo
c 
c    end loop windgeschw
c
      close(UNIT=IOOUT) 
c
 103  format(a15,6f6.1)
 104  format(f6.2,a20)
 105  format(a10,f8.3)
c
 106  format(13f10.3,f10.2,i10)
 107  format(15a10)
c
 108    format(a10,i3,a15,f6.2)
 109    format(a30)
c
      END






