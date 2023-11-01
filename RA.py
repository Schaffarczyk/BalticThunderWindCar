#
#  python version of R(acing)A(eolus) - a blade design code for Head Wind Driven Cars
#
#
#     PROGRAMM RA  (C)  A.P. SCHAFFARCZYK Kiel University of Applied Scieces
#     based on versions from
#     march 1998, nov 2004, oct 2009, 23. august 2011/2019 (chinook), 2023

#     basic equations are from deVRIES AGARDOgraph 243 (1979), chapter 4.4 ff
#     and  
#     Gaunaa, Oeye and Mikkelsen, proc. ewec 2009, Equation  (36)  
#

import math

# Global variables

# numerical constants

PI       = 4. * math.atan(1.)
thetamax = PI / 2.

EPSall = 1.E-06
EPSF = EPSall
EPSG = EPSall
EPST = EPSall

EPSA = 1.E-03
relax= 1.

MAXI = 100
MAXJ = 100

#   nx number of blade - sections
nx    = 20

# parameters for car

B      = 3.0
# design TSR
LAMBDA = 5.
#tip radius (m)
Rtip   = 0.95

# airfoil data
GZ     = 80.0
drag = 1./GZ
#design lift(-coefficient) cL-des
CLdes  = 1.00
# angle of attach at cL-des
aoa    = 6.

# assumed velocity ratio of car
vratio = 0.5

#assumed drive train efficiency
etadr = 0.88


#  convert TSR-des to string (for use in file-names)

intla = int(10*LAMBDA)
i1 = intla//10
i2 = intla-10*i1
s1 = chr(48+i1)
s2 = chr(48+i2)

# declaration of functions

def shen():
    gc = 0.1 + math.exp(-0.125*(B*LAMBDA-21.))
    return gc

def F(X, THETA):
#    print(f"enter F")
    gc = shen()
    gc = 1.
    F_val = gc * 0.5 * B * (1.-X) / (X*math.sin(THETA))
    if F_val > 0.:
       F_val = (2./PI) * math.acos(math.exp(-F_val))
#       print(f"F= {F_val}")
    else:
       F_val = 1.
    return F_val

def APHI(X, SIG, THETA):
    APHI_val = SIG / (4.*F(X, THETA)*math.cos(THETA) - SIG)
    if APHI_val < 0:
#        print(f"wrong APHI= {APHI_val}")
        APHI_val = EPSA
    return APHI_val

def A(SI, TH):
    A_val = SI * math.cos(TH) / (4.*math.sin(TH)*math.sin(TH))
    if A_val < 0:
#        print(f"wrong A= {A_val}")
        A_val = EPSA
    return A_val

def AAX(X, SIG, THETA):
    AAX_val = math.sqrt(F(X,THETA)*F(X,THETA)+4.*A(SIG,THETA)*F(X,THETA)*(1.-F(X,THETA)))
    AAX_val = (2.*A(SIG,THETA)+F(X,THETA)-AAX_val)/(2.*(A(SIG,THETA)+F(X,THETA)**2))
    return AAX_val  


def CCL(X, SIG, THETA):
    return 2.*PI*X*SIG/B

def GRX(X, S, T):
    go = AAX (X,S,T)*(1.-AAX(X, S, T)*F(X, T))*math.tan(T)
    gu = APHI(X,S,T)*(1.-AAX(X,S,T))
    return LAMBDA*X - go/gu

def DELTACT(X, S, T):
    xx = LAMBDA*X
    return (math.tan(T)/GZ + 1.) * AAX(X, S, T)*F(X, T) * (1.-AAX(X, S, T)*F(X, T)) / xx

def DELTACP(X, S, T):
    return (-1./GZ+math.tan(T)) * AAX(X, S, T)*F(X, T) * (1.-AAX(X, S, T)*F(X, T))

def DELTACPauto(X, S, T):
    return etadr*(1.+1./vratio)*DELTACP(X, S, T) - DELTACT(X, S, T)

def DGDTHETA(X, S, T):
    DT = T/100.
    return (GRX(X, S, T+DT) - GRX(X, S, T)) / DT

def DGDSIGMA(X, S, T):
    DS   = S/100.
    DG   = GRX(X, S+DS, T) - GRX(X, S, T)
    if DG < EPSF:
        DG = EPSF
    return DG/DS

#********************************************************
# main iteration loops:
# (1)  sections
# (2)  outer = via theta (flow angle): maximum of cP-car by some kind of binary search
# (3)  inner = via sigma (silidity)  : induction factors a and a' by Newton's method


def main():
    print(f"***** START python code *****")
    print(f" ")
    print(f"vra= {vratio:.3f} eta= {etadr:.3f}")
    print(f"TSR= {LAMBDA:.0f} L2D= {GZ:.0f} cL-des= {CLdes:.2}")
    print(f" ")
    cptot = 0.0
    cttot = 0.0

    dx    = 1./nx

    with open('rotor.txt', 'w') as outfile:
        print        (f"r/rtip  r/mm    chord TH(°)  a       aphi  F\n")
        outfile.write(f"r/rtip  r/mm    chord TH(°)  a       aphi  F\n")
 
#       section loop 1 .. nx

        for IX in range(nx):
            X  = (0.5+IX)*dx
            XX = LAMBDA*X
            THETANULL = 0.5*math.atan(1./XX)

#            print(f"IX X= {IX} {X:.3f}")
#            print(f"XX= {XX:.3f}")
#            print(f"THETANULL= {THETANULL:.3f}")

#           estimate from simple AD Theory

            SIGMANULL = 4.0*(1.-math.cos(THETANULL))
            CCLACT    = 2.0*PI*X*SIGMANULL/B
            DCPALT    = DELTACP(X, SIGMANULL, THETANULL)
            DCTALT    = DELTACT(X, SIGMANULL, THETANULL)
            DELTATHETA = 0.1*THETANULL

            Jouter= 0
 
            THETA = THETANULL
            SIGMA = SIGMANULL
            FALT  = DELTACPauto(X,SIGMA,THETA)
            FNEU  = 0.
            DFEPS = 1.

#           outer (THETA) Iteration loop
#
#           maximize
#           cP-car = etadr*(1 + v-wind/v-car)*cP -cT
#

            while THETA <= thetamax and Jouter <= MAXJ and DFEPS > EPST:
#                print(f"outer itno theta DFEPS {Jouter} {THETA:2f} {DFEPS:5f}")
                Jouter += 1
                if Jouter == 1:
                    THETA = THETANULL
                    FALT  = DELTACPauto(X, SIGMA, THETA)

#-----------------------------------------------------------------------------------
#           inner (sigma, solidity) iteration by binary search
#----------------------------------------------------------------------------------
# 
                Jinner= 0
                while abs(GRX(X,SIGMA,THETA)) > EPSG and Jinner <= MAXI:
                    Jinner += 1
                    dSIGMA = -relax*GRX(X,SIGMA,THETA)/DGDSIGMA(X,SIGMA,THETA)
                    SIGMA += dSIGMA

#                    print(f"inner GRX =  {Jinner} {abs(GRX(X,SIGMA,THETA)):.5f}")
#                    print(f"inner s ds= {SIGMA:.5f} {dSIGMA:.5f}  ")
#                    print(f"")

#               end inner iteration

                FNEU  = DELTACPauto(X,SIGMA,THETA)
                DF    = FNEU - FALT
                DFEPS = abs(DF/FNEU)

#               beyond maximum if fneu < falt
#               go back then

                if DF < 0.:
#                    print(f" beyond max jounter= {Jouter}")
                    DELTATHETA = -0.5*DELTATHETA

                THETA += DELTATHETA
                FALT  = FNEU


#           end outer (theta) iteration

            CCL_out  = CCL (X,SIGMA,THETA)
            AX_out   = AAX (X,SIGMA,THETA)
            APHI_out = APHI(X,SIGMA,THETA)

#           increments of cP and cT

            sth  = math.sin(THETA)
            cth  = math.cos(THETA)
            tth  = sth/cth
            la   = LAMBDA
            ftip = F(X,THETA)
            AA   = AAX(X,SIGMA,THETA)

            DCP = AA*(1.-AA)*(tth  -   drag)*X*X*dx
            DCT = AA*(1.-AA)*(1. + tth*drag)  *X*dx

            rr    = 1.e3*X*Rtip
            chout = 1.e3*CCL_out*Rtip/CLdes
            THdeg = THETA*180./PI

#           Writing output data
#                           r/rtip   r/mm     chord TH(°)    theta         a              aphi

            print        (f"{X:.4f}  {rr:5.1f} {chout:5.1f}  {THdeg:6.3f}  {AX_out:.3f}   {APHI_out:.3f} {ftip:5.3f} ")
            outfile.write(f"{X:.4f}  {rr:5.1f} {chout:5.1f}  {THdeg:6.3f}  {AX_out:.3f}   {APHI_out:.3f} {ftip:5.3f}\n")

            if X <= Rtip:
                cptot += DCP
                cttot += DCT

#----------------------------------------------------------------------------------------------------------
#   normalization

    cptot  = 8.*la*cptot
    cttot  = 8.*   cttot
    cpauto = etadr*(1.+ 1./vratio)*cptot - cttot

    aaux   = 1. - cttot
    aaux   = 1. -math.sqrt(aaux)
    aa     = 0.5*aaux
     
# Print-out total CP and CT values
    print(f"")
    print(f"***** integral values *****")
    print(f"cP= {cptot:.3f} cT= {cttot:.3f} cP-car/eta-dr= {cpauto:.3f} a= {aa:.3f}")
    print(f"")
    print(f"***** END *****")

#outfile.write(f"cttot= {cttot:.3f} cptot= {cptot:.3f}")
 
#   Interal python
if __name__ == "__main__":
    main()
