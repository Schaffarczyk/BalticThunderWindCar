import numpy as np

NO1 = "VP-BTS-cp-max.dat"
NO2 = "dummy"
PI = 3.1415926
rho = 1.25
g = 9.81
R1 = 0.0
R2 = 0.95
R = R2
A = PI * (R1**2 + R2**2)
btt = False

if btt:
    print(f"{'**** BTT **** ':<15}{'Rotor-Area = ':<15}{A:6.2f}{' m**2':<5}")
    a = 4.0

lamax = 5.24
cpmax = 0.40
aa = 3.0 * cpmax / lamax**2
bb = 2.0 * cpmax / lamax**3
ks = [-5.386e-3, -6.454, 2.468, -7.007e-2]
kct = [-3.642e-1, 3.571e-1, -4.200e-2, 1.471e-3]
upmmax = 1500.0
pomax = 20000.0
mradmax = 500.0
rrad = 0.0256 * 26.0 / 2.0
effdr = 0.85
mass = 285.0
croll = 0.01
cw = 0.62
as_ = 0.4
optflag = True
cpmaxflag = False
gearflag = False
gear = [9.0, 8.0, 7.3, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 1.0]
ngmax = 9
gearup = 1.1
laopt = 5.24
cpopt = 0.4
ctopt = 0.57

IO1 = open(NO1, 'w')
IO2 = open(NO2, 'w')
vwa = 2.0
vwe = 12.0
dvw = 0.25
nv = int((vwe - vwa) / dvw)

if gearflag:
    print("driving on gear-mode")
    for i in range(1, ngmax + 1):
        print(f"{'gear = ':<10}{i:<3}{'Ueberset.: ':<15}{gear[i-1]:6.2f}")

if cpmaxflag:
    ngear = 0
    print("driving on cp-max")

if optflag:
    ngear = 0
    print("driving on P-rest max")

print()
header = ['v-wi(m/s)', 'upm-r', 'ro-P(W)', 'car-P(W)', 'schub(N)', 'Rollw', 'Drag', 'lambda', 'cP', 'cT', 'M-rad(Nm)', 'v-au(m/s)', 'ratio ', 'Uebers.', ' Gang']
print("\t".join(header))
IO1.write("\t".join(header) + "\n")

for i in range(nv + 1):
    pautoalt = 0.0
    pautoneu = 0.0
    porotalt = 0.0
    porotneu = 0.0
    Vw = vwa + i * dvw
    vautomin = 0.05 * Vw
    vautomax = 15.0
    vauto = vautomin
    dvauto = 0.05 * vautomin
    
    while True:
        vauto += dvauto
        if vauto >= vautomax:
            print("***** vc Max")
            break

        vwrel = Vw + vauto
        if gearflag:
            ngear = 1
            while True:
                upma = 30.0 * vauto / (PI * rrad)
                upm = gear[ngear - 1] * upma
                om = upm * PI / 30.0
                la = om * R / vwrel
                cp = aa * la**2 - bb * la**3
                if la > gearup * laopt and ngear <= ngmax:
                    ngear += 1
                else:
                    break
            
            ct = sum([kct[j] * la**j for j in range(4)])
            ct = max(ct, 0.0)

        elif cpmaxflag:
            om = lamax * vwrel / R
            upm = 30.0 * om / PI
            la = om * R / vwrel
            cp = aa * la**2 - bb * la**3
            ct = sum([kct[j] * la**j for j in range(4)])
            ct = max(ct, 0.0)

        elif optflag:
            la = laopt
            om = la * vwrel / R
            upm = 30.0 * om / PI
            cp = cpopt
            ct = ctopt

        if upm > upmmax:
            upm = upmmax
            om = PI * upm / 30.0
            la = om * R / vwrel
            print("***** upm Max")

        porotalt = porotneu
        porotneu = effdr * cp * 0.5 * rho * A * vwrel**3
        th = ct * 0.5 * rho * A * vwrel**2

        if porotneu > pomax:
            print("***** Power Max")
            porotneu = pomax
            cp = porotneu / (0.5 * rho * A * vwrel**3)
            maxiit = 300
            iit = 0
            laneu = la
            laalt = laneu
            la = laneu + 0.05
            while True:
                la += 0.01
                iit += 1
                if iit > maxiit:
                    print("***** Iterations Max")
                    break
                laalt = laneu
                laneu = la
                la = (laalt + laneu) / 2.0
                cpit = aa * la**2 - bb * la**3
                ct = sum([kct[j] * la**j for j in range(4)])
                ct = max(ct, 0.0)
                if abs(cpit - cp) < 1e-6:
                    break

            om = la * vwrel / R
            upm = 30.0 * om / PI

        mrad = 0.5 * rho * A * vwrel**2 * sum([ks[j] * la**j for j in range(4)])
        mrad = min(mrad, mradmax)
        porad = om * mrad
        pautoalt = pautoneu
        pautoneu = porotneu + porad
        pauto = mradmax * om + cw * 0.5 * rho * A * vauto**2 + croll * mass * g
        acc = (pautoneu - pauto) / (mass * vauto)

        if acc > 0.0:
            vautomax = vauto + (porotneu - porad) / (cw * 0.5 * rho * A * vauto**2)
            vautomax = min(vautomax, 15.0)
            dvauto = 0.2
        else:
            dvauto = 0.05

        if vauto >= vautomax:
            break

        # Output
        row = [Vw, upm, porotneu, pauto, th, croll * mass * g, cw * 0.5 * rho * A * vauto**2, la, cp, ct, mrad, vauto, (pautoneu - porad) / pauto, gear[ngear - 1], ngear]
        print("\t".join(map(str, row)))
        IO1.write("\t".join(map(str, row)) + "\n")

IO1.close()
IO2.close()
