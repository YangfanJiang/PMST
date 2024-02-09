# coefficients for W74
import numpy as np

KEDGE = 69.53 # Unit KeV K edge energy
L1EDGE = 12.1
L2EDGE = 11.54
L3EDGE = 10.21
MINENERGY = 1 #min energy
# MAXENERGY = 10000 #max energy
MAXENERGY = 1000 #max energy

KTOTALWITHCO = 11.23 # cm^2/g
L1TOTALWITHCO = 238.2 #the point after K edge
L2TOTALWITHCO = 231.2
L3TOTALWITHCO = 233.4
KBETOTALWITHCO = 2.552# cm^2/g the point before K edge (or entering K edge)
L1BETOTALWITHCO = 206.5
L2BETOTALWITHCO = 168.9
L3BETOTALWITHCO = 92
BEGTOTALWITHCO = 3683 #starting coefficient for total attentuatio with coherent scaterring
# ENDTOTALWITHCO = 0.04747 #end coefficient for total attentuatio with coherent scaterring for 10000 keV
ENDTOTALWITHCO = 0.06618 #end coefficient for total attentuatio with coherent scaterring for 1000 keV

KPE = 10.8 # cm^2/g
L1PE = 234.4
L2PE = 227.3
L3PE= 229
KBEPE = 2.118 # cm^2/g
L1BEPE = 202.8
L2BEPE = 165
L3BEPE= 87.59
BEGPE = 3672 #starting coefficient for photoelectric effect
# ENDPE = 3.747*10**(-4) #end coefficient for photoelectric effect for 10000 keV
ENDPE = 0.01283 #end coefficient for photoelectric effect for 1000 keV

KCS = 0.1022 # cm^2/g
L1CS = 0.05441
L2CS = 0.0528
L3CS= 0.04861
BEGCS = 0.00434 #starting coefficient for Compton Scaterring
# ENDCS = 0.0124 #end coefficient for Compton Scaterring for 10000 keV
ENDCS = 0.05087 #end coefficient for Compton Scaterring for 1000 keV

KRS = 0.3318 # cm^2/g
L1RS = 3.675
L2RS = 3.859
L3RS= 4.363
BEGRS = 11.44 #starting coefficient for Rayleigh (Coherent) Scaterring
# ENDRS = 2.533*10**(-5) #end coefficient for Rayleigh Scaterring for 10000keV
ENDRS = 2.477*10**(-3) #end coefficient for Rayleigh Scaterring for 1000keV
def interpolate(energy):
    # from max to K-edge
    if energy >= MINENERGY and energy < L3EDGE:
        peCo = BEGPE + (energy - MINENERGY)*(L3BEPE - BEGPE)/(L3EDGE - MINENERGY)
        csCo = BEGCS + (energy - MINENERGY) * (L3CS - BEGCS) / (L3EDGE - MINENERGY)
        rsCO = BEGRS + (energy - MINENERGY) * (L3RS - BEGRS) / (L3EDGE - MINENERGY)
        totalCo = BEGTOTALWITHCO + (energy - MINENERGY)*(L3BETOTALWITHCO - BEGTOTALWITHCO)/(L3EDGE - MINENERGY)
    elif energy >= L3EDGE and energy < L2EDGE:
        peCo = L3PE + (energy - L3EDGE) * (L2BEPE - L3PE) / (L2EDGE - L3EDGE)
        csCo = L3CS + (energy - L3EDGE) * (L2CS - L3CS) / (L2EDGE - L3EDGE)
        rsCo = L3RS + (energy - L3EDGE) * (L2RS - L3RS) / (L2EDGE - L3EDGE)
        totalCo = L3TOTALWITHCO + (energy - L3EDGE) * (L2BETOTALWITHCO - L3TOTALWITHCO) / (L2EDGE - L3EDGE)
    elif energy >= L2EDGE and energy < L1EDGE:
        peCo = L2PE + (energy - L2EDGE) * (L1BEPE - L2PE) / (L1EDGE - L2EDGE)
        csCo = L2CS + (energy - L2EDGE) * (L1CS - L2CS) / (L1EDGE - L2EDGE)
        rsCo = L2RS + (energy - L2EDGE) * (L1RS - L2RS) / (L1EDGE - L2EDGE)
        totalCo = L2TOTALWITHCO + (energy - L2EDGE) * (L1BETOTALWITHCO - L2TOTALWITHCO) / (L1EDGE - L2EDGE)
    elif energy >= L1EDGE and energy < KEDGE:
        peCo = L1PE + (energy - L1EDGE) * (KBEPE - L1PE) / (KEDGE - L1EDGE)
        csCo = L1CS + (energy - L1EDGE) * (KCS - L1CS) / (KEDGE - L1EDGE)
        rsCo = L1RS + (energy - L1EDGE) * (KRS - L1RS) / (KEDGE - L1EDGE)
        totalCo = L1TOTALWITHCO + (energy - L1EDGE) * (KBETOTALWITHCO - L1TOTALWITHCO) / (KEDGE - L1EDGE)
    elif energy >= KEDGE:
        peCo = KPE + (energy - KEDGE) * (ENDPE - KPE) / (MAXENERGY - KEDGE)
        csCo = KCS + (energy - KEDGE) * (ENDCS - KCS) / (MAXENERGY - KEDGE)
        rsCo = KRS + (energy - KEDGE) * (ENDRS - KRS) / (MAXENERGY - KEDGE)
        totalCo = KTOTALWITHCO + (energy - KEDGE) * (ENDTOTALWITHCO - KTOTALWITHCO) / (MAXENERGY - KEDGE)
    return peCo, csCo, rsCo, totalCo

# re = interpolate(50)
# for co in re:
#     print(co)