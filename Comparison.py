import random
import numpy as np
## Comparison method for W
def compareShell(K, L1, L2, L3):
        r = random.random()
        A = K
        B = A + L1
        C = B + L2
        D = C + L3
        # E = D + M
        if r <= A :
            return 'K'
        if r > A and r <= B:
            return 'L1'
        if r > B and r <= C:
            return 'L2'
        if r > C and r <=D:
            return 'L3'
        else:
            return False

def compareAttenuation(pPE, pCS, pRS):
        r = random.uniform(0, 1)
        A = pPE
        B = A + pCS
        C = B + pRS
        # D = C + self.L3
        if r <= A :
            # print('r', r)
            return 'PE' # Photoelectric Effect
        if r > A and r <= B:
            return 'CS' # Compton Scattering
        if r > B:
            return 'RS' # Reyleigh Scattering
        # if r > C:
        #     return 'P' # passing through


def compareKtransition(A1, A2, B1, B3, B2):
    r = random.random()
    A = A1
    B = A + A2
    C = B + B1
    D = C + B3
    E = D + B2
    if r <= A:
        return 'Kalpha1'
    if r > A and r <= B:
        return 'Kalpha2'
    if r > B and r <= C:
        return 'Kbeta1'
    if r > C and r <= D:
        return 'Kbeta3'
    if r > D:
        return 'Kbeta2'

def compareL1transition(G3, B3, B4):
    r = random.random()
    A = G3
    B = A + B3
    C = B + B4
    if r <= A:
        return 'L1gamma3'
    if r > A and r <= B:
        return 'L1beta3'
    if r > B :
        return 'L1beta4'

def compareL2transition(B1, G1):
    r = random.random()
    A = B1
    B = A + G1
    if r <= A:
        return 'L2beta1'
    if r > A:
        return 'L2gamma1'

def compareL3transition(L, A2, A1, B2):
    r = random.random()
    A = L
    B = A + A2
    C = B + A1
    D = C + B2
    if r <= A:
        return 'L3l'
    if r > A and r <=B:
        return 'L3alpha2'
    if r > B and r <= C:
        return 'L3alpha1'
    if r > C:
        return "L3beta2"



def compareNon_RaKtransition(KL1L1, KL1L2, KL1L3, KL2L2, KL2L3, KL3L3, KL1X, KL2X, KL3X, KXX):
    r = random.random()
    A = KL1L1
    B = A + KL1L2
    C = B + KL1L3
    D = C + KL2L2
    E = D + KL2L3
    F = E + KL3L3
    G = F + KL1X
    H = G + KL2X
    I = H + KL3X
    J = I + KXX

    if r <= A:
        return 'KL1L1'
    if r > A and r <= B:
        return 'KL1L2'
    if r > B and r <= C:
        return 'KL1L3'
    if r > C and r <= D:
        return 'KL2L2'
    if r > D and r <= E:
        return 'KL2L3'
    if r > E and r <= F:
        return 'KL3L3'
    if r > F and r <= G:
        return 'KL1X'
    if r > G and r <= H:
        return 'KL2X'
    if r > H and r <= I:
        return 'KL3X'
    if r > I:
        return 'KXX'

def compareNon_RaL1transition(L1L2X, L1L3X, L1XX):
    r = random.random()
    A = L1L2X
    B = A + L1L3X
    C = B + L1XX
    if r <= A:
        return 'L1L2X'
    if r > A and r <= B:
        return 'L1L3X'
    if r > B :
        return 'L1XX'

def compareNon_RaL2transition(L2L3X, L2XX):
    r = random.random()
    A = L2L3X
    B = A + L2XX
    if r <= A:
        return 'L2L3X'
    if r > A:
        return 'L2XX'

def compareNon_RaL3transition():
    return "L3XX"

def comparecsAngle(angles, dist):
    # rangeLs = []
    # s = 0
    # rangeLs.append(s)
    # for i in dist:
    #     s += i
    #     rangeLs.append(s)
    # # maxarg = np.argmax(np.array(dist))
    rangeLs = np.cumsum(dist)
    rangeLs = rangeLs / max(rangeLs)
    # test = np.add.accumulate(dist)
    #while r < rangeLs[0] or r > rangeLs[-1]:
    r = random.random()
    for j in range(len(rangeLs)):
        if r <=rangeLs[j]:
            return angles[j]

