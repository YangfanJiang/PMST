density = 4.933  # g/cm^3
mass = 126.90  # g/mol

##Jump Factors
JK = 6.039
JL1 = 1.149
JL2 = 1.400
JL3 = 2.950

##Absorption Edges
KEDGE = 33.169  # Unit KeV K edge energy
L1EDGE = 5.188
L2EDGE = 4.852
L3EDGE = 4.557

##fluo yields
fyK = 0.884
fyL1 = 0.044
fyL2 = 0.079
fyL3 = 0.079
## fluo transition probabilities
kA1 = 0.5260
kA2 = 0.2830
kB1 = 0.1010
kB2 = 0.0380
kB3 = 0.0520
l1G3 = 0.4140
l1B3 = 0.4830
l1B4 = 0.1030
l2B1 = 0.8650
l2G1 = 0.1350
l3L = 0.0320
l3A2 = 0.1100
l3A1 = 0.7480
l3B2 = 0.1100

##fluo transition energies
fluoTranEnergies = {
'Kalpha1': 26.612,
'Kalpha2': 28.318,
'Kbeta1' : 32.295,
'Kbeta3': 32.240,
'Kbeta2': 33.054,

'L1gamma3': 5.072,
'L1beta3': 4.313,
'L1beta4': 4.257,

'L2beta1': 4.221,
'L2gamma1': 4.801,

'L3l': 3.485,
'L3alpha2': 3.926,
'L3alpha1': 3.938,
"L3beta2": 4.509}


##non-radiative transitions
## Original Vacancy K shell
kL1L1 = 0.078
kL1L2 = 0.094
kL1L3 = 0.115
kL2L2 = 0.012
kL2L3 = 0.250
kL3L3 = 0.126
kL1X = 0.087
kL2X = 0.079
kL3X = 0.134
kXX = 0.025
## L1 shell
l1L2X = 0.188
l1L3X = 0.293
l1XX = 0.519
## L2 shell
l2L3X = 0.167
l2XX = 0.833