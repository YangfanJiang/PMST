density = 1.873  # g/cm^3
mass = 132.91  # g/mol

##Jump Factors
JK = 5.949
JL1 = 1.150
JL2 = 1.400
JL3 = 2.955

##Absorption Edges
KEDGE = 35.985  # Unit KeV K edge energy
L1EDGE = 5.714
L2EDGE = 5.359
L3EDGE = 5.012

##fluo yields
fyK = 0.898
fyL1 = 0.049
fyL2 = 0.090
fyL3 = 0.091
## fluo transition probabilities
kA1 = 0.5240
kA2 = 0.2830
kB1 = 0.1020
kB2 = 0.0390
kB3 = 0.0520
l1G3 = 0.4220
l1B3 = 0.4720
l1B4 = 0.1060
l2B1 = 0.8310
l2G1 = 0.1690
l3L = 0.0270
l3A2 = 0.1170
l3A1 = 0.6230
l3B2 = 0.2330

##fluo transition energies
fluoTranEnergies = {
'Kalpha1': 30.973,
'Kalpha2': 30.625,
'Kbeta1' : 34.987,
'Kbeta3': 34.920,
'Kbeta2': 35.833,

'L1gamma3': 5.567,
'L1beta3': 4.717,
'L1beta4': 4.652,

'L2beta1': 4.620,
'L2gamma1': 5.280,

'L3l': 3.795,
'L3alpha2': 4.272,
'L3alpha1': 4.286,
"L3beta2": 4.936}


##non-radiative transitions
## Original Vacancy K shell
kL1L1 = 0.080
kL1L2 = 0.098
kL1L3 = 0.115
kL2L2 = 0.012
kL2L3 = 0.245
kL3L3 = 0.123
kL1X = 0.092
kL2X = 0.082
kL3X = 0.136
kXX = 0.017
## L1 shell
l1L2X = 0.120
l1L3X = 0.294
l1XX = 0.586
## L2 shell
l2L3X = 0.169
l2XX = 0.831