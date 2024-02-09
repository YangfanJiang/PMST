import xraydb, xraylib
import vedo
import numpy as np
import ast
from contextlib import suppress

RES = 120 #120
class Absorber: # testing with W (tungsten)

        def __init__(self, atomNum, property_file=None):

            self.atomNum = atomNum
            self.mass = xraylib.AtomicWeight(atomNum)
            self.pos = None # distance to source
            self.x, self.y, self.z = None, None, None
            self.sep = 0
            self.density = xraylib.ElementDensity(self.atomNum)
            self.w = None
            self.l = None
            self.hole = None
            self.r = None
            self.t = None
            self.pts = []
            self.topPts = []
            self.botPts = []
            self.energies = None
            self.csCoLs = None
            self.rsCoLs = None
            self.peCoLs = None
            self.totalCoLs = None
            self.ffs = None
            self.hole_pos = []

            if property_file == None:

                ##Jump Factors
                self.JK = xraylib.JumpFactor(atomNum, xraylib.K_SHELL)
                self.JL1 = xraylib.JumpFactor(atomNum, xraylib.L1_SHELL)
                self.JL2 = xraylib.JumpFactor(atomNum, xraylib.L2_SHELL)
                self.JL3 = xraylib.JumpFactor(atomNum, xraylib.L3_SHELL)
                # self.JM1 = xraylib.JumpFactor(atomNum, xraylib.M1_SHELL)

                ##Absorption Edges
                self.KEDGE = xraylib.EdgeEnergy(atomNum, xraylib.K_SHELL)
                self.L1EDGE = xraylib.EdgeEnergy(atomNum, xraylib.L1_SHELL)
                self.L2EDGE = xraylib.EdgeEnergy(atomNum, xraylib.L2_SHELL)
                self.L3EDGE = xraylib.EdgeEnergy(atomNum, xraylib.L3_SHELL)
                # self.M1EDGE = xraylib.EdgeEnergy(atomNum, xraylib.M1_SHELL)

                ##fluo yields
                self.fyK = xraylib.FluorYield(atomNum, xraylib.K_SHELL)
                self.fyL1 = xraylib.FluorYield(atomNum, xraylib.L1_SHELL)
                self.fyL2 = xraylib.FluorYield(atomNum, xraylib.L2_SHELL)
                self.fyL3 = xraylib.FluorYield(atomNum, xraylib.L3_SHELL)

                fluor_Ka1_energy,  self.kA1 = xraydb.fluor_yield(atomNum, 'K', 'Ka1', 1)[1] / 1000, \
                                                   xraydb.fluor_yield(atomNum, 'K', 'Ka1', 1)[2]

                fluor_Ka2_energy, self.kA2 = xraydb.fluor_yield(atomNum, 'K', 'Ka2', 1)[1] / 1000, \
                                             xraydb.fluor_yield(atomNum, 'K', 'Ka2', 1)[2]

                fluor_Kb1_energy, self.kB1 = xraydb.fluor_yield(atomNum, 'K', 'Kb1', 1)[1] / 1000, \
                                                   xraydb.fluor_yield(atomNum, 'K', 'Kb1', 1)[2]

                fluor_Kb2_energy, self.kB2 = xraydb.fluor_yield(atomNum, 'K', 'Kb2', 1)[1] / 1000, \
                                                   xraydb.fluor_yield(atomNum, 'K', 'Kb2', 1)[2]

                fluor_Kb3_energy, self.kB3 = xraydb.fluor_yield(atomNum, 'K', 'Kb3', 1)[1] / 1000, \
                                             xraydb.fluor_yield(atomNum, 'K', 'Kb3', 1)[2]
                sum = self.kA1+ self.kA2+ self.kB1+ self.kB2+ self.kB3
                self.kA1,  self.kA2, self.kB1, self.kB2, self.kB3 = self.kA1/sum,  self.kA2/sum, self.kB1/sum, self.kB2/sum, self.kB3/sum


                fluor_L1g3_energy, self.l1G3 = xraydb.fluor_yield(atomNum, 'L1', 'Lg3', 1)[1] / 1000, \
                                             xraydb.fluor_yield(atomNum, 'L1', 'Lg3', 1)[2]

                fluor_L1b3_energy, self.l1B3 = xraydb.fluor_yield(atomNum, 'L1', 'Lb3', 1)[1] / 1000, \
                                               xraydb.fluor_yield(atomNum, 'L1', 'Lb3', 1)[2]

                fluor_L1b4_energy, self.l1B4 = xraydb.fluor_yield(atomNum, 'L1', 'Lb4', 1)[1] / 1000, \
                                               xraydb.fluor_yield(atomNum, 'L1', 'Lb4', 1)[2]

                sum = self.l1G3 + self.l1B3 + self.l1B4
                self.l1G3, self.l1B3, self.l1B4 = self.l1G3 / sum, self.l1B3 / sum, self.l1B4 / sum

                fluor_L2b1_energy, self.l2B1 = xraydb.fluor_yield(atomNum, 'L2', 'Lb1', 1)[1] / 1000, \
                                               xraydb.fluor_yield(atomNum, 'L2', 'Lb1', 1)[2]

                fluor_L2g1_energy, self.l2G1 = xraydb.fluor_yield(atomNum, 'L2', 'Lg1', 1)[1] / 1000, \
                                               xraydb.fluor_yield(atomNum, 'L2', 'Lg1', 1)[2]

                sum = self.l2B1 + self.l2G1

                self.l2B1, self.l2G1 =  self.l2B1/sum, self.l2G1/sum

                fluor_L3l_energy, self.l3L = xraydb.fluor_yield(atomNum, 'L3', 'Ll', 1)[1] / 1000, \
                                               xraydb.fluor_yield(atomNum, 'L3', 'Ll', 1)[2]

                fluor_L3a2_energy, self.l3A2 = xraydb.fluor_yield(atomNum, 'L3', 'La2', 1)[1] / 1000, \
                                             xraydb.fluor_yield(atomNum, 'L3', 'La2', 1)[2]

                fluor_L3a1_energy, self.l3A1 = xraydb.fluor_yield(atomNum, 'L3', 'La1', 1)[1] / 1000, \
                                               xraydb.fluor_yield(atomNum, 'L3', 'La1', 1)[2]

                fluor_L3b2_energy, self.l3B2 = xraydb.fluor_yield(atomNum, 'L3', 'Lb2', 1)[1] / 1000, \
                                               xraydb.fluor_yield(atomNum, 'L3', 'Lb2', 1)[2]

                sum = self.l3L+self.l3A2+self.l3A1+self.l3B2

                self.l3L, self.l3A2, self.l3A1, self.l3B2 =  self.l3L/sum, self.l3A2/sum, self.l3A1/sum, self.l3B2/sum


                if fluor_Ka1_energy == 0:
                    self.kA1 = 0

                if fluor_Ka2_energy == 0:
                    self.kA2 = 0

                if fluor_Kb1_energy == 0:
                    self.kB1 = 0

                if fluor_Kb2_energy == 0:
                    self.kB2 = 0

                if fluor_Kb3_energy == 0:
                    self.kB3 = 0

                if fluor_L1g3_energy == 0:
                    self.l1G3 = 0

                if fluor_L1b3_energy == 0:
                    self.l1B3 = 0

                if fluor_L1b4_energy == 0:
                    self.l1B4 = 0

                if fluor_L2b1_energy == 0:
                    self.l2B1 = 0

                if fluor_L2g1_energy == 0:
                    self.l2G1 = 0

                if fluor_L3l_energy == 0:
                    self.l3L = 0

                if fluor_L3a2_energy == 0:
                    self.l3A2 = 0

                if fluor_L3a1_energy == 0:
                    self.l3A1 = 0

                if fluor_L3b2_energy == 0:
                    self.l3B2 = 0






                ## fluo transition probabilities
                # self.kA1 = xraydb.fluor_yield(atomNum, 'K', 'Ka1', 1)[2]
                # self.kA2 = xraydb.fluor_yield(atomNum, 'K', 'Ka2', 1)[2]
                # self.kB1 = xraydb.fluor_yield(atomNum, 'K', 'Kb1', 1)[2]
                # self.kB2 = xraydb.fluor_yield(atomNum, 'K', 'Kb2', 1)[2]
                # self.kB3 = xraydb.fluor_yield(atomNum, 'K', 'Kb3', 1)[2]
                # self.l1G3 = xraydb.fluor_yield(atomNum, 'L1', 'Lg3', 1)[2]
                # self.l1B3 = xraydb.fluor_yield(atomNum, 'L1', 'Lb3', 1)[2]
                # self.l1B4 = xraydb.fluor_yield(atomNum, 'L1', 'Lb4', 1)[2]
                # self.l2B1 = xraydb.fluor_yield(atomNum, 'L2', 'Lb1', 1)[2]
                # self.l2G1 = xraydb.fluor_yield(atomNum, 'L2', 'Lg1', 1)[2]
                # self.l3L = xraydb.fluor_yield(atomNum, 'L3', 'Ll', 1)[2]
                # self.l3A2 = xraydb.fluor_yield(atomNum, 'L3', 'La2', 1)[2]
                # self.l3A1 = xraydb.fluor_yield(atomNum, 'L3', 'La1', 1)[2]
                # self.l3B2 = xraydb.fluor_yield(atomNum, 'L3', 'Lb2', 1)[2]




                ##fluo transition energies
                self.fluoTranEnergies = {
                    'Kalpha1': fluor_Ka1_energy, #xraydb.fluor_yield(atomNum, 'K', 'Ka1', 1)[1]/1000,
                    'Kalpha2': fluor_Ka2_energy, #xraydb.fluor_yield(atomNum, 'K', 'Ka2', 1)[1]/1000,
                    'Kbeta1': fluor_Kb1_energy, #xraydb.fluor_yield(atomNum, 'K', 'Kb1', 1)[1]/1000,
                    'Kbeta3': fluor_Kb3_energy, #xraydb.fluor_yield(atomNum, 'K', 'Kb3', 1)[1]/1000,
                    'Kbeta2': fluor_Kb2_energy, #xraydb.fluor_yield(atomNum, 'K', 'Kb2', 1)[1]/1000,

                    'L1gamma3': fluor_L1g3_energy, #xraydb.fluor_yield(atomNum, 'L1', 'Lg3', 1)[1]/1000,
                    'L1beta3': fluor_L1b3_energy, #xraydb.fluor_yield(atomNum, 'L1', 'Lb3', 1)[1]/1000,
                    'L1beta4': fluor_L1b4_energy, #xraydb.fluor_yield(atomNum, 'L1', 'Lb4', 1)[1]/1000,

                    'L2beta1': fluor_L2b1_energy, #xraydb.fluor_yield(atomNum, 'L2', 'Lb1', 1)[1]/1000,
                    'L2gamma1': fluor_L2g1_energy, #xraydb.fluor_yield(atomNum, 'L2', 'Lg1', 1)[1]/1000,

                    'L3l': fluor_L3l_energy, #xraydb.fluor_yield(atomNum, 'L3', 'Ll', 1)[1]/1000,
                    'L3alpha2': fluor_L3a2_energy, #xraydb.fluor_yield(atomNum, 'L3', 'La2', 1)[1]/1000,
                    'L3alpha1': fluor_L3a1_energy, #xraydb.fluor_yield(atomNum, 'L3', 'La1', 1)[1]/1000,
                    "L3beta2": fluor_L3b2_energy}#xraydb.fluor_yield(atomNum, 'L3', 'Lb2', 1)[1]/1000}

                # self.KAlpha1 = 'Kalpha1', 59.318
                # self.KAlpha2 = 'Kalpha2', 57.982
                # self.KBeta1 = 'Kbeta1', 67.245
                # self.KBeta3 = 'Kbeta3', 66.952
                # self.KBeta2 = 'Kbeta2', 69.067
                #
                # self.L1amma3 = 'L1gamma3', 11.675
                # self.L1Beta3 = 'L1beta3', 9.819
                # self.L1Beta4 = 'L1beta4', 9.526
                #
                # self.L2Beta1 = 'L2beta1', 9.672
                # self.L2Gamma1 = 'L2gamma1', 11.286
                #
                # self.L3L = 'L3l', 7.388
                # self.L3Alpha2 = 'L3alpha2', 8.335
                # self.L3Alpha1 = 'L3alpha1', 8.398
                # self.L3Beta2 = "L3beta2", 9.962

                ##non-radiative transitions
                ## Original Vacancy K shell inilisation
                self.kL1L1 = 0
                self.kL1L2 = 0
                self.kL1L3 = 0
                self.kL2L2 = 0
                self.kL2L3 = 0
                self.kL3L3 = 0
                self.kL1M1 = 0
                self.kL1M2 = 0
                self.kL1M3 = 0
                self.kL1M4 = 0
                self.kL1M5 = 0
                self.kL1N1 = 0
                self.kL1N2 = 0
                self.kL1N3 = 0
                self.kL1N4 =0
                self.kL1N5 = 0
                self.kL1N6 = 0
                self.kL1N7 = 0

                self.kL2M1 = 0
                self.kL2M2 = 0
                self.kL2M3 = 0
                self.kL2M4 = 0
                self.kL2M5 = 0
                self.kL2N1 = 0
                self.kL2N2 = 0
                self.kL2N3 = 0
                self.kL2N4 = 0
                self.kL2N5 = 0
                self.kL2N6 = 0
                self.kL2N7 = 0

                self.kL3M1 = 0
                self.kL3M2 =0
                self.kL3M3 = 0
                self.kL3M4 = 0
                self.kL3M5 = 0
                self.kL3N1 = 0
                self.kL3N2 = 0
                self.kL3N3 = 0
                self.kL3N4 =0
                self.kL3N5 = 0
                self.kL3N6 = 0
                self.kL3N7 = 0

                ## L1 shell initilisation
                self.l1L2X = 0
                self.l1L3X = 0

                ## L2 shell initilisation
                self.l2L3X =0

                # try:
                try:
                    self.kL1L1 = xraylib.AugerRate(atomNum, xraylib.K_L1L1_AUGER)
                except:
                    pass
                try:
                    self.kL1L2 = xraylib.AugerRate(atomNum, xraylib.K_L1L2_AUGER)
                except:
                    pass
                try:
                    self.kL1L3 = xraylib.AugerRate(atomNum, xraylib.K_L1L3_AUGER)
                except:
                    pass
                    self.kL2L2 = xraylib.AugerRate(atomNum, xraylib.K_L2L2_AUGER)
                try:
                    self.kL2L3 = xraylib.AugerRate(atomNum, xraylib.K_L2L3_AUGER)
                except:
                    pass
                try:
                    self.kL3L3 = xraylib.AugerRate(atomNum, xraylib.K_L3L3_AUGER)
                except:
                    pass
                try:
                    self.kL1M1 = xraylib.AugerRate(atomNum, xraylib.K_L1M1_AUGER)
                except:
                    pass
                try:
                    self.kL1M2 = xraylib.AugerRate(atomNum, xraylib.K_L1M2_AUGER)
                except:
                    pass
                try:
                    self.kL1M3 = xraylib.AugerRate(atomNum, xraylib.K_L1M3_AUGER)
                except:
                    pass
                try:
                    self.kL1M4 = xraylib.AugerRate(atomNum, xraylib.K_L1M4_AUGER)
                except:
                    pass
                try:
                    self.kL1M5 = xraylib.AugerRate(atomNum, xraylib.K_L1M5_AUGER)
                except:
                    pass
                try:
                    self.kL1N1 = xraylib.AugerRate(atomNum, xraylib.K_L1N1_AUGER)
                except:
                    pass
                try:
                    self.kL1N2 = xraylib.AugerRate(atomNum, xraylib.K_L1N2_AUGER)
                except:
                    pass
                try:
                    self.kL1N3 = xraylib.AugerRate(atomNum, xraylib.K_L1N3_AUGER)
                except:
                    pass
                try:
                    self.kL1N4 = xraylib.AugerRate(atomNum, xraylib.K_L1N4_AUGER)
                except:
                    pass
                try:
                    self.kL1N5 = xraylib.AugerRate(atomNum, xraylib.K_L1N5_AUGER)
                except:
                    pass
                try:
                    self.kL1N6 = xraylib.AugerRate(atomNum, xraylib.K_L1N6_AUGER)
                except:
                    pass
                try:
                    self.kL1N7 = xraylib.AugerRate(atomNum, xraylib.K_L1N7_AUGER)
                except:
                    pass
                self.kL1X = (self.kL1M1 + self.kL1M2 + self.kL1M3 + self.kL1M4 + self.kL1M5
                             + self.kL1N1 + self.kL1N2 + self.kL1N3 + self.kL1N4 + self.kL1N5+self.kL1N6 + self.kL1N7)
                try:
                    self.kL2M1 = xraylib.AugerRate(atomNum, xraylib.K_L2M1_AUGER)
                except:
                    pass
                try:
                    self.kL2M2 = xraylib.AugerRate(atomNum, xraylib.K_L2M2_AUGER)
                except:
                    pass
                try:
                    self.kL2M3 = xraylib.AugerRate(atomNum, xraylib.K_L2M3_AUGER)
                except:
                    pass
                try:
                    self.kL2M4 = xraylib.AugerRate(atomNum, xraylib.K_L2M4_AUGER)
                except:
                    pass
                try:
                    self.kL2M5 = xraylib.AugerRate(atomNum, xraylib.K_L2M5_AUGER)
                except:
                    pass
                try:
                    self.kL2N1 = xraylib.AugerRate(atomNum, xraylib.K_L2N1_AUGER)
                except:
                    pass
                try:
                    self.kL2N2 = xraylib.AugerRate(atomNum, xraylib.K_L2N2_AUGER)
                except:
                    pass
                try:
                    self.kL2N3 = xraylib.AugerRate(atomNum, xraylib.K_L2N3_AUGER)
                except:
                    pass
                try:
                    self.kL2N4 = xraylib.AugerRate(atomNum, xraylib.K_L2N4_AUGER)
                except:
                    pass
                try:
                    self.kL2N5 = xraylib.AugerRate(atomNum, xraylib.K_L2N5_AUGER)
                except:
                    pass
                try:
                    self.kL2N6 = xraylib.AugerRate(atomNum, xraylib.K_L2N6_AUGER)
                except:
                    pass
                try:
                    self.kL2N7 = xraylib.AugerRate(atomNum, xraylib.K_L2N7_AUGER)
                except:
                    pass
                self.kL2X = (self.kL2M1 + self.kL2M2 + self.kL2M3 + self.kL2M4 + self.kL2M5
                                 + self.kL2N1 + self.kL2N2 + self.kL2N3 + self.kL2N4 + self.kL2N5 + self.kL2N6 + self.kL2N7)
                try:
                    self.kL3M1 = xraylib.AugerRate(atomNum, xraylib.K_L3M1_AUGER)
                except:
                    pass
                try:
                    self.kL3M2 = xraylib.AugerRate(atomNum, xraylib.K_L3M2_AUGER)
                except:
                    pass
                try:
                    self.kL3M3 = xraylib.AugerRate(atomNum, xraylib.K_L3M3_AUGER)
                except:
                    pass
                try:
                    self.kL3M4 = xraylib.AugerRate(atomNum, xraylib.K_L3M4_AUGER)
                except:
                    pass
                try:
                    self.kL3M5 = xraylib.AugerRate(atomNum, xraylib.K_L3M5_AUGER)
                except:
                    pass
                try:
                    self.kL3N1 = xraylib.AugerRate(atomNum, xraylib.K_L3N1_AUGER)
                except:
                    pass
                try:
                    self.kL3N2 = xraylib.AugerRate(atomNum, xraylib.K_L3N2_AUGER)
                except:
                    pass
                try:
                    self.kL3N3 = xraylib.AugerRate(atomNum, xraylib.K_L3N3_AUGER)
                except:
                    pass
                try:
                    self.kL3N4 = xraylib.AugerRate(atomNum, xraylib.K_L3N4_AUGER)
                except:
                    pass
                try:
                    self.kL3N5 = xraylib.AugerRate(atomNum, xraylib.K_L3N5_AUGER)
                except:
                    pass
                try:
                    self.kL3N6 = xraylib.AugerRate(atomNum, xraylib.K_L3N6_AUGER)
                except:
                    pass
                try:
                    self.kL3N7 = xraylib.AugerRate(atomNum, xraylib.K_L3N7_AUGER)
                except:
                    pass
                self.kL3X = (self.kL3M1 + self.kL3M2 + self.kL3M3 + self.kL3M4 + self.kL3M5
                                 + self.kL3N1 + self.kL3N2 + self.kL3N3 + self.kL3N4 + self.kL3N5 + self.kL3N6 + self.kL3N7)
                self.kXX = 1-(self.kL1L1 + self.kL1L2+self.kL1L3+self.kL2L2+self.kL2L3+self.kL3L3+self.kL1X+self.kL2X+self.kL3X)
                    ## L1 shell
                try:
                    self.l1L2X = xraylib.CosKronTransProb(atomNum, xraylib.FL12_TRANS) # obtained from xraylib online calculator, original 0.2
                except:
                    pass
                try:
                    self.l1L3X = xraylib.CosKronTransProb(atomNum, xraylib.FL13_TRANS)  # obtained from xraylib online calculator, original 0.382
                except:
                    pass
                self.l1XX = 1 - (self.l1L2X  + self.l1L3X)  # obtained from xraylib online calculator, original 0.472
                    ## L2 shell
                try:
                    self.l2L3X = xraylib.CosKronTransProb(atomNum, xraylib.FL23_TRANS)
                except:
                    pass
                self.l2XX = 1 - self.l2L3X
            else:
                import importlib
                property_file = importlib.import_module(property_file)
                # import property_file
                self.density = property_file.density
                self.mass = property_file.mass
                ##Jump Factors
                self.JK = property_file.JK
                self.JL1 = property_file.JL1
                self.JL2 = property_file.JL2
                self.JL3 = property_file.JL3
                # self.JM1 = property_file.JM1

                ##Absorption Edges
                self.KEDGE = property_file.KEDGE  # Unit KeV K edge energy
                self.L1EDGE = property_file.L1EDGE
                self.L2EDGE = property_file.L2EDGE
                self.L3EDGE = property_file.L3EDGE
                # self.M1EDGE = property_file.M1EDGE

                ##fluo yields
                self.fyK = property_file.fyK
                self.fyL1 = property_file.fyL1
                self.fyL2 = property_file.fyL2
                self.fyL3 = property_file.fyL3

                ## fluo transition probabilities
                self.kA1 = property_file.kA1
                self.kA2 = property_file.kA2
                self.kB1 = property_file.kB1
                self.kB2 = property_file.kB2
                self.kB3 = property_file.kB3
                self.l1G3 = property_file.l1G3
                self.l1B3 = property_file.l1B3
                self.l1B4 = property_file.l1B4
                self.l2B1 = property_file.l2B1
                self.l2G1 = property_file.l2G1
                self.l3L = property_file.l3L
                self.l3A2 = property_file.l3A2
                self.l3A1 = property_file.l3A1
                self.l3B2 = property_file.l3B2

                ##fluo transition energies
                self.fluoTranEnergies = property_file.fluoTranEnergies #{
                    # 'Kalpha1': 59.318,
                    # 'Kalpha2': 57.982,
                    # 'Kbeta1': 67.245,
                    # 'Kbeta3': 66.952,
                    # 'Kbeta2': 69.067,
                    #
                    # 'L1gamma3': 11.675,
                    # 'L1beta3': 9.819,
                    # 'L1beta4': 9.526,
                    #
                    # 'L2beta1': 9.672,
                    # 'L2gamma1': 11.286,
                    #
                    # 'L3l': 7.388,
                    # 'L3alpha2': 8.335,
                    # 'L3alpha1': 8.398,
                    # "L3beta2": 9.962}

                ##non-radiative transitions
                ## Original Vacancy K shell
                self.kL1L1 = property_file.kL1L1
                self.kL1L2 = property_file.kL1L2
                self.kL1L3 = property_file.kL1L3
                self.kL2L2 = property_file.kL2L2
                self.kL2L3 = property_file.kL2L3
                self.kL3L3 = property_file.kL3L3
                self.kL1X = property_file.kL1X
                self.kL2X = property_file.kL2X
                self.kL3X = property_file.kL3X
                self.kXX = property_file.kXX
                ## L1 shell
                self.l1L2X = property_file.l1L2X
                self.l1L3X = property_file.l1L3X
                self.l1XX = property_file.l1XX
                ## L2 shell
                self.l2L3X = property_file.l2L3X
                self.l2XX = property_file.l2XX

        # def setLength(self, length):
        #     self.length = length

        def getLength(self):
            return self.l
        # def setWidth(self, width):
        #     self.width = width
        def getWidth(self):
            return self.w
        def setDensity(self):
            self.density = xraylib.ElementDensity(self.atomNum)
        def getDensity(self):
            return self.density
        def setThickness(self, t): # 0.001-10mm
            self.t = t
        def setRadius(self, r):
            self.r = r
        def setPositon(self, pos): # initialise at z=0 position
            self.pos = pos
            self.x, self.y, self.z = pos
        def setType(self, type):
            self.type = type
        def setSepration(self, seperation):
            self.sep = seperation
        def getSepration(self):
            return self.sep
        # def setHoleRePosition(self, hole_relative_z):
        #     self.hole_relative_z = hole_relative_z
        # def getHoleRePosition(self):
        #     return self.hole_relative_z
        def getType(self):
            return self.type
        def setLength(self, l): # need to be confirmed
            self.l = l
        def setWidth(self, w): # need to be confirmed
            self.w = w

        def getFormFactor(self, q):
            # xraydb.f0(self.atomNum, q)
            # f0 = self.offest
            # for s, e in zip(self.scale,
            #                 self.exponents):
            #     f0 += s * np.exp(-e * q * q)
            return xraydb.f0(self.atomNum, q)

        # def setAngle(self, angle): #knife edged pinhole angle of the horizontal plane and one side of the pinhole, pi/2 for half clinder hole
        #     self.angle = angle
        # def getAngle(self):
        #     return self.angle
        # def getVolume(self):
        #     if self.type == "cylinder":
        #         self.voulume = ParametricRegion((r*cos(theta), r*sin(theta), z), (theta, 0, 2*pi), (z, 0, self.t), (r, 0, self.sep/2)),

        # def getBoundaries(self):
        #     self.boundaries = []
        #     if self.type == 'cylinder':
        #         self.rightBoundary = Line(Point(self.sep/2, 0), Point(self.sep, self.t)) # x-z plane
        #         self.leftBoundary = Line(Point(-self.sep/2, 0), Point(-self.sep, self.t)) # x-z plane
        #     if self.type == 'knife': # for x-z plane
        #         self.rightUpperBoundary = Line(Point(self.sep/2, self.t/2), Point(self.r, self.t))
        #         self.leftUpperBoundary = Line(Point(-self.sep/2, self.t/2), Point(-self.r, self.t))
        #         self.rightLowerBoundary = Line(Point(self.r, 0), Point(self.sep/2, self.t/2))
        #         self.leftLowerBoundary = Line(Point(-self.r, 0), Point(-self.sep/2, self.t/2))

        # def getEntryPlane(self):
        #     # self.planes = []
        #     # self.upperPlane = vedo.Plane(pos=(0,0, z+self.t), normal=(0, 0, -1), sx=self.l, sy=self.w)
        #     self.entryPlane = vedo.Plane(pos=(0, 0, self.z), normal=(0, 0, 1), sx=self.length, sy=self.length)
        #     # self.planes.append(self.upperPlane)
        #     # self.planes.append(self.lowerPlane)
        #     return self.entryPlane

        # def getPlanes(self, z):
        #     self.planes = []
        #     self.upperPlane = vedo.Plane(pos=(0, 0, z + self.t), normal=(0, 0, -1), sx=self.l, sy=self.w)
        #     self.lowerPlane = vedo.Plane(pos=(0, 0, z), normal=(0, 0, 1), sx=self.l, sy=self.w)
        #     self.planes.append(self.upperPlane)
        #     self.planes.append(self.lowerPlane)
        #     return self.planes

        # def getPlanes(self): # for z coordinates
        #     # self.planesDict = {}
        #     self.planes = []
        #     self.upperPlane = Plane(Point3D(0,0,self.t), normal_vector=(0,0,1))
        #     self.lowerPlane = Plane(Point3D(0, 0, 0), normal_vector=(0, 0, 1))
        #     if self.type == 'cylinder':
        #         self.cylinder = Eq(x**2+y**2, self.sep/2)
        #     # if self.type == 'cylinder':
        #     #     self.rightPlane = Plane(Point3D(0.5*self.sep, 0, 0), normal_vector=(1,0,0))
        #     #     self.leftPlane = Plane(Point3D(-0.5*self.sep,0,0), normal_vector=(1,0,0))
        #     #     self.planes.append(self.leftPlane)
        #     #     self.planes.append(self.rightPlane)
        #     # elif self.type == 'knife':
        #     #     self.rightUpper = Plane(Point3D(0.5 * self.sep, 0, 0.5*self.t), normal_vector=(np.sin(self.angle), 0, np.cos(self.angle)))
        #     #     self.leftUpper = Plane(Point3D(-0.5 * self.sep, 0, 0.5*self.t), normal_vector=(np.sin(-self.angle), 0, np.cos(-self.angle)))
        #     #     self.rightLower = Plane(Point3D(0.5 * self.sep, 0, 0.5*self.t), normal_vector=(np.sin(-self.angle), 0, np.cos(-self.angle)))
        #     #     self.leftLower = Plane(Point3D(-0.5 * self.sep, 0, 0.5*self.t), normal_vector=(np.sin(self.angle), 0, np.cos(self.angle)))
        #     #     self.planes.append(self.rightUpper)
        #     #     self.planes.append(self.leftUpper)
        #     #     self.planes.append(self.rightLower)
        #     #     self.planes.append(self.leftLower)
        #     self.planes.append(self.upperPlane)
        #     self.planes.append(self.lowerPlane)
        #     return self.planes
        def getBox(self, color = 'green'):
            # x, y, z = self.pos
            box = vedo.Box(pos=self.pos,length=self.l, width=self.w, height=self.t, c= color)
            return box
        def getCone(self, up, color = 'green'):
            x, y, z = self.pos
            # height = self.t
            if up:
                cone = vedo.Cone(alpha=0.5, r=self.r, pos=(x, y, z), res=RES, height=self.t, axis=(0, 0, 1), c=color)
            else:
                cone = vedo.Cone(alpha=0.5, r=self.r, pos=(x, y, z), res=RES, height=self.t, axis=(0, 0, -1), c=color)

            return cone
        def getCones(self, pos = None, color='yellow'): # knife-edged hole in the absorber, the hole is throughout the absorber fixed in the centre

            ##comment the block ceatre cones with radius not angle
            # halfH = self.r*self.t/(4*self.r-2*self.sep)
            # halfH1 = (3 * self.r * self.t - 2 * self.sep * self.t) / (4 * self.r - 2 * self.sep)
            # bcone = vedo.Cone(alpha=0.5, r= self.r, pos=(0, 0, halfH), res=120, height=2*halfH, axis=(0, 0, 1), c='yellow') #bottom cone
            # tcone = vedo.Cone(alpha=0.5, r= self.r, pos=(0, 0, halfH1), res=120, height=2*halfH, axis=(0, 0, -1), c='yellow') #top cone
            if self.angle:
                temp = self.sep / 2 / (np.tan(self.angle))
                # temp = (np.tan(np.pi / 6)) * self.t / 2
                # self.r = temp + self.sep / 2

            if pos is not None:
                x, y, z = pos
                self.hole_relative_z = z - self.z
                z = self.z - self.t / 2
                if abs(self.hole_relative_z) <  self.t/2:
                    self.hole_side_dia_t = np.tan(self.angle) * (self.t / 2 - self.hole_relative_z + temp)  ## same result r
                    self.hole_side_dia_b = np.tan(self.angle) * (self.t / 2 + self.hole_relative_z + temp)
                    # z = z-self.t/2
                    height_t = self.hole_side_dia_t / np.tan(self.angle)
                    height_b = self.hole_side_dia_b / np.tan(self.angle)

                    bcone = vedo.Cone(alpha=0.5, r=self.hole_side_dia_b, pos=(x, y, height_b / 2 + z), res=RES, height=height_b,
                                      axis=(0, 0, 1),
                                      c=color)

                    tcone = vedo.Cone(alpha=0.5, r=self.hole_side_dia_t, pos=(x, y, self.t - height_t / 2 + z), res=RES,
                                      height=height_t, axis=(0, 0, -1),
                                      c=color)
                    return tcone, bcone

            self.hole_side_dia = np.tan(self.angle) * (self.t / 2 + temp)  ## same result r
            x, y, z = self.pos
            z = z - self.t / 2
            # x,y,z = pos
            height = self.hole_side_dia / np.tan(self.angle)

            bcone = vedo.Cone(alpha=0.5, r=self.hole_side_dia, pos=(x, y, height / 2 + z), res=RES, height=height,
                              axis=(0, 0, 1),
                              c=color)

            tcone = vedo.Cone(alpha=0.5, r=self.hole_side_dia, pos=(x, y, self.t - height / 2 + z), res=RES,
                              height=height, axis=(0, 0, -1),
                              c=color)


            return tcone, bcone
        def getCylinder(self, rad, Capping=True, pos =None, color = 'green'): #cylindrical holes in the absorber, the hole is throughout the absorber fixed in the centre or absorber
            if pos is None:
                cp1, cp2 = (self.x, self.y, self.z-self.t), (self.x, self.y, self.z+self.t)
                cylinder = vedo.Cylinder(pos=[cp1, cp2], r=rad, height=self.t, axis=(0, 0, 1), res=RES, cap=Capping,
                                         c=color)
            else:
                # x,y,z =pos
                cylinder = vedo.Cylinder(pos=pos, r = rad, height = self.t, axis=(0,0,1), res=RES, cap=Capping, c=color)

            return cylinder

        def getSphere(self, color='green'):
            sphere = vedo.Sphere(pos=self.pos, r=self.r, c=color, alpha=0.5, res=RES)
            return sphere

        def getPlane(self, color='green'):  # detector plane for imaging recorded events
            plane = vedo.Plane(pos=self.pos, normal=(0, 0, -1), sx=self.l, sy=self.w, c=color)
            return plane
        def getCircle(self,  color='green'):
            circle = vedo.Circle(pos=self.pos, r=self.sep, res=RES, alpha=0.5, c=color)
            return circle
        def setHalfAngle(self, angle):
            self.angle = angle

        def getHole(self, pos=None, color='yellow'):
            if self.hole != None:
                self.hole_pos = pos
                if self.hole == 'Kn':
                    if pos is None:
                        self.holeShape = self.getCones(color = color)
                    else:
                        # x,y,z =pos
                        # pos = x, y, z+self.t/2
                        self.holeShape = self.getCones(pos, color=color)
                elif self.hole == 'Cy':
                    self.hole_side_dia = self.sep/2
                    if pos is None:
                        self.holeShape = self.getCylinder(self.sep/2, Capping=False, color = color)
                    else:
                        self.holeShape = self.getCylinder(self.sep/2, Capping=False, pos = pos, color=color)
                else:
                    self.holeShape = None

                return self.holeShape
            else:
                self.holeShape = None
                return None

        # def setRadiusByHalfAngle(self):
        #     if self.angle:
        #         temp = self.sep / 2 / (np.tan(self.angle))
        #         # temp = (np.tan(np.pi / 6)) * self.t / 2
        #         self.r = np.tan(self.angle) * (self.t / 2 + temp)  ## same result r
        #         # self.r = temp + self.sep / 2
        #         return self.r

        # def getCircle(self, z): # for x and y coordinates
        #     if self.type == 'cylinder':
        #         self.circle = Circle(Point(0,0,z), self.sep/2)
        #     if self.type == 'knife':
        #         if z > self.t/2:
        #             l = self.r # upper bottom
        #             L = self.sep/2 #lower bottom
        #             H = self.t/2 #total height
        #             h = z - self.t/2 #height to the lower bottom
        #             radius = (l-L)*h/H +l
        #             # radius = 2*self.z*(self.sep/2)/self.t
        #             self.circle = Circle(Point(0,0,z), radius)
        #
        #         else:
        #             L = self.r
        #             l = self.sep/2
        #             H = self.t/2
        #             h = z
        #             radius = (H-h)/H*(L-l) + l
        #             self.circle = Circle(Point(0, 0, z), radius)

class Cy(Absorber):
    def __init__(self, atomNum, pos, t, r, hole_angle_half,  hole_diameter,  hole=None, property_file=None, color = 'green'):
        super().__init__(atomNum, property_file)
        self.hole = hole
        self.setPositon(pos)
        self.setThickness(t)
        self.setRadius(r)
        if hole !=None:
            self.setHalfAngle(hole_angle_half)
            self.setSepration(hole_diameter)
            # self.setHoleRePosition(hole_relative_z)
        # self.setLength(length)
        # self.setWidth(width)
        self.shape = super().getCylinder(self.r, pos = pos, color=color)
        # self.holeShape = None

    # def setRegion(self):
    #     if self.type == "block":
    #         self.r = [[i for i in range(-self.width, self.width+1)], [i for i in range(-self.height, self.height+1)], [i for i in range(0, self.t+1)]]
    #     if self.type == "pinhole":
    #         self.r = [[i for i in range(-self.width, self.width + 1)],
    #                   [i for i in range(-self.height, self.height + 1)], [i for i in range(0, self.t + 1)]]
class Bo(Absorber):
    def __init__(self, atomNum, pos, t, w, l, hole_angle_half, hole_diameter, hole=None, property_file=None, color = 'green'):
        super().__init__(atomNum, property_file)
        self.hole = hole
        self.setPositon(pos)
        self.setThickness(t)
        self.setLength(l)
        self.setWidth(w)
        if hole != None:
            self.setHalfAngle(hole_angle_half)
            self.setSepration(hole_diameter)
        self.shape = super().getBox(color)
        # self.holeShape = None
    # def getHole(self, pos=None, color='yellow'):
    #     if self.hole != None:
    #         if self.hole == 'Kn':
    #             self.holeShape =  self.getCones(pos = pos, color = color)
    #         elif self.hole == 'Cy':
    #             self.holeShape = self.getCylinder(pos = pos, color = color)
    #     return self.holeShape

class Co(Absorber):
    def __init__(self, atomNum, pos, t, r, up, property_file=None, color = 'green'):
        super().__init__(atomNum, property_file)
        self.r = r
        self.t = t
        self.setPositon(pos)
        self.shape = super().getCone(up, color)
        self.holeShape = None

class Sp(Absorber):
    def __init__(self, atomNum, pos, r, property_file=None, color = 'green'):
        super().__init__(atomNum, property_file)
        self.r = r
        self.setPositon(pos)
        self.shape = super().getSphere(color)
        self.holeShape = None
        self.t = self.r

class Pl(Absorber):
    def __init__(self, atomNum, pos, w, l, hole_angle_half, hole_diameter, hole=None, property_file=None, color = 'green'):
        super().__init__(atomNum, property_file)
        self.hole = hole
        self.setPositon(pos)
        self.w = w
        self.l = l
        self.shape = super().getPlane()
        if self.hole != None:
            self.setHalfAngle(hole_angle_half)
            self.setSepration(hole_diameter)
            self.holeShape = self.getCircle(color)

class Detector():
    def __init__(self, detector_material, density, width, length, thickness, pixel_size, pos):
        # self.x = 0
        # self.y = 0
        # self.z = 1  # unit of cm, pinhole to detector distance 10mm
        detector_dict = xraylib.CompoundParser(detector_material)
        self.detector_material = detector_material
        self.x, self.y, self.z = pos
        self.pos = pos
        self.density = density  # g/cm^3
        self.width = width  # cm, same size as CCD active area
        self.length = length  # cm, same size as CCD active area
        if detector_dict['nElements'] == 2:
            self.atom1, self.atom2 = detector_dict['Elements']
            self.atom1_weight, self.atom2_weight = detector_dict['nAtoms']
        elif detector_dict['nElements'] == 1:
            self.atom1 = detector_dict['Elements'][0]
            self.atom2 = self.atom1
            self.atom1_weight, self.atom2_weight = 1, 1
        # self.sep = 1e-5 #cm vaccum gap
        # self.diamter = 3e-3 #cm scintillator needle diameter 5-15 micro meters
        self.t = thickness  # cm 1.5mm->0.15cm, or 0.06cm
        self.mass = detector_dict['molarMass'] #xraylib.AtomicWeight(atomNumber1)*atom1_weight + xraylib.AtomicWeight(atomNumber2)*atom2_weight
        self.pixel_size = pixel_size
        self.pts = []
        # self.mass = 259.81 #g/mol
        self.scintillatorBox =  vedo.Box(pos=pos, length=self.length, width=self.width, height=self.t, c='r')
        self.energies = None
        self.csCoLs = None
        self.rsCoLs = None
        self.peCoLs = None
        self.totalCoLs = None
        # self.scintillatorPlanes = []

    def setPosition(self, pos):
        self.pos = pos
        self.x, self.y, self.z = pos
    def getScitillator(self):
        # scintillatorCylinders = dict()
        # scintillatorPlanes = []
        # xstepsize = self.sep + self.diamter
        # xstart = -self.width/2+ xstepsize
        # xend = self.width/2
        #
        # ystepsize = self.sep + self.diamter
        # ystart = -self.length / 2 + ystepsize
        # yend = self.length / 2
        # xs = np.arange(xstart, xend, xstepsize, dtype=float)
        # ys = np.arange(ystart,yend,ystepsize, dtype = float)


        ### Option 1: modelling every crystal needle as a cylinder
        # for i in range(len(xs)):
        #     for j in range(len(ys)):
        #         cp1, cp2 = (xs[i], ys[j], self.z), (xs[i],ys[j], self.z+self.t)
        #         cylinder = vedo.Cylinder(pos=[cp1, cp2], r=self.diamter/2, height=self.t, axis=(0, 0, 1), res=120, cap=True)
        #         columnPosition = (xs[i], ys[j], self.z)
        #         scintillatorCylinders[columnPosition] = cylinder
        ###Option 2: modelling the scintillator as a singe solid block with 6 planes:
        # scintillatorBox = vedo.Box(pos=(0,0,self.z+self.t/2), length=self.length, width=self.width, height=self.t, c='r')

        # scintillatorTop = vedo.Plane(pos=(0, 0, self.z+self.t), normal=(0, 0, -1), sx=self.width, sy=self.length)
        # scintillatorBottom = vedo.Plane(pos=(0, 0, self.z), normal=(0, 0, -1), sx=self.width, sy=self.length)
        #
        #
        #
        #
        # scintillatorLeft = vedo.Plane(pos=(0, -self.width/2, self.z + self.t/2), normal=(0, -1, 0), sy=self.length, sx=self.t, c='r')
        # scintillatorRight = vedo.Plane(pos=(0, self.width/2, self.z + self.t/2), normal=(0, 1, 0), sy=self.length, sx=self.t, c='r')
        # scintillatorFront = vedo.Plane(pos=(self.width/2, 0, self.z + self.t/2), normal=(1, 0, 0), sy=self.length, sx=self.t, c= 'b')
        # scintillatorBack = vedo.Plane(pos=(-self.width/2, 0, self.z + self.t/2), normal=(-1, 0, 0), sy=self.length, sx=self.t, c= 'b')

        #comment below for a faster speed
        # self.scintillatorPlanes.extend([scintillatorTop, scintillatorBottom, scintillatorLeft, scintillatorRight, scintillatorFront, scintillatorBack])
        # scintillatorMerged =  vedo.merge(self.scintillatorPlanes)#vedo.merge(scintillatorPlanes)

        # vedo.show(scintillatorTop, scintillatorBottom, scintillatorLeft, scintillatorRight, scintillatorFront, scintillatorBack)
        return self.scintillatorBox#, scintillatorMerged

    def intersectWithScintillatorPlanes(self, currentPosition, nextPosition):
        scintillatorPts = self.scintillatorBox.intersectWithLine(currentPosition, nextPosition, tol=1e-04)
        return scintillatorPts

    def getPlane(self): #detector plane for imaging recorded events
        x,y,z =self.pos
        z = z + self.t/2
        plane = vedo.Plane(pos=(x,y,z), normal=(0,0,-1), sx= self.length, sy=self.width)
        # points = vedo.Point((self.x,self.y,self.z), r=10).color('yellow')
        # vedo.show(plane, points, axes=True)
        return plane
    def getPixelArray(self):
        ccdXArray = np.arange(-self.length/2, self.length/2+self.pixel_size, self.pixel_size)
        self.ccdXRanges = []
        for i in range(len(ccdXArray)-1):
            self.ccdXRanges.append((ccdXArray[i], ccdXArray[i+1]))
        ccdYArray = np.arange(-self.width/ 2, self.width/ 2+self.pixel_size, self.pixel_size)
        self.ccdYRanges = []
        for i in range(len(ccdYArray)-1):
            self.ccdYRanges.append((ccdYArray[i], ccdYArray[i+1]))
    def checkPixel(self, px, py):
        x = None
        y= None
        for i in range(len(self.ccdXRanges)):
            r=self.ccdXRanges[i]
            if px >= r[0] and px < r[1]:
                x = i
                break
        if x is None:
            if abs(px - self.ccdXRanges[0][0]) < 1e-04:
                x = 0
            if abs(px - self.ccdXRanges[-1][-1]) < 1e-04:
                x = len(self.ccdXRanges) - 1
        for i in range(len(self.ccdYRanges)):
            r = self.ccdYRanges[i]
            if py >= r[0] and py < r[1]:
                y = i
                break
        if y is None:
            if abs(py- self.ccdYRanges[0][0]) < 1e-04:
                y = 0
            if abs(py - self.ccdYRanges[-1][-1]) < 1e-04:
                y = len(self.ccdYRanges) - 1
        return x, y
    def getCCDArray(self):
        x = int(self.length/self.pixel_size)
        y = int(self.width/self.pixel_size)
        return np.zeros((x, y))

class Detector_sub(Absorber, Detector):
        def __init__(self, atomNum, property_file=None):
            super().__init__(atomNum, property_file)
            # self.scale = [17.418674, 8.314444, 10.323193, 1.383834, 19.876251]
            # self.exponents = [0.399828, 0.016872, 25.605827, 233.339676, 3.826915]
            # self.offset = -2.322802
            # self.z = 0  # distance to source
            # self.density = 1.873  # g/cm^3
            # # self.width = 0.8192  # cm, same size as CCD active area
            # # self.length = 0.8192  # cm, same size as CCD active area
            # # self.t = 0.06 #cm
            # self.mass = 132.91  # g/mol
            #
            # ##Jump Factors
            # self.JK = 5.949
            # self.JL1 = 1.150
            # self.JL2 = 1.400
            # self.JL3 = 2.955
            #
            # ##Absorption Edges
            # self.KEDGE = 35.985  # Unit KeV K edge energy
            # self.L1EDGE = 5.714
            # self.L2EDGE = 5.359
            # self.L3EDGE = 5.012
            #
            # ##fluo yields
            # self.fyK = 0.898
            # self.fyL1 = 0.049
            # self.fyL2 = 0.090
            # self.fyL3 = 0.091
            # ## fluo transition probabilities
            # self.kA1 = 0.5240
            # self.kA2 = 0.2830
            # self.kB1 = 0.1020
            # self.kB2 = 0.0390
            # self.kB3 = 0.0520
            # self.l1G3 = 0.4220
            # self.l1B3 = 0.4720
            # self.l1B4 = 0.1060
            # self.l2B1 = 0.8310
            # self.l2G1 = 0.1690
            # self.l3L = 0.0270
            # self.l3A2 = 0.1170
            # self.l3A1 = 0.6230
            # self.l3B2 = 0.2330
            #
            # ##fluo transition energies
            # self.fluoTranEnergies = {
            # 'Kalpha1': 30.973,
            # 'Kalpha2': 30.625,
            # 'Kbeta1' : 34.987,
            # 'Kbeta3': 34.920,
            # 'Kbeta2': 35.833,
            #
            # 'L1gamma3': 5.567,
            # 'L1beta3': 4.717,
            # 'L1beta4': 4.652,
            #
            # 'L2beta1': 4.620,
            # 'L2gamma1': 5.280,
            #
            # 'L3l': 3.795,
            # 'L3alpha2': 4.272,
            # 'L3alpha1': 4.286,
            # "L3beta2": 4.936}
            #
            #
            # ##non-radiative transitions
            # ## Original Vacancy K shell
            # self.kL1L1 = 0.080
            # self.kL1L2 = 0.098
            # self.kL1L3 = 0.115
            # self.kL2L2 = 0.012
            # self.kL2L3 = 0.245
            # self.kL3L3 = 0.123
            # self.kL1X = 0.092
            # self.kL2X = 0.082
            # self.kL3X = 0.136
            # self.kXX = 0.017
            # ## L1 shell
            # self.l1L2X = 0.120
            # self.l1L3X = 0.294
            # self.l1XX = 0.586
            # ## L2 shell
            # self.l2L3X = 0.169
            # self.l2XX = 0.831
        #
        # def intersectWithScintillatorPlanes(self, currentPosition, nextPosition):
        #     return super().intersectWithScintillatorPlanes(currentPosition, nextPosition)


# class I(CsI):
#         def __init__(self):
#             # self.scale = [19.884502, 6.736593, 8.110516, 1.170953, 17.548716]
#             # self.exponents = [4.628591, 0.027754, 31.849096, 84.406387, 0.46355]
#             # self.offset = -0.448811
#
#             # self.z = 0  # distance to source
#             self.density = 4.933  # g/cm^3
#             # self.width = 0.8192  # cm, same size as CCD active area
#             # self.length = 0.8192  # cm, same size as CCD active area
#             # self.t = 0.06 #cm
#             self.mass = 126.90  # g/mol
#
#             ##Jump Factors
#             self.JK = 6.039
#             self.JL1 = 1.149
#             self.JL2 = 1.400
#             self.JL3 = 2.950
#
#             ##Absorption Edges
#             self.KEDGE = 33.169  # Unit KeV K edge energy
#             self.L1EDGE = 5.188
#             self.L2EDGE = 4.852
#             self.L3EDGE = 4.557
#
#             ##fluo yields
#             self.fyK = 0.884
#             self.fyL1 = 0.044
#             self.fyL2 = 0.079
#             self.fyL3 = 0.079
#             ## fluo transition probabilities
#             self.kA1 = 0.5260
#             self.kA2 = 0.2830
#             self.kB1 = 0.1010
#             self.kB2 = 0.0380
#             self.kB3 = 0.0520
#             self.l1G3 = 0.4140
#             self.l1B3 = 0.4830
#             self.l1B4 = 0.1030
#             self.l2B1 = 0.8650
#             self.l2G1 = 0.1350
#             self.l3L = 0.0320
#             self.l3A2 = 0.1100
#             self.l3A1 = 0.7480
#             self.l3B2 = 0.1100
#
#             ##fluo transition energies
#             self.fluoTranEnergies = {
#             'Kalpha1': 26.612,
#             'Kalpha2': 28.318,
#             'Kbeta1' : 32.295,
#             'Kbeta3': 32.240,
#             'Kbeta2': 33.054,
#
#             'L1gamma3': 5.072,
#             'L1beta3': 4.313,
#             'L1beta4': 4.257,
#
#             'L2beta1': 4.221,
#             'L2gamma1': 4.801,
#
#             'L3l': 3.485,
#             'L3alpha2': 3.926,
#             'L3alpha1': 3.938,
#             "L3beta2": 4.509}
#
#
#             ##non-radiative transitions
#             ## Original Vacancy K shell
#             self.kL1L1 = 0.078
#             self.kL1L2 = 0.094
#             self.kL1L3 = 0.115
#             self.kL2L2 = 0.012
#             self.kL2L3 = 0.250
#             self.kL3L3 = 0.126
#             self.kL1X = 0.087
#             self.kL2X = 0.079
#             self.kL3X = 0.134
#             self.kXX = 0.025
#             ## L1 shell
#             self.l1L2X = 0.188
#             self.l1L3X = 0.293
#             self.l1XX = 0.519
#             ## L2 shell
#             self.l2L3X = 0.167
#             self.l2XX = 0.833

        # def intersectWithScintillatorPlanes(self, currentPosition, nextPosition):
        #     return super().intersectWithScintillatorPlanes(currentPosition, nextPosition)


# class Al(Absorber): # entrance window
#     def __init__(self, w):
#         # self.scale = [4.730796, 2.313951, 1.54198, 1.117564, 3.154754]
#         # self.exponents = [3.628931, 43.051167, 0.09596, 108.932388, 1.555918]
#         # self.offset = 0.139509
#         self.x = 0
#         self.y = 0
#         self.z = w.z+w.t  # distance to source
#         self.density = 2.71  # cm/cm^3
#
#
#         ##Jump Factors
#         self.JK = 10.959
#         self.JL1 = 1.089
#         self.JL2 = 1
#         self.JL3 = 1
#
#         ##Absorption Edges
#         self.KEDGE = 1.56  # Unit KeV K edge energy
#         self.L1EDGE = 0.118
#         self.L2EDGE = 0.073
#         self.L3EDGE = 0.072
#
#         ##fluo yields
#         self.fyK = 0.036
#         self.fyL1 = 0
#         self.fyL2 = 0
#         self.fyL3 = 0
#         ## fluo transition probabilities
#         self.kA1 = 0.6580 #updated with xraylib, orginal was 0.6633
#         self.kA2 = 0.3290 #0.3330
#         self.kB1 = 0.0087#0.0009
#         self.kB2 = 0
#         self.kB3 = 0.0043#0.0004
#
#         self.l1G3 = 0
#         self.l1B3 = 0
#         self.l1B4 = 0
#         self.l2B1 = 0
#         self.l2G1 = 0
#         self.l3L = 0
#         self.l3A2 = 0
#         self.l3A1 = 0
#         self.l3B2 = 0
#
#         ##fluo transition energies
#         #Kbeta3 updated with xraylib calculator
#         self.fluoTranEnergies = {
#             'Kalpha1': 1.487,
#             'Kalpha2': 1.486,
#             'Kbeta1': 1.558,
#             'Kbeta3': 1.554,
#             'Kbeta2': 1.588,
#
#            }
#
#
#         ##non-radiative transitions
#         ## Original Vacancy K shell
#         self.kL1L1 = 1
#         self.kL1L2 = 0
#         self.kL1L3 = 0
#         self.kL2L2 = 0
#         self.kL2L3 = 0
#         self.kL3L3 = 0
#         self.kL1X = 0
#         self.kL2X = 0
#         self.kL3X = 0
#         self.kXX = 0
#         ## L1 shell
#         self.l1L2X = 1
#         self.l1L3X = 0
#         self.l1XX = 0
#         ## L2 shell
#         self.l2L3X = 1
#         self.l2XX = 0
# class Vaccum(Absorber):
#     def __init__(self):
#         super().__init__()

class Air(Absorber):
    def __init__(self):
        self.z = 0  # distance to source
        self.sep = 0
        self.density = 0.001225
        self.KEDGE = 3.2  # unit Kev
        self.fyK = 1
