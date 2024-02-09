import numpy as np
import random
# import xraylib
from scoop import futures

# import scoop

# import os
import sys
# from vedo import *
# import sympy
# module_path = os.path.abspath(os.getcwd())
# if module_paht not in sys.path:
#  sys.path.append(module_path)
import vedo
import Absorber_flex
import Comparison
import dataGrid
from contextlib import redirect_stdout
# from methodtools import lru_cache
# from functools import lru_cache
# from memoization import cached
# import functools
# import Coefficient
# from threading import Thread
# import matplotlib.pyplot as plt
import math
import matplotlib.pyplot as plt
# import operator
import time
from scipy.spatial import KDTree
import configparser
import ast
import xraylib
import xraydb
# import xraydb
# import multiprocessing as mp
# from ray.util.multiprocessing.pool import Pool
# from mpi4py import MPI
# from mpi4py.futures import MPIPoolExecutor
# from scipy.linalg import hadamard, subspace_angles
# from scipy.spatial.transform import Rotation as R

# from sympy import Plane, Point3D, Ray3D, Ray, Point, solve, sympify, solveset, Poly
# from sympy.geometry import intersection, Line3D
# from sympy.vector import CoordSys3D, AxisOrienter
# from sympy import Symbol, symbols
# from sympy import Eq, S
# from sympy.abc import r, x, y, z, theta, phi
import argparse

# parser = argparse.ArgumentParser(description='input camera papermeters')
# parser.add_argument('--distance', required=True, type=float,
#                         help='source to collimator distance')
# parser.add_argument('--angle',  required=True, type=float,
#                         help='source angle in degree')
# parser.add_argument('--pinhole', required=True, type=float,
#                         help='pinhole diameter')
# parser.add_argument('--secondary', required=True,
#                         help='if simulate secondary effects')
# parser.add_argument('--energy', required=True, type=float,
#                         help='photon energy')
# parser.add_argument('--ppNum', required=True, type=float,
#                         help='photon number per processor')
# parser.add_argument('--position', required=False, type=str, default='x',
#                         help='photon position')

# parser.add_argument('--para_file', required=True, type=str, help='parameter file name')
# args = parser.parse_args()
# args = parser.parse_args()

# source to collimator distance: 1cm, 5cm, 10cm
# collimator to detector distance: 10mm -> 1cm





Airenergies, AircsCoLs, AirrsCoLs, AirpeCoLs, AirtotalCoLs = dataGrid.loadFile('Airdata.txt')
config = configparser.ConfigParser()
param_f = "parameters.ini"
config.read(param_f)
source_pos = ast.literal_eval(config['Source']['source_pos'])
source_energy = float(config['Source']['source_energy'])
source_ppNum = int(float(config['Source']['source_ppNum']))
source_type = config['Source']['source_type']
source_radius = float(config['Source']['source_radius'])
source_height = float(config['Source']['source_height'])
source_length = float(config['Source']['source_length'])
source_width = float(config['Source']['source_width'])

MIN_X = ast.literal_eval(config['Simulation_Boundaries']['min_x'])
MAX_X = ast.literal_eval(config['Simulation_Boundaries']['max_x'])
MIN_Y = ast.literal_eval(config['Simulation_Boundaries']['min_y'])
MAX_Y = ast.literal_eval(config['Simulation_Boundaries']['max_y'])
VAC_Z = ast.literal_eval(config['Simulation_Boundaries']['vacuum_z_boundary'])

ffactor_file = config['Form_Factor_File']['ff_file']
if ffactor_file != 'None':
    FFs = np.genfromtxt(ffactor_file, delimiter='|', names=True)
    XS = FFs['X']

intermediate_absorber_num = int(config['Absorbers']['intermediate_absorber_num'])
def_absorbers = [None] * intermediate_absorber_num
secondary = eval(config['If_Secondary']['secondary'])
process_num = int(config['Process_Num']['process_num'])
for i in range(intermediate_absorber_num):
    # absorber_set = set()
    color = config['Absorbers']['absorber_' + str(i + 1) + '_color']
    pos = ast.literal_eval(config['Absorbers']['absorber_' + str(i + 1) + '_position'])
    atom = int(config['Absorbers']['absorber_' + str(i + 1) + '_material'])
    thickness = float(config['Absorbers']['absorber_' + str(i + 1) + '_thickness'])
    radius = float(config['Absorbers']['absorber_' + str(i + 1) + '_radius'])
    length = float(config['Absorbers']['absorber_' + str(i + 1) + '_length'])
    width = float(config['Absorbers']['absorber_' + str(i + 1) + '_width'])
    attenuation_file = config['Absorbers']['absorber_' + str(i + 1) + '_attenuation_file']
    property_file = config['Absorbers']['absorber_' + str(i + 1) + '_property_file']
    if '.py' in property_file:
        property_class = property_file.split('.py')[0]
    else:
        property_class = None
    if config['Absorbers']['absorber_' + str(i + 1) + '_shape'] != 'None':
        shape = config['Absorbers']['absorber_' + str(i + 1) + '_shape']
    else:
        shape = None
    if config['Absorbers']['absorber_' + str(i + 1) + '_hole_type'] != 'None':
        hole_type = config['Absorbers']['absorber_' + str(i + 1) + '_hole_type']
    else:
        hole_type = None

    hole_diameter = float(config['Absorbers']['absorber_' + str(i + 1) + '_hole_diameter'])
    hole_angle_half = float(config['Absorbers']['absorber_' + str(i + 1) + '_hole_angle']) * np.pi / (2 * 180)
    hole_pos= ast.literal_eval(config['Absorbers']['absorber_' + str(i + 1) + '_hole_pos'])
    if shape == 'Cy':
        absorber = Absorber_flex.Cy(atom, pos, thickness, radius, hole_angle_half, hole_diameter, hole_type,
                                    property_class, color)
        # absorber.getHole(pos=hole_pos, color='yellow')
        absorber.getHole(pos=hole_pos, color='yellow')
    elif shape == 'Bo':
        absorber = Absorber_flex.Bo(atom, pos, thickness, width, length, hole_angle_half, hole_diameter, hole_type,
                                    property_class, color)
        absorber.getHole(pos=hole_pos, color='yellow')
    elif shape == 'Co':
        up = eval(config['Absorbers']['absorber_' + str(i + 1) + '_cone_upwards'])
        absorber = Absorber_flex.Co(atom, pos, thickness, radius, up, property_class, color)
    elif shape == 'Sp':
        absorber = Absorber_flex.Sp(atom, pos, radius, property_class, color)
    elif shape == 'Pl':
        absorber = Absorber_flex.Pl(atom, pos, width, length, hole_angle_half, hole_diameter, hole_type, property_class,
                                    color)
    if ffactor_file != 'None':
        absorber.ffs = FFs[xraylib.AtomicNumberToSymbol(atom)]

    # absorber_shape = absorber.shape
    # if type(absorber_shape) is vedo.shapes.Cylinder:
    #     pos = (absorber.x, absorber.y, absorber.z - absorber.t / 2)
    #     absorber_hole = absorber.getHole(color='yellow')
    # else:
    # hole_pos = (absorber.x, absorber.y, absorber.z)
    # absorber_hole = absorber.getHole(pos=hole_pos, color='yellow')
    if attenuation_file != 'None':
        absorber.energies, absorber.csCoLs, absorber.rsCoLs, absorber.peCoLs, absorber.totalCoLs = dataGrid.loadFile(
            attenuation_file)
    # absorber_set.add(absorber)
    # if absorber_hole !=None:
    #     absorber_set.add(absorber_hole)

    # absorber.add(pos)
    # absorber.add(atom)
    # absorber.add(thickness)
    # absorber.add(shape)
    # absorber.add(hole_type)
    # absorber.add(hole_diameter)
    # absorber.add(hole_angle)
    def_absorbers[i] = absorber

# print(absorbers)


detector_pos = ast.literal_eval(config['Detector']['detector_position'])
detector_material = config['Detector']['detector_material']
# detector_atom1, detector_atom2 = detector_dict['Elements']
# atom1_weight, atom2_weight = ast.literal_eval(config['Detector']['detector_atom_weights'])

detector_thickness = float(config['Detector']['detector_thickness'])
detector_length = float(config['Detector']['detector_length'])
detector_width = float(config['Detector']['detector_width'])
detector_pixel_size = float(config['Detector']['detector_pixel_size'])
detector_density = float(config['Detector']['detector_density'])
detector_attenuation_file = config['Detector']['detector_attenuation_file']
detector_sub1_attenuation_file = config['Detector']['detector_atom1_attenuation_file']
detector_sub2_attenuation_file = config['Detector']['detector_atom2_attenuation_file']
detector_sub1_property_file = config['Detector']['detector_atom1_property_file']
detector_sub2_property_file = config['Detector']['detector_atom2_property_file']
detector = Absorber_flex.Detector(detector_material, detector_density, detector_width, detector_length,
                                  detector_thickness, detector_pixel_size, detector_pos)
if '.py' in detector_sub1_property_file:
    detector_sub1_property_class = detector_sub1_property_file.split('.py')[0]
else:
    detector_sub1_property_class = None
if '.py' in detector_sub2_property_file:
    detector_sub2_property_class = detector_sub2_property_file.split('.py')[0]
else:
    detector_sub2_property_class = None
detector_sub1 = Absorber_flex.Detector_sub(detector.atom1, detector_sub1_property_class)
detector_sub2 = Absorber_flex.Detector_sub(detector.atom2, detector_sub2_property_class)

if ffactor_file != 'None':
    detector_sub1.ffs = FFs[xraylib.AtomicNumberToSymbol(detector.atom1)]
    detector_sub2.ffs = FFs[xraylib.AtomicNumberToSymbol(detector.atom2)]
# scintllatorBlock = detector.getScitillator()

if detector_attenuation_file != 'None':
    detector.energies, detector.csCoLs, detector.rsCoLs, detector.peCoLs, detector.totalCoLs = dataGrid.loadFile(
        detector_attenuation_file)
if detector_sub1_attenuation_file != 'None':
    detector_sub1.energies, detector_sub1.csCoLs, detector_sub1.rsCoLs, detector_sub1.peCoLs, detector_sub1.totalCoLs = dataGrid.loadFile(
        detector_sub1_attenuation_file)
if detector_sub2_attenuation_file != 'None':
    detector_sub2.energies, detector_sub2.csCoLs, detector_sub2.rsCoLs, detector_sub2.peCoLs, detector_sub2.totalCoLs = dataGrid.loadFile(
        detector_sub2_attenuation_file)

air = Absorber_flex.Air()

# Alenergies, AlcsCoLs, AlrsCoLs, AlpeCoLs, AltotalCoLs= dataGrid.loadFile('Aldata.txt')
# CsIenergies, CsIcsCoLs, CsIrsCoLs, CsIpeCoLs, CsItotalCoLs = dataGrid.loadFile('CsIdata.txt')
# Wenergies, WcsCoLs, WrsCoLs, WpeCoLs, WtotalCoLs= dataGrid.loadFile('data.txt')
# Csenergies, CscsCoLs, CsrsCoLs, CspeCoLs, CstotalCoLs = dataGrid.loadFile('Csdata.txt')
# Ienergies, IcsCoLs, IrsCoLs, IpeCoLs, ItotalCoLs= dataGrid.loadFile('Idata.txt')
# ANGLES = np.arange(-180, 180.1, 0.1)/180 * np.pi


# w = def_absorbers[1]


# w = Absorber.W() # 19.3 or 19.35 g/cm^3
# # w.setSepration(0.1) #unit of cm  0.5mm, 1mm -hole diameter
# w.setSepration(args.pinhole)
# # air = Absorber.Air() #0.001225 g/cm^3
# w.setPositon(0) # cm
# # air.setPositon(0)
# #t = float("{:.3f}".format(random.uniform(0.001, 10.000)))
# w.setThickness(0.6) #unit/cm but iterate the thickness from 0.001mm-10mm #0.01
#
# w.setType('knife')  # knife cylinder
# # W.setRadius(0.2232) #
# w.setHalfAngle(np.pi/6) #half of acceptance angle 30 degrees
# w.setRadiusByHalfAngle()
# al = Absorber.Al(w)
# al.setThickness(0.1) #cm -> 1mm
# al.setRadius(w.r) #-> 20 mm diameter -> 1cm radius set as w.r for validation temporarily
# window = al.getCylinder(True).color('pink')

# wDisc = def_absorbers[1]

# wDisc = Absorber.W() # 19.3 or 19.35 g/cm^3
# wDisc.setThickness(0.6) #6mm -> 0.6cm
# wDisc.setRadius(2.25) # diameter = 45mm -> r = 2.25 cm
# wDisc.setType('cylinder')

# al = def_absorbers[0]
# window = def_absorbers[0].shape

# al = Absorber.Al(w)
# al.setThickness(0.1) #cm -> 1mm
# al.setRadius(wDisc.r) #-> 20 mm diameter -> 1cm radius set as w.r for validation temporarily
# window = al.getCylinder(True).color('pink')
# def setSimulationRegion(w, l, h):

# vaccum_z = None
# for ab in def_absorbers:
#     if vaccum_z is None:
#         vaccum_z = ab.z + ab.t/2
#     else:
#         if vaccum_z < ab.z+ab.t/2:
#             vaccum_z = ab.z + + ab.t/2


scintillator = detector
# scintillator.setPositon(1+wDisc.t/2-scintillator.t/2)
# scintillator = Absorber.CsI()
# scintillator.setPositon(1+wDisc.t/2-scintillator.t/2)
# ScintillatorBlock, ScintillatorBlockInside = scintillator.getScitillator()
ScintillatorBlock = scintillator.getScitillator()
scintillator.getPixelArray()
pixelArray = detector.getCCDArray()
CCDPlane = detector.getPlane()
# Detector = Detector(0,0,scintillator.z+scintillator.t)  #10mm -> 1cm
# Detector.getPixelArray() # actually the detector is 0.6 mm thick but currently simplified to a plane
# pixelArray = Detector.getCCDArray()
# ScintilatorEntryPlane = Scintilator.getEntryPlane()
# CCDPlane = Detector.getPlane()#unit cm
CCDPlotsX = []
CCDPlotsY = []
display = False
vis = eval(config['If_Visualisation']['visualisation'])

# discCylinder = wDisc.shape
# if np.array(w.holeShape).size == 2:
#     tcone, bcone = w.holeShape
#     # if visualisation:
#     #     vedo.show(discCylinder, tcone, bcone, window, ScintillatorBlock, CCDPlane, axes=True)
# else:
#     tcone, bcone = None, None


    # if visualisation:
    #     if type(w.holeShape) is vedo.shapes.Cylinder:
    #         vedo.show(discCylinder, vedo.Cylinder(pos=[(w.x, w.y, w.z), (w.x, w.y, w.z + w.t)], r =w.hole_side_dia , height = w.t, axis=(0,0,1), res=120, cap=True, c='yellow'),  ScintillatorBlock, CCDPlane, axes=True)
    #     else:
    #         vedo.show(discCylinder, w.holeShape, ScintillatorBlock, CCDPlane, axes=True)


# vaccum = Absorber.Vaccum()

# FFs = np.genfromtxt('FFTest1971.txt', delimiter='|', names=True)
# XS = FFs['X']
# WATOMICFFs = FFs['W']
# AlATOMICFFs = FFs['Al']
# IATOMICFFs = FFs['I']
# CsATOMICFFs = FFs['Cs']

# materialChange = False
# CurrentMaterial = Absorber.Air()

## Left-handed coordinates rotatiob matrics, clock-wise for positive rotation
def Rx(theta):
    ## right handed coordibates
    return np.matrix([[1, 0, 0],
                      [0, np.cos(theta), -np.sin(theta)],
                      [0, np.sin(theta), np.cos(theta)]])


def Ry(theta):
    ## right handed coordibates
    return np.matrix([[np.cos(theta), 0, np.sin(theta)],
                      [0, 1, 0],
                      [-np.sin(theta), 0, np.cos(theta)]])


def Rz(theta):
    ## left_handed coordinates
    # return np.matrix([[-np.sin(theta), -np.cos(theta), 0],
    #                   [np.cos(theta), -np.sin(theta), 0],
    #                   [0, 0, 1]])

    ## right handed coordibates
    return np.matrix([[np.cos(theta), -np.sin(theta), 0],
                      [np.sin(theta), np.cos(theta), 0],
                      [0, 0, 1]])


class Photon():
    def __init__(self, x, y, z, energy):
        self.x = x
        self.y = y
        self.z = z
        self.energy = energy
        self.theta, self.phi = self.samplingAngles()
        # self.theta = 0 #for the model validation
        self.s = None  # cm
        self.peCo = None
        self.csCo = None
        self.rsCo = None
        self.totalCo = None  # cm^2/g
        self.materialChange = False
        self.currentAbsorber = air
        # if self.currentAbsorber != air and self.currentAbsorber.energies == None:
        #    self.getCoefficients(self.currentAbsorber)
        # else:
        #    self.getCoefficients_file(self.currentAbsorber)
        if (self.currentAbsorber == air or self.currentAbsorber.energies != None) and energy >= 1:
            self.getCoefficients_file(self.currentAbsorber)
        else:
            self.getCoefficients(self.currentAbsorber)
        self.nextAbsorber = air
        self.Cscattered = False
        self.RScattered = False
        self.debug = False
        self.notIntersectOnEdge = False
        self.initialCheck = False
        self.collimatorCounted = False
        # self.withinScintillator = False
        # super(Photon, self).__init__()

    def checkIfOnAbsorber(self, absorber, point):
        result = False
        x = point[0]
        y = point[1]

        # z = point[2]

        if (type(absorber) is not Absorber_flex.Air) and (type(absorber) is not Absorber_flex.Detector):
            dclo = absorber.shape.closestPoint(point)
            dlen = vedo.mag(dclo - list(point))
            if absorber.hole == None:
                if (dlen < 1e-04) or absorber.shape.isInside(point):
                    result = True
            elif np.array(absorber.holeShape).size == 2:
                relative_length = np.sqrt((x - absorber.hole_pos[0]) ** 2 + (y - absorber.hole_pos[1]) ** 2)
                tcone, bcone = absorber.holeShape
                tclo = tcone.closestPoint(point)
                bclo = bcone.closestPoint(point)
                tlen = vedo.mag(tclo - list(point))
                blen = vedo.mag(bclo - list(point))
                if (tlen < 1e-04) or (blen < 1e-04) or (dlen < 1e-04) or (
                        absorber.shape.isInside(point) and (not tcone.isInside(point)) and (not bcone.isInside(point))):
                    # if abs(point[2] - (absorber.t + absorber.z)) < 1e-04 or abs(point[2] -absorber.z) < 1e-04:
                    # if abs(point[2] - (absorber.t + absorber.z)) < 1e-04 or abs(point[2] - absorber.z) < 1e-04:
                    if abs(point[2] - (absorber.t / 2 + absorber.z)) < 1e-04 or abs(
                            point[2] - (absorber.z - absorber.t / 2)) < 1e-04:


                        if ((abs(point[2] - (absorber.z + absorber.t / 2)) < 1e-04 and
                             abs(relative_length - absorber.hole_side_dia_t) < 1e-04) or
                            (abs(point[2] - (absorber.z + absorber.t / 2)) < 1e-04 and
                             (relative_length - absorber.hole_side_dia_t) >= 1e-04)) or \
                                ((abs(point[2] - (absorber.z - absorber.t / 2)) < 1e-04 and
                                  abs(relative_length - absorber.hole_side_dia_b) < 1e-04) or
                                 (abs(point[2] - (absorber.z - absorber.t / 2)) < 1e-04 and
                                  (relative_length - absorber.hole_side_dia_b) >= 1e-04)):  # interact with top circle (edge) of the top cone
                            result = True
                        else:
                            result = False
                    # elif abs(point[2] - absorber.z - absorber.t / 2) < 1e-04:
                    elif abs(point[2] - absorber.z) < 1e-04:
                        if abs(relative_length - absorber.sep / 2) < 1e-04 \
                                or (relative_length - absorber.sep / 2) >= 1e-04:  # interact with bottom circle of the top cone
                            result = True
                        else:
                            result = False
                    # elif point[2] > w.t / 2 and point[2] < w.t:  # interact with side of the top cone
                    #     result = True
                    #
                    # # if abs(point[2]) < 1e-04:
                    # #     if abs(np.sqrt(x * x + y * y) - w.r) < 1e-04:
                    # #         result = True
                    # # elif abs(point[2] - w.t / 2) < 1e-04:
                    # #     if abs(np.sqrt(x * x + y * y) - w.sep / 2) < 1e-04:
                    # #         result = True
                    # elif point[2] > 0 and point[2] < w.t / 2:
                    #     result = True
                    #
                    # if abs(point[2]) < 1e-04 or abs(point[2] - w.t) < 1e-04 or abs(
                    #         np.sqrt(x * x + y * y) - w.r) < 1e-04:
                    #     if (np.sqrt(x * x + y * y) - w.r) >= 1e-04:
                    #         result = True
                    #     # elif abs(
                    #     #         np.sqrt(x * x + y * y) - w.r) < 1e-04:
                    #     #     line = vedo.Line([pt, currentPosition])
                    #     #     distanceBetPts[pt] = line.length()
                    #     #     discPts.append(pt)

                    else:
                        result = True
            elif np.array(absorber.holeShape).size == 1:
                relative_length = np.sqrt((x - absorber.hole_pos[0]) ** 2 + (y - absorber.hole_pos[1]) ** 2)
                cy_hole = absorber.holeShape
                hole_clo = cy_hole.closestPoint(point)
                hole_len = vedo.mag(hole_clo - list(point))
                if (hole_len < 1e-04) or (dlen < 1e-04) or (
                        absorber.shape.isInside(point) and (not cy_hole.isInside(point))):
                    # if abs(point[2] - (absorber.t + absorber.z)) < 1e-04 or abs(point[2] - absorber.z) < 1e-04:
                    if abs(point[2] - (absorber.t / 2 + absorber.z)) < 1e-04 or abs(
                            point[2] - (absorber.z - absorber.t / 2)) < 1e-04:
                        if abs(relative_length - absorber.hole_side_dia) < 1e-04 or (relative_length - absorber.hole_side_dia) >= 1e-04:  # interact with top circle (edge) of the top cone
                            result = True
                        else:
                            result = False
                    # elif point[2] > w.t / 2 and point[2] < w.t:  # interact with side of the top cone
                    #     result = True
                    #
                    # # if abs(point[2]) < 1e-04:
                    # #     if abs(np.sqrt(x * x + y * y) - w.r) < 1e-04:
                    # #         result = True
                    # # elif abs(point[2] - w.t / 2) < 1e-04:
                    # #     if abs(np.sqrt(x * x + y * y) - w.sep / 2) < 1e-04:
                    # #         result = True
                    # elif point[2] > 0 and point[2] < w.t / 2:
                    #     result = True
                    #
                    # if abs(point[2]) < 1e-04 or abs(point[2] - w.t) < 1e-04 or abs(
                    #         np.sqrt(x * x + y * y) - w.r) < 1e-04:
                    #     if (np.sqrt(x * x + y * y) - w.r) >= 1e-04:
                    #         result = True
                    #     # elif abs(
                    #     #         np.sqrt(x * x + y * y) - w.r) < 1e-04:
                    #     #     line = vedo.Line([pt, currentPosition])
                    #     #     distanceBetPts[pt] = line.length()
                    #     #     discPts.append(pt)

                    else:
                        result = True

        return result

    # def checkIfOnCollimator(self, point): #including inside the collimator
    #     result = False
    #     x = point[0]
    #     y = point[1]
    #     # z = point[2]
    #
    #     tclo = tcone.closestPoint(point)
    #     bclo = bcone.closestPoint(point)
    #     dclo = discCylinder.closestPoint(point)
    #     tlen = vedo.mag(tclo-list(point))
    #     blen = vedo.mag(bclo-list(point))
    #     dlen = vedo.mag(dclo-list(point))
    #
    #
    #     if (tlen < 1e-04) or (blen < 1e-04) or (dlen < 1e-04) or (discCylinder.isInside(point) and (not tcone.isInside(point)) and (not bcone.isInside(point))):
    #         if abs(point[2] - w.t) < 1e-04 or abs(point[2]) < 1e-04:
    #             if abs(np.sqrt(x * x + y * y) - w.hole_side_dia) < 1e-04 or (np.sqrt(x * x + y * y) - w.hole_side_dia) >= 1e-04:  # interact with top circle (edge) of the top cone
    #                 result = True
    #             else:
    #                 result = False
    #         elif abs(point[2] - w.t / 2) < 1e-04:
    #             if abs(np.sqrt(x * x + y * y) - w.sep / 2) < 1e-04 or (np.sqrt(x * x + y * y) - w.sep / 2) >= 1e-04:  # interact with bottom circle of the top cone
    #                 result = True
    #             else:
    #                 result = False
    #         # elif point[2] > w.t / 2 and point[2] < w.t:  # interact with side of the top cone
    #         #     result = True
    #         #
    #         # # if abs(point[2]) < 1e-04:
    #         # #     if abs(np.sqrt(x * x + y * y) - w.r) < 1e-04:
    #         # #         result = True
    #         # # elif abs(point[2] - w.t / 2) < 1e-04:
    #         # #     if abs(np.sqrt(x * x + y * y) - w.sep / 2) < 1e-04:
    #         # #         result = True
    #         # elif point[2] > 0 and point[2] < w.t / 2:
    #         #     result = True
    #         #
    #         # if abs(point[2]) < 1e-04 or abs(point[2] - w.t) < 1e-04 or abs(
    #         #         np.sqrt(x * x + y * y) - w.r) < 1e-04:
    #         #     if (np.sqrt(x * x + y * y) - w.r) >= 1e-04:
    #         #         result = True
    #         #     # elif abs(
    #         #     #         np.sqrt(x * x + y * y) - w.r) < 1e-04:
    #         #     #     line = vedo.Line([pt, currentPosition])
    #         #     #     distanceBetPts[pt] = line.length()
    #         #     #     discPts.append(pt)
    #
    #
    #         else:
    #             result = True
    #     # else:
    #     #     return False
    #     return result

    # def checkIfOnWindow(self, point): #including inside the collimator
    #
    #     sclo = window.closestPoint(point)
    #     slen = vedo.mag(sclo-list(point))
    #
    #
    #     if (slen < 1e-04) or window.isInside(point):
    #         return True
    #     else:
    #         return False

    def checkIfOnScintillator(self, point):  # including inside the collimator
        x, y, z = point
        if (z >= scintillator.z - scintillator.t / 2) and (z <= scintillator.z + scintillator.t / 2) and (
                y >= scintillator.y - scintillator.width / 2) and (y <= scintillator.width / 2 + scintillator.y) and (
                x >= scintillator.x - scintillator.length / 2) and (x <= scintillator.length / 2 + scintillator.x):
            inside = True
        else:
            inside = False
        if inside:
            return True
        else:

            # filtered_list = [plane for plane in ScintillatorBlock if
            #                  (vedo.mag(plane.closestPoint(point)-list(point))< 1e-04)]
            # if len(filtered_list) !=0:
            #     return True
            # return False

            # for plane in ScintillatorBlock:

            sclo = ScintillatorBlock.closestPoint(point)
            slen = vedo.mag(sclo - list(point))

            if (slen < 1e-04):  # or plane.isInside(point):
                # if testRe!=True:
                #     print(testRe, 'error')
                return True
            # print(testRe)
            return False

    # def checkIfOnDetector(self, point, detector):
    #     dclo = detector.closestPoint(point)
    #     dlen = vedo.mag(dclo - list(point))
    #
    #     if dlen < 1e-04:
    #         return True
    def checkIfInSimulationRegion(self, point):
        x, y, z = point
        source_x, source_y, source_z = source_pos
        max_x = MAX_X
        min_x = MIN_X
        max_y = MAX_Y
        min_y = MIN_Y
        for ab in def_absorbers:
            if ab.t > 0.0001:
                if type(ab.shape) is not vedo.shapes.Box and type(ab.shape) is not vedo.shapes.Plane:
                    if ab.r is not None:
                        if min_x is None:
                            min_x = ab.x - ab.r
                        else:
                            if ab.x - ab.r < min_x:
                                min_x = ab.x - ab.r
                        if max_x is None:
                            max_x = ab.r + ab.x
                        else:
                            if ab.r + ab.x > max_x:
                                max_x = ab.r + ab.x
                        if min_y is None:
                            min_y = ab.y - ab.r
                        else:
                            if ab.y - ab.r < min_y:
                                min_y = ab.y - ab.r
                        if max_y is None:
                            max_y = ab.r + ab.y
                        else:
                            if ab.r + ab.y > max_y:
                                max_y = ab.r + ab.y
                else:
                    if ab.l is not None:
                        if min_x is None:
                            min_x = ab.x - ab.l / 2
                        else:
                            if ab.x - ab.l / 2 < min_x:
                                max_x = ab.l / 2 + ab.x
                        if max_x is None:
                            max_x = ab.l / 2 + ab.x
                        else:
                            if ab.l / 2 + ab.x > max_x:
                                max_x = ab.l / 2 + ab.x
                    if ab.w is not None:
                        if min_y is None:
                            min_y = ab.y - ab.w / 2
                        else:
                            if ab.y - ab.w / 2 < min_y:
                                min_y = ab.y - ab.w / 2
                        if max_y is None:
                            max_y = ab.w / 2 + ab.y
                        else:
                            if ab.w / 2 + ab.y > max_y:
                                max_y = ab.w / 2 + ab.y
        if max_y is None or max_y == 0 or max_y < detector.width / 2 + detector.y:
            max_y = detector.width / 2 + detector.y
        if min_y is None or min_y == 0 or min_y > detector.y - detector.width / 2:
            min_y = detector.y - detector.width / 2
        if max_x is None or max_x == 0 or max_x < detector.length / 2 + detector.x:
            max_x = detector.length / 2 + detector.x
        if min_x is None or min_x == 0 or min_x > detector.x - detector.length / 2:
            min_x = detector.x - detector.length / 2

        if x > max_x or x < min_x or y > max_y or y < min_y or z < source_z or z > (
                detector.z + detector.t / 2):  # 0.06 is the thickness of the detector
            return False
        else:
            return True

    def checkCurrentAbsorber(self, currentPosition):
        # currentPosition = (self.x, self.y, self.z)
        absobersDict = dict()
        if self.materialChange == False and self.notIntersectOnEdge:
            for i in range(len(def_absorbers)):
                absorber = def_absorbers[i]
                if self.checkIfOnAbsorber(absorber, currentPosition):
                    source_cloPt = absorber.shape.closestPoint(source_pos)
                    source_len = vedo.mag(source_cloPt - list(source_pos))
                    absobersDict[absorber] = source_len

            if len(absobersDict) != 0:
                currentAbsorber = sorted(absobersDict.items(), key=lambda x: x[1], reverse=True)[0][
                    0]  # get the windows when the point is on the collimator and window joint surface
                self.currentAbsorber = currentAbsorber
                return
            #
            # if self.checkIfOnCollimator(currentPosition) and not self.checkIfOnWindow(currentPosition):
            #     self.currentAbsorber = w
            elif self.checkIfOnScintillator(currentPosition):
                self.currentAbsorber = scintillator
            # elif currentPosition[2] <= (al.z+al.t):
            #     self.currentAbsorber = air
            else:
                self.currentAbsorber = air

    def checkInitialAbsorber(self, currentPosition):
        # currentPosition = (self.x, self.y, self.z)
        absobersDict = dict()
        for i in range(len(def_absorbers)):
            absorber = def_absorbers[i]
            if self.checkIfOnAbsorber(absorber, currentPosition):
                source_cloPt = absorber.shape.closestPoint(source_pos)
                source_len = vedo.mag(source_cloPt - list(source_pos))
                absobersDict[absorber] = source_len

        if len(absobersDict) != 0:
            currentAbsorber = sorted(absobersDict.items(), key=lambda x: x[1], reverse=True)[0][
                0]  # get the windows when the point is on the collimator and window joint surface
            self.currentAbsorber = currentAbsorber
            return
        #
        # if self.checkIfOnCollimator(currentPosition) and not self.checkIfOnWindow(currentPosition):
        #     self.currentAbsorber = w
        elif self.checkIfOnScintillator(currentPosition):
            self.currentAbsorber = scintillator
        # elif currentPosition[2] <= (al.z+al.t):
        #     self.currentAbsorber = air
        else:
            self.currentAbsorber = air

    def getInWhichAbsorber(self, currentPosition):
        absobersDict = dict()
        # if self.materialChange == False and self.notIntersectOnEdge:
        for i in range(len(def_absorbers)):
            absorber = def_absorbers[i]
            if self.checkIfOnAbsorber(absorber, currentPosition):
                source_cloPt = absorber.shape.closestPoint(source_pos)
                source_len = vedo.mag(source_cloPt - list(source_pos))
                absobersDict[absorber] = source_len
        if len(absobersDict) != 0:
            currentAbsorber = sorted(absobersDict.items(), key=lambda x: x[1], reverse=True)[
                0][0]  # get th
            return currentAbsorber
        #
        # if self.checkIfOnCollimator(currentPosition) and not self.checkIfOnWindow(currentPosition):
        #     self.currentAbsorber = w
        elif self.checkIfOnScintillator(currentPosition):
            return scintillator
        # elif currentPosition[2] <= (al.z+al.t):
        #     self.currentAbsorber = air
        else:
            return air

    # def checkCurrentAbsorber(self, currentPosition):
    #     # currentPosition = (self.x, self.y, self.z)
    #     if self.materialChange == False and self.notIntersectOnEdge:
    #         if self.checkIfOnCollimator(currentPosition) and not self.checkIfOnWindow(currentPosition):
    #             self.currentAbsorber = w
    #         elif self.checkIfOnScintillator(currentPosition):
    #             self.currentAbsorber = scintillator
    #         elif self.checkIfOnWindow(currentPosition):
    #             self.currentAbsorber = al
    #         # elif currentPosition[2] <= (al.z+al.t):
    #         #     self.currentAbsorber = air
    #         else:
    #             self.currentAbsorber = air

    # def checkInitialAbsorber(self, currentPosition): # same as checking absorber
    #     # currentPosition = (self.x, self.y, self.z)
    #     # if self.materialChange == False and self.notIntersectOnEdge:
    #     if self.checkIfOnCollimator(currentPosition) and not self.checkIfOnWindow(currentPosition):
    #         self.currentAbsorber = w
    #     elif self.checkIfOnScintillator(currentPosition):
    #         self.currentAbsorber = scintillator
    #     elif self.checkIfOnWindow(currentPosition):
    #         self.currentAbsorber = al
    #         # elif currentPosition[2] <= (al.z+al.t):
    #         #     self.currentAbsorber = air
    #     else:
    #         self.currentAbsorber = air

    def checkNextAbsorber(self, pt, intersectionPointsWithScintillator):
        # currentPosition = (self.x, self.y, self.z)
        intersectionPointsWithScintillatorTuples = [tuple(x) for x in intersectionPointsWithScintillator]

        # intersectionPointsWithWindow = [tuple(x) for x in intersectionPointsWithWindow]
        self.notIntersectOnEdge = False
        absobersDict = dict()
        # if self.checkIfOnCollimator(pt) and not self.checkIfOnWindow(pt):
        for i in range(len(def_absorbers)):
            absorber = def_absorbers[i]
            if type(absorber) is not Absorber_flex.Air and type(absorber) is not Absorber_flex.Detector:
                # if absorber.hole != None:
                #     if len(absorber.holeShape == 2):
                topPts = [tuple(x) for x in absorber.topPts]
                botPts = [tuple(x) for x in absorber.botPts]
                discPts = [tuple(x) for x in absorber.pts]
                if self.checkIfOnAbsorber(absorber, pt):
                    # and ((len(topPts)+len(botPts)+ len(discPts))>1):
                    if pt in discPts and len(discPts) > 1:
                        self.notIntersectOnEdge = True
                    if pt in topPts and len(topPts) > 1:
                        self.notIntersectOnEdge = True
                    if pt in botPts and len(botPts) > 1:
                        self.notIntersectOnEdge = True
                    source_cloPt = absorber.shape.closestPoint(source_pos)
                    source_len = vedo.mag(source_cloPt - list(source_pos))
                    absobersDict[absorber] = source_len

        if len(absobersDict) != 0:
            nextAbsorber = sorted(absobersDict.items(), key=lambda x: x[1], reverse=True)[
                0][0]  # get the windows when the point is on the collimator and window joint surface
            # self.currentAbsorber = currentAbsorber
            self.nextAbsorber = nextAbsorber
            return



        elif pt in intersectionPointsWithScintillatorTuples:
            if len(intersectionPointsWithScintillatorTuples) > 1:
                self.notIntersectOnEdge = True
            self.nextAbsorber = scintillator

            # elif self.checkIfOnWindow(pt):
            #     if len(intersectionPointsWithWindow)>1:
            #         self.notIntersectOnEdge = True
            #     self.nextAbsorber = al
            # elif pt[2] <= (al.t + al.z):
            #     self.nextAbsorber = air
        else:
            self.nextAbsorber = air

    # def checkNextAbsorber(self, pt,
    #                       topPts, botPts, intersectionPointsWithScintillator,
    #                       discPts, intersectionPointsWithWindow):
    #     # currentPosition = (self.x, self.y, self.z)
    #     intersectionPointsWithScintillatorTuples = [tuple(x) for x in intersectionPointsWithScintillator]
    #     topPts = [tuple(x) for x in topPts]
    #     botPts = [tuple(x) for x in botPts]
    #     discPts = [tuple(x) for x in discPts]
    #     intersectionPointsWithWindow = [tuple(x) for x in intersectionPointsWithWindow]
    #     self.notIntersectOnEdge = False
    #     if self.checkIfOnCollimator(pt) and not self.checkIfOnWindow(pt):
    #         # and ((len(topPts)+len(botPts)+ len(discPts))>1):
    #         if pt in discPts and len(discPts)>1:
    #             self.notIntersectOnEdge = True
    #         if pt in topPts and len(topPts)>1:
    #             self.notIntersectOnEdge = True
    #         if pt in botPts and len(botPts)>1:
    #             self.notIntersectOnEdge = True
    #
    #         self.nextAbsorber = w
    #
    #     elif pt in intersectionPointsWithScintillatorTuples:
    #         if len(intersectionPointsWithScintillatorTuples)>1:
    #             self.notIntersectOnEdge = True
    #         self.nextAbsorber = scintillator
    #
    #     elif self.checkIfOnWindow(pt):
    #         if len(intersectionPointsWithWindow)>1:
    #             self.notIntersectOnEdge = True
    #         self.nextAbsorber = al
    #     # elif pt[2] <= (al.t + al.z):
    #     #     self.nextAbsorber = air
    #     else:
    #         self.nextAbsorber = air

    # def listToTuple(function):
    #     def wrapper(*args):
    #         args = [tuple(x) if type(x) == list else x for x in args]
    #         result = function(*args)
    #         result = tuple(result) if type(result) == list else result
    #         return result
    #
    #     return wrapper
    #
    # # your cached function
    # @listToTuple
    # @cached
    # def collimatorCountCheck(self, windowPoints, window, pt, collimatorCount): #validate collimator sensitivity by the window
    #     if len(windowPoints) != 0:
    #         if self.checkIfOnAbsorber(window, pt) and abs( pt[2] - (al.z - al.t / 2)) < 1e-04 and self.collimatorCounted == False:
    #             # if self.checkIfOnAbsorber(al, pt) and abs(pt[2] - al.z) <1e-04 and self.collimatorCounted == False:
    #             # if self.checkIfOnWindow(pt) and abs(pt[2] - al.z) <1e-04 and self.collimatorCounted == False:
    #             # if currentPosition != (0,0,distanceSA):
    #             #     print('additional collimator count')
    #             collimatorCount += 1
    #             self.collimatorCounted = True
    #             return collimatorCount

    def collimatorCountCheck(self, collimator, pt, collimatorCount): #validate collimator sensitivity by the window
        # if len(windowPoints) != 0:
            if abs(pt[2] - (collimator.z + collimator.t / 2)) <1e-04 and self.collimatorCounted == False:
                collimatorCount += 1
                self.collimatorCounted = True
                return collimatorCount
            else:
                return collimatorCount

    def main(self, traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
             eventSequence, depenergies):  # density = air.density
        # global pixelArray
        # global collimatorCount

        if self.initialCheck == False:
            self.checkInitialAbsorber((self.x, self.y, self.z))
            self.initialCheck = True

        # line = vedo.Line([(0,0,-0.5,),(0.09661432974832347, 0.056447345966953805, 1.1893493800941657)]).color('p')
        # # sctpt = vedo.Points([currentPosition, nextPosition], r=10).color('red')
        # sctpt1 = vedo.Points([(0,0,-0.5,),(0.09661432974832347, 0.056447345966953805, 1.1893493800941657)], r=10).color('black')
        # line = vedo.Line([(0, 0, 0,), (0,0,1)]).color('p')
        # point1= tcone.intersectWithLine((0, 0, 0,), (0,0,1))
        # point2 = bcone.intersectWithLine((0, 0, 0,), (0,0,1))
        # sctpt2 = vedo.Points(point2,
        #                      r=10).color('black')
        # sctpt1 = vedo.Points(point1,
        #                     r=10).color('black')
        # vedo.show(tcone, bcone, line, sctpt1, sctpt2,axes=True)
        while True:
            currentPosition = (self.x, self.y, self.z)

            self.checkCurrentAbsorber(currentPosition)
            # if self.materialChange == False:
            #     if self.checkIfOnCollimator(currentPosition):
            #         self.currentAbsorber = w
            #     elif self.checkIfOnScintillator(currentPosition):
            #         self.currentAbsorber = scintillator
            #
            #     else:
            #         self.currentAbsorber = air

            if len(pathPts) > 1:
                if not np.array_equal(currentPosition, pathPts[-1]):
                    traLengths.append(self.s)
                    pathPts.append(currentPosition)
                    absorbers.append(self.currentAbsorber)

                # else:
                #     debug = True

            else:
                traLengths.append(self.s)
                pathPts.append(currentPosition)
                absorbers.append(self.currentAbsorber)

            self.s = self.getNextPostion(self.currentAbsorber)
            # self.s = self.getNextPostionDebug(self.currentAbsorber, self.debug)

            nextPosition = (self.x, self.y, self.z)

            distanceBetPts = dict()
            for absorber in def_absorbers:
                    tcone, bcone = None, None
                # if absorber != air and absorber != detector:
                    if absorber.holeShape != None:
                        if np.array(absorber.holeShape).size == 2:
                            tcone, bcone = absorber.holeShape
                        # if visualisation:
                        #     vedo.show(discCylinder, tcone, bcone, window, ScintillatorBlock, CCDPlane, axes=True)
                        else:
                            tcone, bcone = None, None
                    if absorber.t > 0.0001:

                        if tcone is not None and bcone is not None:
                            absorber.topPts = tcone.intersectWithLine(currentPosition, nextPosition, tol=1e-04)
                            absorber.botPts = bcone.intersectWithLine(currentPosition, nextPosition, tol=1e-04)

                        absorber.pts = absorber.shape.intersectWithLine(currentPosition, nextPosition, tol=1e-04)
                    if len(absorber.topPts) != 0:
                        absorber.topPts = [pt for pt in absorber.topPts if vedo.mag(pt - currentPosition) >= 1e-04]
                    if len(absorber.botPts) != 0:
                        absorber.botPts = [pt for pt in absorber.botPts if vedo.mag(pt - currentPosition) >= 1e-04]
                    if len(absorber.pts) != 0:
                        absorber.pts = [pt for pt in absorber.pts if vedo.mag(pt - currentPosition) >= 1e-04]

                    if tcone is not None and bcone is not None:
                        topPts = []
                        for pt in absorber.topPts:
                            pt = tuple(e for e in pt)
                            x = pt[0]
                            y = pt[1]
                            relative_length = np.sqrt((x - absorber.hole_pos[0]) ** 2 + (y - absorber.hole_pos[1]) ** 2)
                            # if abs(pt[2] - w.t) < 1e-04:
                            if abs(pt[2] - (absorber.t / 2 + absorber.z)) < 1e-04:
                                if abs(relative_length - absorber.hole_side_dia_t) < 1e-04:  # interact with top circle (edge) of the top cone
                                    line = vedo.Line([pt, currentPosition])
                                    distanceBetPts[pt] = line.length()
                                    topPts.append(pt)
                            # elif abs(pt[2] - w.t/2) < 1e-04:
                            # elif abs(pt[2] - absorber.z) < 1e-04:
                            elif abs(pt[2] - absorber.hole_pos[2]) < 1e-04:
                                if abs(relative_length - absorber.sep / 2) < 1e-04:  # interact with bottom circle of the top cone
                                    line = vedo.Line([pt, currentPosition])
                                    distanceBetPts[pt] = line.length()
                                    topPts.append(pt)
                            elif pt[2] > absorber.z and pt[2] < absorber.z + absorber.t / 2:
                                # elif pt[2] > w.t/2 and pt[2] < w.t:
                                line = vedo.Line([pt, currentPosition])  # interact with side of the top cone
                                distanceBetPts[pt] = line.length()
                                topPts.append(pt)
                        absorber.topPts = topPts
                        botPts = []
                        for pt in absorber.botPts:
                            pt = tuple(e for e in pt)
                            x = pt[0]
                            y = pt[1]
                            # if abs(pt[2]) < 1e-04:
                            relative_length = np.sqrt((x - absorber.hole_pos[0]) ** 2 + (y - absorber.hole_pos[1]) ** 2)
                            if abs(pt[2] - (absorber.z - absorber.t / 2)) < 1e-04:
                                if abs(relative_length - absorber.hole_side_dia_b) < 1e-04:
                                    line = vedo.Line([pt, currentPosition])
                                    distanceBetPts[pt] = line.length()
                                    botPts.append(pt)
                            # elif abs(pt[2] - w.t / 2) < 1e-04:
                            # elif abs(pt[2] - absorber.z) < 1e-04:
                            elif abs(pt[2] - absorber.hole_pos[2]) < 1e-04:
                                if abs(relative_length - absorber.sep / 2) < 1e-04:
                                    line = vedo.Line([pt, currentPosition])
                                    distanceBetPts[pt] = line.length()
                                    botPts.append(pt)
                            # elif pt[2]>0 and pt[2]<w.t/2:
                            elif pt[2] > absorber.z - absorber.t / 2 and pt[2] < absorber.z:
                                line = vedo.Line([pt, currentPosition])
                                distanceBetPts[pt] = line.length()
                                botPts.append(pt)
                        absorber.botPts = botPts

                    pts = []
                    # testDebugMode =False
                    for pt in absorber.pts:

                        pt = tuple(e for e in pt)
                        x = pt[0]
                        y = pt[1]


                        # if abs(pt[2]) < 1e-04 or abs(pt[2] - w.t) < 1e-04:
                        if (abs(pt[2] - (absorber.z - absorber.t / 2)) < 1e-04 or abs(
                                pt[2] - (absorber.z + absorber.t / 2)) < 1e-04) and absorber.hole is not None and \
                                (type(absorber.shape) == vedo.shapes.Box or type(absorber.shape) == vedo.shapes.Cylinder):
                            relative_length = np.sqrt((x - absorber.hole_pos[0]) ** 2 + (y - absorber.hole_pos[1]) ** 2)
                            # if (abs(pt[2] - (absorber.z - absorber.t / 2)) < 1e-04 and (np.sqrt(x * x + y * y) - absorber.hole_side_dia_b) >= 1e-04) or \
                            #         ((abs(pt[2] - (absorber.z + absorber.t / 2)) < 1e-04 and (np.sqrt(x * x + y * y) - absorber.hole_side_dia_t) >= 1e-04)):
                            if (abs(pt[2] - (absorber.z - absorber.t / 2)) < 1e-04 and (
                                    relative_length - absorber.hole_side_dia_b) >= 1e-04) or \
                                    ((abs(pt[2] - (absorber.z + absorber.t / 2)) < 1e-04 and
                                      (relative_length - absorber.hole_side_dia_t) >= 1e-04)):
                                line = vedo.Line([pt, currentPosition])
                                distanceBetPts[pt] = line.length()
                                pts.append(pt)
                            elif len(absorber.botPts) == 0 and len(absorber.topPts) == 0 and \
                                    ((abs(relative_length - absorber.hole_side_dia_b) < 1e-04
                                      and abs(pt[2] - (absorber.z - absorber.t / 2)) < 1e-04) or
                                        (abs(pt[2] - (absorber.z + absorber.t / 2)) < 1e-04
                                        and (relative_length - absorber.hole_side_dia_t) >= 1e-04)):
                                line = vedo.Line([pt, currentPosition])
                                distanceBetPts[pt] = line.length()
                                pts.append(pt)
                        else:
                            line = vedo.Line([pt, currentPosition])
                            distanceBetPts[pt] = line.length()
                            pts.append(pt)

                        # absorber.topPts = np.array(topPts)
                        # absorber.botPts = np.array(botPts)
                    absorber.pts = np.array(pts)









            # intersectionPointsWithTop, w.tpts = [], []
            # intersectionPointsWithBottom, w.bpts = [], []
            #
            # intersectionPointsWithDisc, w.pts = [], []
            # intersectionPointsWithWindow, al.pts = [], []
            # intersectionPointsWithWindow = []

            intersectionPointsWithDetector = []
            intersectionPointsWithScintillator = []
            # if w.t > 0.0001:
            #     if tcone is not None and bcone is not None:
            #         w.tpts = tcone.intersectWithLine(currentPosition, nextPosition, tol=1e-04)
            #
            #         w.bpts = bcone.intersectWithLine(currentPosition, nextPosition, tol=1e-04)
            #
            #     w.pts = discCylinder.intersectWithLine(currentPosition, nextPosition, tol=1e-04)

            # if al.t > 0.0001:
            #     al.pts = window.intersectWithLine(currentPosition, nextPosition, tol=1e-04)

            scintillatorPts = scintillator.intersectWithScintillatorPlanes(currentPosition, nextPosition)

            # scintillatorPts = ScintillatorBlock.intersectWithLine(currentPosition, nextPosition, tol=1e-04)

            detectorPts = CCDPlane.intersectWithLine(currentPosition, nextPosition, tol=1e-04)

            # vedo.show(discCylinder, window, ScintillatorBlock, CCDPlane, tcone, bcone, axes=True)
            # if len(w.tpts) != 0:
            #     intersectionPointsWithTop += [pt for pt in w.tpts if vedo.mag(pt - currentPosition) >= 1e-04]
            #     # tpts = vedo.Points(intersectionPointsWithTop, r=10).color('red')
            #     # vedo.show(discCylinder, tcone, bcone,tpts, axes=True)

            # if len(w.bpts) != 0:
            #     intersectionPointsWithBottom += [pt for pt in w.bpts if vedo.mag(pt - currentPosition) >= 1e-04]
            #     # bpts = vedo.Points(intersectionPointsWithBottom, r=10).color('red')
            #     # vedo.show(discCylinder, tcone, bcone,bpts, axes=True)

            # if len(w.pts) != 0:
            #     intersectionPointsWithDisc += [pt for pt in w.pts if vedo.mag(pt - currentPosition) >= 1e-04]
            #     # discPts = vedo.Points(intersectionPointsWithDisc, r=10).color('red')
            #     # vedo.show(discCylinder, tcone, bcone,discPts, axes=True)
            # if len(al.pts) != 0:
            #     intersectionPointsWithWindow += [pt for pt in al.pts if vedo.mag(pt - currentPosition) >= 1e-04]

            if len(scintillatorPts) != 0:
                intersectionPointsWithScintillator += [pt for pt in scintillatorPts if
                                                       vedo.mag(pt - currentPosition) >= 1e-04]

            if len(detectorPts) != 0:
                intersectionPointsWithDetector += [pt for pt in detectorPts if vedo.mag(pt - currentPosition) >= 1e-04]
            # for pt in intersectionPointsWithDetector:
            #     if vedo.mag(pt - currentPosition) >= 1e-04 and vedo.mag(pt - currentPosition)< 1e-03:
            #         print(self.s)

            # distanceBetPts = dict()
            # w.topPts = []
            # for pt in intersectionPointsWithTop:
            #     pt = tuple(e for e in pt)
            #     x = pt[0]
            #     y = pt[1]
            #     # if abs(pt[2] - w.t) < 1e-04:
            #     if abs(pt[2] - (w.t / 2 + w.z)) < 1e-04:
            #         if abs(np.sqrt(
            #                 x * x + y * y) - w.hole_side_dia) < 1e-04:  # interact with top circle (edge) of the top cone
            #             line = vedo.Line([pt, currentPosition])
            #             distanceBetPts[pt] = line.length()
            #             w.topPts.append(pt)
            #     # elif abs(pt[2] - w.t/2) < 1e-04:
            #     elif abs(pt[2] - w.z) < 1e-04:
            #         if abs(np.sqrt(x * x + y * y) - w.sep / 2) < 1e-04:  # interact with bottom circle of the top cone
            #             line = vedo.Line([pt, currentPosition])
            #             distanceBetPts[pt] = line.length()
            #             w.topPts.append(pt)
            #     elif pt[2] > w.z and pt[2] < w.z + w.t / 2:
            #         # elif pt[2] > w.t/2 and pt[2] < w.t:
            #         line = vedo.Line([pt, currentPosition])  # interact with side of the top cone
            #         distanceBetPts[pt] = line.length()
            #         w.topPts.append(pt)
            # w.botPts = []
            # for pt in intersectionPointsWithBottom:
            #     pt = tuple(e for e in pt)
            #     x = pt[0]
            #     y = pt[1]
            #     # if abs(pt[2]) < 1e-04:
            #     if abs(pt[2] - (w.z - w.t / 2)) < 1e-04:
            #         if abs(np.sqrt(x * x + y * y) - w.hole_side_dia) < 1e-04:
            #             line = vedo.Line([pt, currentPosition])
            #             distanceBetPts[pt] = line.length()
            #             w.botPts.append(pt)
            #     # elif abs(pt[2] - w.t / 2) < 1e-04:
            #     elif abs(pt[2] - w.z) < 1e-04:
            #         if abs(np.sqrt(x * x + y * y) - w.sep / 2) < 1e-04:
            #             line = vedo.Line([pt, currentPosition])
            #             distanceBetPts[pt] = line.length()
            #             w.botPts.append(pt)
            #     # elif pt[2]>0 and pt[2]<w.t/2:
            #     elif pt[2] > w.z - w.t / 2 and pt[2] < w.z:
            #         line = vedo.Line([pt, currentPosition])
            #         distanceBetPts[pt] = line.length()
            #         w.botPts.append(pt)
            #
            # w.pts = []
            # # testDebugMode =False
            # for pt in intersectionPointsWithDisc:
            #
            #     pt = tuple(e for e in pt)
            #     x = pt[0]
            #     y = pt[1]
            #     # if abs(pt[2]) < 1e-04 or abs(pt[2] - w.t) < 1e-04:
            #     if (abs(pt[2] - (w.z - w.t / 2)) < 1e-04 or abs(
            #             pt[2] - (w.z + w.t / 2)) < 1e-04) and w.hole is not None and \
            #             (type(w.shape) == vedo.shapes.Box or type(w.shape) == vedo.shapes.Cylinder):
            #         if (np.sqrt(x * x + y * y) - w.hole_side_dia) >= 1e-04:
            #             line = vedo.Line([pt, currentPosition])
            #             distanceBetPts[pt] = line.length()
            #             w.pts.append(pt)
            #         elif len(intersectionPointsWithBottom) == 0 and len(intersectionPointsWithTop) == 0 and abs(
            #                 np.sqrt(x * x + y * y) - w.hole_side_dia) < 1e-04:
            #             line = vedo.Line([pt, currentPosition])
            #             distanceBetPts[pt] = line.length()
            #             w.pts.append(pt)
            #     else:
            #         line = vedo.Line([pt, currentPosition])
            #         distanceBetPts[pt] = line.length()
            #         w.pts.append(pt)
            #
            # for pt in intersectionPointsWithWindow:
            #     pt = tuple(e for e in pt)
            #     line = vedo.Line([pt, currentPosition]).color('r')
            #     distanceBetPts[pt] = line.length()
                # sctpt = vedo.Points(intersectionPointsWithScintillator, r=10).color('red')

            for pt in intersectionPointsWithScintillator:
                pt = tuple(e for e in pt)
                line = vedo.Line([pt, currentPosition]).color('p')
                distanceBetPts[pt] = line.length()
                # sctpt = vedo.Points(intersectionPointsWithScintillator, r=10).color('red')
                # ScintillatorBlock.color('b')
                # vedo.show(bcone, tcone, line, sctpt, ScintillatorBlock, axes=
                #           True)

            for pt in intersectionPointsWithDetector:
                pt = tuple(e for e in pt)
                line = vedo.Line([pt, currentPosition])
                distanceBetPts[pt] = line.length()

            if len(distanceBetPts.keys()) != 0:

                pt = min(distanceBetPts, key=distanceBetPts.get)
                # if abs(pt[2] - al.z) < 1e-04 and len(pathPts) >3:
                #     print('collimate', pathPts, al.z, pt[2])
                if (self.currentAbsorber != air) and (
                        not self.checkIfOnAbsorber(self.currentAbsorber, pt)) and self.currentAbsorber != scintillator:

                    # if (self.currentAbsorber == w or self.currentAbsorber==wDisc) and not self.checkIfOnCollimator(pt):
                    self.materialChange = True

                # if self.checkIfOnScintillator(pt):
                #     print(pt)

                # if self.checkIfOnWindow(pt):
                #     self.nextAbsorber = al
                # elif self.checkIfOnScintillator(pt):
                #     self.nextAbsorber = scintillator
                # else:
                #     self.nextAbsorber = air
                # elif (self.currentAbsorber == al) and not self.checkIfOnWindow(pt):
                #     self.materialChange = True
                # if self.checkIfOnCollimator(pt):
                #     self.nextAbsorber = w
                # elif self.checkIfOnScintillator(pt):
                #     self.nextAbsorber = scintillator
                # else:
                #     self.nextAbsorber = air
                elif (self.currentAbsorber == scintillator) and not self.checkIfOnScintillator(pt):
                    self.materialChange = True
                    # if self.checkIfOnCollimator(pt):
                    #     self.nextAbsorber = w
                    # elif self.checkIfOnWindow(pt):
                    #     self.nextAbsorber = al
                    # else:
                    #     self.nextAbsorber = air
                elif self.currentAbsorber == air and (self.checkIfOnScintillator(pt) or
                                                      any([self.checkIfOnAbsorber(ab, pt) for ab in def_absorbers])):
                    # elif self.currentAbsorber == air and (self.checkIfOnScintillator(pt)or
                    #                                       self.checkIfOnCollimator(pt) or self.checkIfOnWindow(pt)):
                    #     if self.checkIfOnScintillator(pt):
                    #         print(pt)

                    self.materialChange = True
                    # if self.checkIfOnWindow(pt):
                    #     self.nextAbsorber = al
                    # elif self.checkIfOnScintillator(pt):
                    #     self.nextAbsorber = scintillator
                    # else:
                    #     self.nextAbsorber = w
                # elif self.currentAbsorber == vacuum and (self.checkIfOnScintillator(pt) or
                #                                       self.checkIfOnCollimator(pt) or pt[2] <= (al.z + al.t)):
                #     self.materialChange = True
                else:
                    self.materialChange = False

                # intersectionPointsWithScintillatorTuples = [tuple(x) for x in intersectionPointsWithScintillator]

                self.checkNextAbsorber(pt, intersectionPointsWithScintillator)

                # if self.checkIfOnCollimator(pt):#and ((len(topPts)+len(botPts)+ len(discPts))>1):
                #     self.nextAbsorber = w
                #
                # elif pt in intersectionPointsWithScintillatorTuples:
                #     self.nextAbsorber = scintillator
                #
                # else:
                #     self.nextAbsorber = air

                # detector
                intersectionPointsWithDetectorTuples = [tuple(x) for x in intersectionPointsWithDetector]
                if self.materialChange == False:
                    if pt in intersectionPointsWithDetectorTuples:
                        # # x, y = Detector.checkPixel(-0.4096, pt[1])
                        # x, y = Detector.checkPixel(pt[0], pt[1])
                        # # if x is None or y is None:
                        # #     return peCount, totalInteractionCount
                        # pixelArray[y][x] += 1
                        # CCDPlotsX.append(pt[0])
                        # CCDPlotsY.append(pt[0])
                        # ptDisplay = vedo.Points([pt, currentPosition], r=10).color('red')
                        # ray = vedo.Line(currentPosition, pt, res=120).color('green')
                        # # vedo.show(CCDPlane, ray, ptDisplay, ScintillatorBlock, axes=True)
                        # totalInteractionCount += 1
                        # passDetectorCount += 1
                        return collimatorCount, RSCount  # peCount, totalInteractionCount


                # commented for the collimator sensitivity testing
                # if len(al.pts) != 0:
                #     if collimatorCount == None:
                #         print("test")
                #     collimatorCount = self.collimatorCountCheck(w, pt, collimatorCount)
                # if len(intersectionPointsWithWindow) != 0:
                #     if self.checkIfOnAbsorber(al, pt) and abs(
                #             pt[2] - (al.z - al.t / 2)) < 1e-04 and self.collimatorCounted == False:
                #         # if self.checkIfOnAbsorber(al, pt) and abs(pt[2] - al.z) <1e-04 and self.collimatorCounted == False:
                #         # if self.checkIfOnWindow(pt) and abs(pt[2] - al.z) <1e-04 and self.collimatorCounted == False:
                #         # if currentPosition != (0,0,distanceSA):
                #         #     print('additional collimator count')
                #         collimatorCount += 1
                #         self.collimatorCounted = True



                # self.x = pt[0]
                # self.y = pt[1]
                # self.z = pt[2]

                self.x = pt[0]
                self.y = pt[1]
                self.z = pt[2]
                if self.materialChange:
                    if self.currentAbsorber != air:
                        self.currentAbsorber = air
                        self.x, self.y, self.z = currentPosition
                        # line = vedo.Line([(0,0,-0.5,),(0.09661432974832347, 0.056447345966953805, 1.1893493800941657)]).color('p')
                        # sctpt = vedo.Points([(self.x, self.y, self.z)], r=10).color('red')
                        # sctpt1 = vedo.Points([(0,0,-0.5,),(0.09661432974832347, 0.056447345966953805, 1.1893493800941657)], r=10).color('black')
                        # vedo.show(tcone, bcone, discCylinder, window, ScintillatorBlock, CCDPlane, sctpt)
                    else:
                        self.currentAbsorber = self.nextAbsorber

                # density = self.currentAbsorber.density
                return self.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                                 eventSequence, depenergies)  # , density)




            else:
                if self.checkIfInSimulationRegion(nextPosition):
                    nextAbsorber = self.getInWhichAbsorber(nextPosition)
                    self.materialChange = False
                    # for nextabsorber in def_absorbers:
                    #     self.checkIfOnAbsorber(nextabsorber, nextPosition)

                    if not self.checkIfOnAbsorber(self.currentAbsorber, nextPosition) and (
                            self.currentAbsorber != air) and (self.currentAbsorber != scintillator):

                        # if (self.currentAbsorber == w or self.currentAbsorber==wDisc) and not self.checkIfOnCollimator(nextPosition):
                        self.materialChange = True
                        # if nextPosition[2] <= (al.z + al.t):
                        self.currentAbsorber = air
                        # else:
                        #     self.currentAbsorber = vaccum
                    elif (self.currentAbsorber == scintillator) and not self.checkIfOnScintillator(nextPosition):
                        self.materialChange = True
                        # if nextPosition[2] <= (al.z + al.t):
                        self.currentAbsorber = air
                        # else:
                        #     self.currentAbsorber = vaccum
                    # elif (self.currentAbsorber == al) and not self.checkIfOnWindow(nextPosition):
                    #     self.materialChange = True
                    #     # if nextPosition[2] <= (al.z + al.t):
                    #     self.currentAbsorber = air
                    #     # else:
                    #     self.currentAbsorber = vaccum

                    elif self.currentAbsorber == air and self.checkIfOnScintillator(nextPosition):
                        self.materialChange = True
                        self.currentAbsorber = scintillator

                    # getInWhichAbsorber(self, currentPosition)
                    elif self.currentAbsorber == air and (air != nextAbsorber):
                        # elif self.currentAbsorber == air and self.checkIfOnCollimator(nextPosition) and not self.checkIfOnWindow(nextPosition):
                        self.materialChange = True
                        self.currentAbsorber = nextAbsorber
                    # elif self.currentAbsorber == air and self.checkIfOnWindow(nextPosition):
                    #     self.materialChange = True
                    #     self.currentAbsorber = al
                    # if self.currentAbsorber == vaccum and self.checkIfOnScintillator(nextPosition):
                    #     self.materialChange = True
                    #     self.currentAbsorber = scintillator
                    #
                    # if self.currentAbsorber == vaccum and self.checkIfOnCollimator(nextPosition):
                    #     self.materialChange = True
                    #     self.currentAbsorber = w

                    if self.materialChange == True:
                        self.x, self.y, self.z = currentPosition
                        # line = vedo.Line([nextPosition, currentPosition]).color('p')
                        # sctpt = vedo.Points([currentPosition, nextPosition], r=10).color('red')
                        # sctpt1 = vedo.Points(pathPts, r=6).color('black')
                        # vedo.show(line, sctpt, sctpt1, tcone, bcone, discCylinder)
                        return self.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                                         eventSequence, depenergies)  # , self.currentAbsorber.density)
                    # if not self.checkIfOnCollimator(nextPosition) and not self.checkIfOnScintillator(nextPosition) and not self.checkIfOnWindow(nextPosition):
                    # outsideCount += 1 ## need to add other condition like in scintillator
                    # if self.checkIfInSimulationRegion(nextPosition) and self.z > 0:
                    # line = vedo.Line([nextPosition, currentPosition]).color('p')
                    # sctpt = vedo.Points([currentPosition, nextPosition], r=10).color('red')
                    # sctpt1 = vedo.Points(pathPts, r=6).color('black')
                    # vedo.show(sctpt, sctpt1, tcone, bcone, discCylinder, ScintillatorBlock, window)
                    # print('insideSimulation but not in the absorbers') # happens in very few possbility but still exists when samping P is very close to 1.
                    # totalInteractionCount += 1
                    # return peCount, totalInteractionCount
                    absorber = self.currentAbsorber
                    if self.currentAbsorber == scintillator:
                        if self.currentAbsorber.energies is None or self.energy < 1:
                            CSpeCo, CScsCo, CSrsCo, CStotalCo, IpeCo, IcsCo, IrsCo, ItotalCo = self.getCoefficients(
                                self.currentAbsorber)
                            # print('energy < 1')
                        else:
                            CSpeCo, CScsCo, CSrsCo, CStotalCo, IpeCo, IcsCo, IrsCo, ItotalCo = self.getCoefficients_file(
                                self.currentAbsorber)
                        Pcs = (detector_sub2.mass * detector.atom2_weight / scintillator.mass) * (
                            CStotalCo) / self.totalCo
                        r = random.uniform(0, 1)
                        if r <= Pcs:
                            absorber = detector_sub2
                            if absorber.energies is None or self.energy < 1:
                                self.getCoefficients(absorber)
                                # print('energy < 1')
                            else:
                                self.getCoefficients_file(absorber)
                        else:
                            absorber = detector_sub1
                            if absorber.energies is None or self.energy < 1:
                                self.getCoefficients(absorber)
                                # print('energy < 1')
                            else:
                                self.getCoefficients_file(absorber)
                        # if absorber.JK == 6.039:
                        #     test_total_cos = []
                        #     test_K_cos = []
                        #     test_L1_cos = []
                        #     test_L2_cos = []
                        #     test_L3_cos = []
                        #     test_M_cos = []
                        #     for test_en in np.arange(1, 142):
                        #         totalCos = self.getCoefficients_file_test(absorber, test_en)
                        #         test_total_cos.append(totalCos)
                        #         if test_en > absorber.KEDGE:
                        #             kCo = (absorber.JK - 1) / absorber.JK * totalCos
                        #             l1Co = (absorber.JL1 - 1) / absorber.JL1 * (totalCos - kCo)
                        #             l2Co = (absorber.JL2 - 1) / absorber.JL2 * (totalCos - kCo -l1Co)
                        #             l3Co = (absorber.JL3 - 1) / absorber.JL3 * (totalCos - kCo -l1Co - l2Co)
                        #             mCo = totalCos - l1Co - l2Co - l3Co - kCo
                        #
                        #             # pK = 1 - math.exp(-kCo * absorber.density * self.s)  # np.exp(-Coefficient.KPE*self.s)
                        #         # else:
                        #         #     kCo = 0
                        #
                        #         elif test_en > absorber.L1EDGE:
                        #             kCo = 0
                        #             l1Co = (absorber.JL1 - 1) / absorber.JL1 * totalCos
                        #             l2Co = (absorber.JL2 - 1) / absorber.JL2 * (totalCos - l1Co)
                        #             l3Co = (absorber.JL3 - 1) / absorber.JL3 * (totalCos - l1Co - l2Co)
                        #             mCo = totalCos - l1Co - l2Co - l3Co
                        #
                        #
                        #         # test_L1_cos.append(l1Co)
                        #         elif test_en > absorber.L2EDGE:
                        #             kCo = 0
                        #             l1Co = 0
                        #             l2Co = (absorber.JL2 - 1) / absorber.JL2 * totalCos
                        #             l3Co = (absorber.JL3 - 1) / absorber.JL3 * (totalCos - l2Co)
                        #             mCo = totalCos - l2Co - l3Co
                        #
                        #         # else:
                        #         #     l2Co = 0
                        #         # test_L2_cos.append(l2Co)
                        #         elif test_en > absorber.L3EDGE:
                        #             kCo = 0
                        #             l1Co = 0
                        #             l2Co = 0
                        #             l3Co = (absorber.JL3 - 1) / absorber.JL3 * totalCos
                        #             mCo = totalCos - l3Co
                        #
                        #         else:
                        #
                        #
                        #             mCo = totalCos
                        #             kCo = 0
                        #             l1Co = 0
                        #             l2Co = 0
                        #             l3Co = 0
                        #
                        #         # test_L3_cos.append(l3Co)
                        #         test_K_cos.append(kCo)
                        #
                        #         test_L1_cos.append(l1Co)
                        #
                        #         test_L2_cos.append(l2Co)
                        #
                        #         test_L3_cos.append(l3Co)
                        #
                        #         test_M_cos.append(mCo)
                        #
                        #     energies = np.arange(1, 142)
                        #     test_K_cos = np.array(test_K_cos)
                        #     test_L1_cos = np.array(test_L1_cos)
                        #     test_L2_cos = np.array(test_L2_cos)
                        #     test_L3_cos = np.array(test_L3_cos)
                        #     test_M_cos = np.array(test_M_cos)
                        #
                        #     plt.plot(energies, test_total_cos, color='red', label='total')
                        #     plt.plot(energies, test_K_cos+test_L1_cos+test_L2_cos+test_L3_cos, color='blue', label='L3+L2+L1+K')
                        #     plt.plot(energies, test_L1_cos+test_L2_cos+test_L3_cos, color='yellow', label='L3+L2+L1')
                        #     plt.plot(energies, test_L2_cos+test_L3_cos, color='green', label='L3+L2')
                        #     plt.plot(energies, test_L3_cos, color='pink', label='L3')
                        #     plt.plot(energies, test_M_cos, color='purple', label='M')
                        #     plt.plot(energies, test_K_cos+test_L1_cos+test_L2_cos+test_L3_cos+test_M_cos , color='black', label='K+L1+L2+L3+M')
                        #     plt.legend()
                        #     plt.show()


                    inten_type = self.predictInteraction(
                        absorber.density)  # self.currentAbsorber not used to avoid CS or I absorber alone

                    if self.currentAbsorber == scintillator and self.checkIfOnScintillator(
                            nextPosition):  ### record gamma events in the scitnialltor
                        # f.write(str(nextPosition) + ' ' + str(self.energy) + ' ' + str(
                        #     type) + '\n') ## wrong energy recorded place, need to recorded after interaction.

                        if inten_type == 'RS':
                            RSCount += 1
                        else:
                            x, y = scintillator.checkPixel(nextPosition[0], nextPosition[1])
                            # print(nextPosition[0], nextPosition[1], nextPosition[2], 'DetectorCount+1', type)
                            pixelArray[y][x] += 1
                            eventSequence.append((x, y))

                    # if absorber == air and nextPosition[2]> al.t+al.z and self.checkIfInSimulationRegion(nextPosition):
                    # print(pathPts)
                    # return self.main(traLengths) # current air actually is in the vaccum
                    ###testing for scintillator CScaterring
                    # if self.currentAbsorber == scintillator:
                    #     r = random.uniform(0, 1)
                    #     if r <= Pcs:
                    #         absorber = Absorber.Cs()
                    #         self.getCoefficients(absorber)
                    #     else:
                    #         absorber = Absorber.I()
                    #         self.getCoefficients(absorber)
                    #     type = 'CS'
                    ## below for collimator sensitivity testing
                    # else:
                    #     print(self.currentAbsorber)
                    if inten_type == 'PE':
                        # peCount += 1
                        # totalInteractionCount += 1
                        if self.currentAbsorber == air or (
                                not any([self.checkIfOnAbsorber(ab, nextPosition) for ab in def_absorbers])
                                and not (self.checkIfOnScintillator(nextPosition))):  # including vaccum?
                            shell = False  # should be self.peInteract(air)  currently the secondary effects in Air are not considered
                            # airPECount += 1
                            absorber = air
                            if self.checkIfInSimulationRegion(nextPosition) and VAC_Z is not None:
                                # airPECount -= airPECount
                                if  nextPosition[2] > VAC_Z:
                                # peCount -= peCount
                                    return self.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                                                     eventSequence, depenergies)
                                else:
                                    return collimatorCount, RSCount

                            # if self.checkIfInSimulationRegion(nextPosition):
                            #     sctpt = vedo.Points([currentPosition, nextPosition], r=7).color('yellow')
                            #     print(len(pathPts))
                            #
                            #     print((self.checkIfOnCollimator(nextPosition)))
                            #     sctpt1 = vedo.Points(pathPts, r=6).color('black')
                            #     line = vedo.Line(pathPts).color('p')
                            #     vedo.show(bcone, tcone, discCylinder, sctpt1, sctpt, ScintillatorBlock,line, axes=
                            #     True)
                            #     shell = # should be self.peInteract(air) but currently the secondary effects in Air are not considered
                            # else:
                            return collimatorCount, RSCount  # peCount, totalInteractionCount
                        # elif self.getInWhichAbsorber(nextPosition):
                        # # elif(self.checkIfOnCollimator(nextPosition)) and not self.checkIfOnWindow(nextPosition):
                        # # else:
                        #     shell = self.peInteract(wDisc) # comment now for collimator checking
                        #     absorber = wDisc
                        # elif self.checkIfOnWindow(nextPosition):
                        #     # else:
                        #     shell = self.peInteract(al)  # comment now for collimator checking
                        #     absorber = al
                        elif self.currentAbsorber == scintillator:
                            # elif self.checkIfOnScintillator(nextPosition):

                            shell = self.peInteract(
                                absorber)  ## self.currentAbsorber not used to avoid CS or I absorber alone
                        else:
                            nextAbs = self.getInWhichAbsorber(nextPosition)
                            if nextAbs != air:
                                # elif(self.checkIfOnCollimator(nextPosition)) and not self.checkIfOnWindow(nextPosition):
                                # else:
                                shell = self.peInteract(nextAbs)  # comment now for collimator checking
                                absorber = nextAbs
                        # else:
                        # print('error checking')
                        # else:
                        #     test = self.checkIfOnCollimator(nextPosition)
                        #     print(test)
                        #     test = self.checkIfOnCollimator(currentPosition)
                        #     print(test)
                        #     # line = vedo.Line([nextPosition, currentPosition]).color('p')
                        #     # sctpt = vedo.Points([currentPosition, nextPosition], r=5).color('red')
                        #     # sctpt1 = vedo.Points(pathPts, r=6).color('black')
                        #     # vedo.show(bcone, tcone, discCylinder, sctpt1, sctpt, ScintillatorBlock, line, axes=
                        #     #     True)
                        #     print(self.currentAbsorber)

                        # comment now for collimator checking
                        if secondary == True:

                            # try:
                            energy = self.energy
                            if shell != False:
                                # if type(absorber) == Absorber_flex.Detector or issubclass(type(absorber),
                                #                                                           Absorber_flex.Detector):
                                #     print(energy)

                                self.secondartEffects(absorber, shell, pixelArray, collimatorCount, RSCount,
                                                             eventSequence, depenergies)
                                if type(absorber) == Absorber_flex.Detector or issubclass(type(absorber),
                                                                                          Absorber_flex.Detector):
                                    # if energy != self.energy:
                                        depenergies.append(self.energy)
                                        # if energy == self.energy:
                                        #     print(energy)
                                    # else:
                                    #     depenergies.append(energy)
                            else:
                                if type(absorber) == Absorber_flex.Detector or issubclass(type(absorber),
                                                                                          Absorber_flex.Detector):
                                    depenergies.append(energy)
                                return collimatorCount, RSCount


                        elif type(absorber) == Absorber_flex.Detector or issubclass(type(absorber), Absorber_flex.Detector):
                            depenergies.append(self.energy)
                        # except Exception as e:
                        #     print(e)

                        return collimatorCount, RSCount  # peCount, totalInteractionCount

                    elif inten_type == "CS":  ## compton scattering
                        # totalInteractionCount += 1
                        # self.Cscattered =True
                        # csCount += 1
                        # totalInteractionCount += 1
                        if self.checkIfInSimulationRegion(nextPosition):
                            if (not any([self.checkIfOnAbsorber(ab, nextPosition) for ab in def_absorbers])
                                    and not (self.checkIfOnScintillator(nextPosition))):
                                # if nextPosition[2] <= al.z+al.t:
                                absorber = air
                                if VAC_Z is not None and  nextPosition[2] > VAC_Z:

                                    # csCount -= csCount
                                    return self.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount,
                                                     RSCount, eventSequence, depenergies)

                            elif self.checkIfOnScintillator(
                                    nextPosition):  ## maybe replace with checking with currentAbsorber and keep absorber from last type prediction
                                absorber = scintillator
                            else:
                                nextAbs = self.getInWhichAbsorber(nextPosition)
                                if nextAbs != air:
                                    absorber = nextAbs
                            return self.csInteract(traLengths, pathPts, absorbers, pixelArray, collimatorCount,
                                                   RSCount, eventSequence, depenergies, absorber)
                        else:
                            # csCount += 1
                            # totalInteractionCount += 1 # moved to the begining for the consistency checking
                            return collimatorCount, RSCount  # peCount, totalInteractionCount

                    else:

                        if self.checkIfInSimulationRegion(nextPosition):
                            if (not any([self.checkIfOnAbsorber(ab, nextPosition) for ab in def_absorbers])
                                    and not (self.checkIfOnScintillator(nextPosition))):
                                absorber = air
                            elif self.checkIfOnScintillator(
                                    nextPosition):  ## maybe replace with checking with currentAbsorber and keep absorber from last type prediction
                                absorber = scintillator
                            else:
                                nextAbs = self.getInWhichAbsorber(nextPosition)
                                if nextAbs != air:
                                    absorber = nextAbs
                            if absorber != air:
                                # absorber = self.currentAbsorber
                                # return self.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount)#, absorber.density)  #if not return doesnt affect the whole function except for the totalInteractionCount
                                return self.rsInteract(traLengths, pathPts, absorbers, pixelArray, collimatorCount,
                                                       RSCount, eventSequence, depenergies)
                            else:
                                return self.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                                                 eventSequence, depenergies)

                        else:
                            # totalInteractionCount += 1
                            # otherCount += 1

                            return collimatorCount, RSCount  # peCount, totalInteractionCount
                else:
                    return collimatorCount, RSCount

    def setDirection(self, theta, phi):
        self.theta = theta
        self.phi = phi

    def setStepSize(self, s):
        self.s = s

    def setEnergy(self, energy):
        self.energy = energy

    def moveStep(self):
        # dr = 0.01
        dx = self.s * np.sin(self.theta) * np.cos(self.phi)
        dy = self.s * np.sin(self.theta) * np.sin(self.phi)
        dz = self.s * np.cos(self.theta)
        self.x = self.x + dx
        self.y = self.y + dy
        self.z = self.z + dz

    # def insideAbsorber(self):  # check if the current photon is inside the material
    #     if self.z <= w.z and self.z >= w.z - w.t and np.abs(self.y) <= 150 and np.abs(self.x) <= 150:
    #         self.s = 0.005
    #         return True
    #     else:
    #         return False

    def getCoefficients_file(self, currentAbsorber):
        # if not ifInAir:
        #     energy = self.energy*1000 #kev->ev
        if (type(currentAbsorber) is not Absorber_flex.Air) and (type(currentAbsorber) is not Absorber_flex.Detector):
            self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, currentAbsorber.energies,
                                                                                 currentAbsorber.csCoLs,
                                                                                 currentAbsorber.rsCoLs,
                                                                                 currentAbsorber.peCoLs,
                                                                                 currentAbsorber.totalCoLs)

        elif type(currentAbsorber) is Absorber_flex.Detector:
            self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, currentAbsorber.energies,
                                                                                 currentAbsorber.csCoLs,
                                                                                 currentAbsorber.rsCoLs,
                                                                                 currentAbsorber.peCoLs,
                                                                                 currentAbsorber.totalCoLs)
            # self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, 'CsIdata.txt') # for W
            # CSpeCo, CScsCo, CSrsCo, CStotalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
            CSpeCo, CScsCo, CSrsCo, CStotalCo = dataGrid.interpolate(self.energy, detector_sub2.energies,
                                                                     detector_sub2.csCoLs, detector_sub2.rsCoLs,
                                                                     detector_sub2.peCoLs, detector_sub2.totalCoLs)
            # dataGrid.interpolate(self.energy, 'Csdata.txt')  # for Cs
            # IpeCo, IcsCo, IrsCo, ItotalCo = xraydb.mu_elam('I', energy, kind='photo'), xraydb.mu_elam('I', energy, kind='incoh'),xraydb.mu_elam('I', energy, kind='coh'), xraydb.mu_elam('I', energy, kind='total')
            IpeCo, IcsCo, IrsCo, ItotalCo = dataGrid.interpolate(self.energy, detector_sub1.energies,
                                                                 detector_sub1.csCoLs, detector_sub1.rsCoLs,
                                                                 detector_sub1.peCoLs, detector_sub1.totalCoLs)
            # dataGrid.interpolate(self.energy, 'Idata.txt')  # for I
            return CSpeCo, CScsCo, CSrsCo, CStotalCo, IpeCo, IcsCo, IrsCo, ItotalCo
        elif issubclass(type(currentAbsorber), Absorber_flex.Detector):
            if type(currentAbsorber) is Absorber_flex.Detector_sub:
                # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
                # self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Csenergies,
                #                                                                  CscsCoLs, CsrsCoLs, CspeCoLs,
                #                                                                  CstotalCoLs)
                self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy,
                                                                                     currentAbsorber.energies,
                                                                                     currentAbsorber.csCoLs,
                                                                                     currentAbsorber.rsCoLs,
                                                                                     currentAbsorber.peCoLs,
                                                                                     currentAbsorber.totalCoLs)

                # dataGrid.interpolate(self.energy, 'Csdata.txt')
                return
        else:
            self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Airenergies, AircsCoLs,
                                                                                 AirrsCoLs, AirpeCoLs, AirtotalCoLs)
            # dataGrid.interpolate(self.energy, 'Airdata.txt') # for air
            return


    def getCoefficients_file_test(self, currentAbsorber, energy):

        if (type(currentAbsorber) is not Absorber_flex.Air) and (type(currentAbsorber) is not Absorber_flex.Detector):
            peCo, csCo, rsCo, totalCo = dataGrid.interpolate(energy, currentAbsorber.energies,
                                                                                 currentAbsorber.csCoLs,
                                                                                 currentAbsorber.rsCoLs,
                                                                                 currentAbsorber.peCoLs,
                                                                                 currentAbsorber.totalCoLs)
            return totalCo
        elif type(currentAbsorber) is Absorber_flex.Detector:
            peCo, csCo, rsCo, totalCo = dataGrid.interpolate(energy, currentAbsorber.energies,
                                                                                 currentAbsorber.csCoLs,
                                                                                 currentAbsorber.rsCoLs,
                                                                                 currentAbsorber.peCoLs,
                                                                                 currentAbsorber.totalCoLs)
            # self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, 'CsIdata.txt') # for W
            # CSpeCo, CScsCo, CSrsCo, CStotalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
            CSpeCo, CScsCo, CSrsCo, CStotalCo = dataGrid.interpolate(energy, detector_sub2.energies,
                                                                     detector_sub2.csCoLs, detector_sub2.rsCoLs,
                                                                     detector_sub2.peCoLs, detector_sub2.totalCoLs)
            # dataGrid.interpolate(self.energy, 'Csdata.txt')  # for Cs
            # IpeCo, IcsCo, IrsCo, ItotalCo = xraydb.mu_elam('I', energy, kind='photo'), xraydb.mu_elam('I', energy, kind='incoh'),xraydb.mu_elam('I', energy, kind='coh'), xraydb.mu_elam('I', energy, kind='total')
            IpeCo, IcsCo, IrsCo, ItotalCo = dataGrid.interpolate(energy, detector_sub1.energies,
                                                                 detector_sub1.csCoLs, detector_sub1.rsCoLs,
                                                                 detector_sub1.peCoLs, detector_sub1.totalCoLs)
            # dataGrid.interpolate(self.energy, 'Idata.txt')  # for I
            return totalCo
        elif issubclass(type(currentAbsorber), Absorber_flex.Detector):
            if type(currentAbsorber) is Absorber_flex.Detector_sub:
                # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
                # self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Csenergies,
                #                                                                  CscsCoLs, CsrsCoLs, CspeCoLs,
                #                                                                  CstotalCoLs)
               peCo, csCo, rsCo, totalCo = dataGrid.interpolate(energy,currentAbsorber.energies,
                                                                                     currentAbsorber.csCoLs,
                                                                                     currentAbsorber.rsCoLs,
                                                                                     currentAbsorber.peCoLs,
                                                                                     currentAbsorber.totalCoLs)

            return totalCo
        else:
            peCo, csCo, rsCo, totalCo = dataGrid.interpolate(energy, Airenergies, AircsCoLs,
                                                                                 AirrsCoLs, AirpeCoLs, AirtotalCoLs)
            # dataGrid.interpolate(self.energy, 'Airdata.txt') # for air
            return totalCo


    def getCoefficients(self, currentAbsorber):
        # if not ifInAir:
        #     energy = self.energy*1000 #kev->ev

        if (type(currentAbsorber) is not Absorber_flex.Air) and (type(currentAbsorber) is not Absorber_flex.Detector):
            # try:
            self.peCo, self.csCo, self.rsCo, self.totalCo = xraylib.CS_Photo(currentAbsorber.atomNum,
                                                                             self.energy), xraylib.CS_Compt(
                currentAbsorber.atomNum, self.energy), xraylib.CS_Rayl(currentAbsorber.atomNum,
                                                                       self.energy), xraylib.CS_Total(
                currentAbsorber.atomNum, self.energy)
            return
        # except:
        #     print(currentAbsorber)
        # if type(currentAbsorber) is Absorber.W:
        #     # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('W', energy, kind='photo'), xraydb.mu_elam('W', energy, kind='incoh'),xraydb.mu_elam('W', energy, kind='coh'), xraydb.mu_elam('W', energy, kind='total')
        #     self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Wenergies,
        #                                                                          WcsCoLs, WrsCoLs, WpeCoLs,
        #                                                                          WtotalCoLs)
        #     #dataGrid.interpolate(self.energy, 'data.txt') # for W
        #     return
        # elif type(currentAbsorber) is Absorber.Al:
        #     # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('Al', energy, kind='photo'), xraydb.mu_elam('Al', energy, kind='incoh'),xraydb.mu_elam('Al', energy, kind='coh'), xraydb.mu_elam('Al', energy, kind='total')
        #     self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Alenergies,
        #                                                                          AlcsCoLs, AlrsCoLs, AlpeCoLs,
        #                                                                          AltotalCoLs)
        #     # dataGrid.interpolate(self.energy, 'Aldata.txt')  # for Al
        #     return

        elif type(currentAbsorber) is Absorber_flex.Detector:
            self.peCo, self.csCo, self.rsCo, self.totalCo = xraylib.CS_Photo_CP(currentAbsorber.detector_material,
                                                                                self.energy), xraylib.CS_Compt_CP(
                currentAbsorber.detector_material, self.energy), xraylib.CS_Rayl_CP(currentAbsorber.detector_material,
                                                                                    self.energy), xraylib.CS_Total_CP(
                currentAbsorber.detector_material, self.energy)
            # dataGrid.interpolate(self.energy, CsIenergies, CsIcsCoLs, CsIrsCoLs, CsIpeCoLs, CsItotalCoLs)
            # self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, 'CsIdata.txt') # for W
            # CSpeCo, CScsCo, CSrsCo, CStotalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
            CSpeCo, CScsCo, CSrsCo, CStotalCo = xraylib.CS_Photo(currentAbsorber.atom2, self.energy), xraylib.CS_Compt(
                currentAbsorber.atom2, self.energy), xraylib.CS_Rayl(currentAbsorber.atom2,
                                                                     self.energy), xraylib.CS_Total(
                currentAbsorber.atom2, self.energy)
            # dataGrid.interpolate(self.energy, 'Csdata.txt')  # for Cs
            # IpeCo, IcsCo, IrsCo, ItotalCo = xraydb.mu_elam('I', energy, kind='photo'), xraydb.mu_elam('I', energy, kind='incoh'),xraydb.mu_elam('I', energy, kind='coh'), xraydb.mu_elam('I', energy, kind='total')
            IpeCo, IcsCo, IrsCo, ItotalCo = xraylib.CS_Photo(currentAbsorber.atom1, self.energy), xraylib.CS_Compt(
                currentAbsorber.atom1, self.energy), xraylib.CS_Rayl(currentAbsorber.atom1,
                                                                     self.energy), xraylib.CS_Total(
                currentAbsorber.atom1, self.energy)
            # dataGrid.interpolate(self.energy, 'Idata.txt')  # for I
            return CSpeCo, CScsCo, CSrsCo, CStotalCo, IpeCo, IcsCo, IrsCo, ItotalCo

        # elif type(currentAbsorber) is Absorber.CsI:
        #     self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, CsIenergies, CsIcsCoLs, CsIrsCoLs, CsIpeCoLs, CsItotalCoLs)
        #     # self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, 'CsIdata.txt') # for W
        #     # CSpeCo, CScsCo, CSrsCo, CStotalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
        #     CSpeCo, CScsCo, CSrsCo, CStotalCo = dataGrid.interpolate(self.energy, Csenergies,
        #                                                                          CscsCoLs, CsrsCoLs, CspeCoLs,
        #                                                                          CstotalCoLs)
        #     #dataGrid.interpolate(self.energy, 'Csdata.txt')  # for Cs
        #     # IpeCo, IcsCo, IrsCo, ItotalCo = xraydb.mu_elam('I', energy, kind='photo'), xraydb.mu_elam('I', energy, kind='incoh'),xraydb.mu_elam('I', energy, kind='coh'), xraydb.mu_elam('I', energy, kind='total')
        #     IpeCo, IcsCo, IrsCo, ItotalCo = dataGrid.interpolate(self.energy, Ienergies,
        #                                                              IcsCoLs, IrsCoLs, IpeCoLs,
        #                                                              ItotalCoLs)
        #     #dataGrid.interpolate(self.energy, 'Idata.txt')  # for I
        #     return CSpeCo, CScsCo, CSrsCo, CStotalCo, IpeCo, IcsCo, IrsCo, ItotalCo

        elif issubclass(type(currentAbsorber), Absorber_flex.Detector):
            if type(currentAbsorber) is Absorber_flex.Detector_sub:
                # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
                # self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Csenergies,
                #                                                                  CscsCoLs, CsrsCoLs, CspeCoLs,
                #                                                                  CstotalCoLs)
                self.peCo, self.csCo, self.rsCo, self.totalCo = xraylib.CS_Photo(currentAbsorber.atomNum,
                                                                                 self.energy), xraylib.CS_Compt(
                    currentAbsorber.atomNum, self.energy), xraylib.CS_Rayl(currentAbsorber.atomNum,
                                                                           self.energy), xraylib.CS_Total(
                    currentAbsorber.atomNum, self.energy)

                # dataGrid.interpolate(self.energy, 'Csdata.txt')
                return
            # else:
            #     # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('I', energy, kind='photo'), xraydb.mu_elam('I', energy, kind='incoh'),xraydb.mu_elam('I', energy, kind='coh'), xraydb.mu_elam('I', energy, kind='total')
            #     self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Ienergies,
            #                                                          IcsCoLs, IrsCoLs, IpeCoLs,
            #                                                          ItotalCoLs)
            #     #dataGrid.interpolate(self.energy, 'Idata.txt')
            #     return
            #
        # elif issubclass(type(currentAbsorber), Absorber.CsI):
        #     if type(currentAbsorber) is Absorber.Cs:
        #         # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
        #         self.peCo, self.csCo, self.rsCo, self.totalCo =dataGrid.interpolate(self.energy, Csenergies,
        #                                                                          CscsCoLs, CsrsCoLs, CspeCoLs,
        #                                                                          CstotalCoLs)
        #         #dataGrid.interpolate(self.energy, 'Csdata.txt')
        #         return
        #     else:
        #         # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('I', energy, kind='photo'), xraydb.mu_elam('I', energy, kind='incoh'),xraydb.mu_elam('I', energy, kind='coh'), xraydb.mu_elam('I', energy, kind='total')
        #         self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Ienergies,
        #                                                              IcsCoLs, IrsCoLs, IpeCoLs,
        #                                                              ItotalCoLs)
        #         #dataGrid.interpolate(self.energy, 'Idata.txt')
        #         return
        else:
            self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Airenergies, AircsCoLs,
                                                                                 AirrsCoLs, AirpeCoLs, AirtotalCoLs)
            # dataGrid.interpolate(self.energy, 'Airdata.txt') # for air
            return

    # def getCoefficients(self, currentAbsorber):
    #     # if not ifInAir:
    #         energy = self.energy*1000 #kev->ev
    #
    #         if type(currentAbsorber) is not Absorber_flex.Air:
    #         if type(currentAbsorber) is Absorber.W:
    #             # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('W', energy, kind='photo'), xraydb.mu_elam('W', energy, kind='incoh'),xraydb.mu_elam('W', energy, kind='coh'), xraydb.mu_elam('W', energy, kind='total')
    #             self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Wenergies,
    #                                                                                  WcsCoLs, WrsCoLs, WpeCoLs,
    #                                                                                  WtotalCoLs)
    #             #dataGrid.interpolate(self.energy, 'data.txt') # for W
    #             return
    #         elif type(currentAbsorber) is Absorber.Al:
    #             # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('Al', energy, kind='photo'), xraydb.mu_elam('Al', energy, kind='incoh'),xraydb.mu_elam('Al', energy, kind='coh'), xraydb.mu_elam('Al', energy, kind='total')
    #             self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Alenergies,
    #                                                                                  AlcsCoLs, AlrsCoLs, AlpeCoLs,
    #                                                                                  AltotalCoLs)
    #             # dataGrid.interpolate(self.energy, 'Aldata.txt')  # for Al
    #             return
    #         elif type(currentAbsorber) is Absorber.CsI:
    #             self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, CsIenergies, CsIcsCoLs, CsIrsCoLs, CsIpeCoLs, CsItotalCoLs)
    #             # self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, 'CsIdata.txt') # for W
    #             # CSpeCo, CScsCo, CSrsCo, CStotalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
    #             CSpeCo, CScsCo, CSrsCo, CStotalCo = dataGrid.interpolate(self.energy, Csenergies,
    #                                                                                  CscsCoLs, CsrsCoLs, CspeCoLs,
    #                                                                                  CstotalCoLs)
    #             #dataGrid.interpolate(self.energy, 'Csdata.txt')  # for Cs
    #             # IpeCo, IcsCo, IrsCo, ItotalCo = xraydb.mu_elam('I', energy, kind='photo'), xraydb.mu_elam('I', energy, kind='incoh'),xraydb.mu_elam('I', energy, kind='coh'), xraydb.mu_elam('I', energy, kind='total')
    #             IpeCo, IcsCo, IrsCo, ItotalCo = dataGrid.interpolate(self.energy, Ienergies,
    #                                                                      IcsCoLs, IrsCoLs, IpeCoLs,
    #                                                                      ItotalCoLs)
    #             #dataGrid.interpolate(self.energy, 'Idata.txt')  # for I
    #             return CSpeCo, CScsCo, CSrsCo, CStotalCo, IpeCo, IcsCo, IrsCo, ItotalCo
    #         elif issubclass(type(currentAbsorber), Absorber.CsI):
    #             if type(currentAbsorber) is Absorber.Cs:
    #                 # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('Cs', energy, kind='photo'), xraydb.mu_elam('Cs', energy, kind='incoh'),xraydb.mu_elam('Cs', energy, kind='coh'), xraydb.mu_elam('Cs', energy, kind='total')
    #                 self.peCo, self.csCo, self.rsCo, self.totalCo =dataGrid.interpolate(self.energy, Csenergies,
    #                                                                                  CscsCoLs, CsrsCoLs, CspeCoLs,
    #                                                                                  CstotalCoLs)
    #                 #dataGrid.interpolate(self.energy, 'Csdata.txt')
    #                 return
    #             else:
    #                 # self.peCo, self.csCo, self.rsCo, self.totalCo = xraydb.mu_elam('I', energy, kind='photo'), xraydb.mu_elam('I', energy, kind='incoh'),xraydb.mu_elam('I', energy, kind='coh'), xraydb.mu_elam('I', energy, kind='total')
    #                 self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Ienergies,
    #                                                                      IcsCoLs, IrsCoLs, IpeCoLs,
    #                                                                      ItotalCoLs)
    #                 #dataGrid.interpolate(self.energy, 'Idata.txt')
    #                 return
    #         else:
    #             self.peCo, self.csCo, self.rsCo, self.totalCo = dataGrid.interpolate(self.energy, Airenergies, AircsCoLs, AirrsCoLs, AirpeCoLs, AirtotalCoLs)
    #             # dataGrid.interpolate(self.energy, 'Airdata.txt') # for air
    #             return

    def ifInteract(self):
        pInteract = np.exp(-self.totalCo * self.s)
        randomCompare = random.random()
        if randomCompare <= pInteract:
            return True
        else:
            return False

    def predictInteraction(self, density):
        # ppe = 1 - math.exp(-self.peCo*W.density*self.s) # cm^2/g * g/cm^3 * cm = 1
        # pcs = 1 - math.exp(-self.csCo*W.density*self.s)
        # prs = 1 - math.exp(-self.rsCo*W.density*self.s)
        # if self.insideCurrentAbsorberUpdated(W):
        ppe = 1 - math.exp(-self.peCo * density * self.s)  # cm^2/g * g/cm^3 * cm = 1
        pcs = 1 - math.exp(-self.csCo * density * self.s)
        prs = 1 - math.exp(-self.rsCo * density * self.s)
        # else:
        #     ppe = 1 - math.exp(-self.peCo*Air.density*self.s) # cm^2/g * g/cm^3 * cm = 1
        #     pcs = 1 - math.exp(-self.csCo*Air.density*self.s)
        #     prs = 1 - math.exp(-self.rsCo*Air.density*self.s)

        # ppe = (self.peCo/ W.density)/ (self.totalCo/W.density)
        # pcs = (self.csCo/ W.density)/ (self.totalCo/W.density)
        # prs = (self.rsCo/ W.density)/ (self.totalCo/W.density)
        psum = ppe + pcs + prs
        ppe = ppe / psum
        pcs = pcs / psum
        prs = prs / psum
        inten_type = Comparison.compareAttenuation(ppe, pcs, prs)
        # if type == 'PE':
        #     print(ppe, pcs, prs, type)
        return inten_type

    # def interact(self):
    #     if self.energy > 69.53:
    #         pPE = W.KAtt / W.totalK
    #         pCS = W.CSK / W.totalK
    #         pRS = W.RSK /W.totalK
    #         sum = pPE + pCS +pRS
    #         pPE = pPE/ sum
    #         pCS = pCS /sum
    #         pRS = pRS /sum
    #         type = Comparison.compareAttenuation(pPE, pCS, pRS)
    #         if type == "PE":
    #             self.PE(0.005)
    #             print('absorbed')
    #         else:
    #             print( "not absorbed")
    def peInteract(self, absorber):  # photoelectrc effect need to be updated with Jump factor
        # if type(absorber) == Absorber_flex.Detector or issubclass(type(absorber),
        #                                                           Absorber_flex.Detector) and self.energy == 141:
        #     print(absorber.JK)
        if self.energy > absorber.KEDGE:
            # kCo = (absorber.JK - 1) / absorber.JK * self.totalCo
            kCo = (absorber.JK - 1) / absorber.JK * self.totalCo
            l1Co = (absorber.JL1 - 1) / absorber.JL1 * (self.totalCo - kCo)
            l2Co = (absorber.JL2 - 1) / absorber.JL2 * (self.totalCo - kCo - l1Co)
            l3Co = (absorber.JL3 - 1) / absorber.JL3 * (self.totalCo - kCo - l1Co - l2Co)
            # m1Co = (absorber.JM1 - 1) / absorber.JM1 * (self.totalCo - kCo - l1Co - l2Co - l3Co)
            mCo = self.totalCo - l1Co - l2Co - l3Co -kCo
            # kCo = kCo  - bias
            # l1Co = l1Co - bias
            # l2Co = l2Co - bias
            # l3Co = l3Co - bias
            pK = 1 - math.exp(-kCo * absorber.density * self.s)  # np.exp(-Coefficient.KPE*self.s)
            pL1 = 1 - math.exp(-l1Co * absorber.density * self.s)  # pL1 = np.exp(-Coefficient.L1PE*self.s)
            pL2 = 1 - math.exp(-l2Co * absorber.density * self.s)
            pL3 = 1 - math.exp(-l3Co * absorber.density * self.s)
            pM = 1 - math.exp(-mCo * absorber.density * self.s)
        # else:
        #     pK = 0
        elif self.energy > absorber.L1EDGE:
            # kCo = 0
            l1Co = (absorber.JL1 - 1) / absorber.JL1 * self.totalCo
            l2Co = (absorber.JL2 - 1) / absorber.JL2 * (self.totalCo - l1Co)
            l3Co = (absorber.JL3 - 1) / absorber.JL3 * (self.totalCo - l1Co - l2Co)
            # m1Co = (absorber.JM1 - 1) / absorber.JM1 * (self.totalCo - l1Co - l2Co - l3Co)
            # bias = self.totalCo - m1Co - l1Co - l2Co - l3Co
            mCo = self.totalCo - l1Co - l2Co - l3Co
            # l1Co = l1Co - bias
            # l2Co = l2Co - bias
            # l3Co = l3Co - bias
            pK = 0
            pL1 = 1 - math.exp(-l1Co * absorber.density * self.s)  # pL1 = np.exp(-Coefficient.L1PE*self.s)
            pL2 = 1 - math.exp(-l2Co * absorber.density * self.s)
            pL3 = 1 - math.exp(-l3Co * absorber.density * self.s)
            pM = 1 - math.exp(-mCo * absorber.density * self.s)

            # l1Co = (absorber.JL1 - 1) / absorber.JL1 * self.totalCo
            # pL1 = 1 - math.exp(-l1Co * absorber.density * self.s)  # pL1 = np.exp(-Coefficient.L1PE*self.s)
        # else:
        #     pL1 = 0
        elif self.energy > absorber.L2EDGE:
            # l2Co = (absorber.JL2 - 1) / absorber.JL2 * self.totalCo
            # pL2 = 1 - math.exp(-l2Co * absorber.density * self.s)  # np.exp(-Coefficient.L2PE*self.s)
            # kCo = 0
            # l1Co = 0
            l2Co = (absorber.JL2 - 1) / absorber.JL2 * self.totalCo
            l3Co = (absorber.JL3 - 1) / absorber.JL3 * (self.totalCo - l2Co)
            mCo = self.totalCo - l2Co - l3Co
            # m1Co = (absorber.JM1 - 1) / absorber.JM1 * (self.totalCo - l2Co - l3Co)
            # bias = self.totalCo - m1Co - l2Co - l3Co
            # l2Co = l2Co - bias
            # l3Co = l3Co - bias
            pK = 0
            pL1 = 0
            pL2 = 1 - math.exp(-l2Co * absorber.density * self.s)
            pL3 = 1 - math.exp(-l3Co * absorber.density * self.s)
            pM = 1 - math.exp(-mCo * absorber.density * self.s)
        # else:
        #     pL2 = 0
        elif self.energy > absorber.L3EDGE:
            # kCo = 0
            # l1Co = 0
            # l2Co = 0
            l3Co = (absorber.JL3 - 1) / absorber.JL3 * self.totalCo
            mCo = self.totalCo - l3Co
            # m1Co = (absorber.JM1 - 1) / absorber.JM1 * (self.totalCo - l3Co)
            # bias = self.totalCo - m1Co - l3Co
            # l3Co = l3Co - bias
            # l3Co = (absorber.JL3 - 1) / absorber.JL3 * self.totalCo
            pK = 0
            pL1 = 0
            pL2 = 0
            pL3 = 1 - math.exp(-l3Co * absorber.density * self.s)  # pL3 = np.exp(-Coefficient.L3PE*self.s)
            pM = 1 - math.exp(-mCo * absorber.density * self.s)
        # else:
        #     pL3 = 0
        else:
            # kCo = 0
            # l1Co = 0
            # l2Co = 0
            # l3Co = 0
            # mCo = self.totalCo
            pK = 0
            pL1 = 0
            pL2 = 0
            pL3 = 0
            pM = 1



        sumPro = (pK + pL1 + pL2 + pL3+pM)
        # if sumPro != 0:
        if pM != 1:
            pK = pK / sumPro
            pL1 = pL1 / sumPro
            pL2 = pL2 / sumPro
            pL3 = pL3 / sumPro
            # pM = pM / sumPro

            # pList = {pK, pL1, pL2, pL3}
            shell = Comparison.compareShell(pK, pL1, pL2, pL3)
            return shell
        else:
            return False

    def predictSecondEffect(self, absorber, shell):
        ranCompare = random.random()
        if shell == 'K':
            # fyK = 0.958
            fyK = absorber.fyK
            if ranCompare <= fyK:
                return 'Fluo'
            else:
                return 'Non-rad'
        elif shell == 'L1':
            # fyL1 = 0.147
            fyL1 = absorber.fyL1
            if ranCompare <= fyL1:
                return 'Fluo'
            else:
                return 'Non-rad'
        elif shell == 'L2':
            # fyL2 = 0.270
            fyL2 = absorber.fyL2
            if ranCompare <= fyL2:
                return 'Fluo'
            else:
                return 'Non-rad'
        elif shell == 'L3':
            # fyL3 = 0.255
            fyL3 = absorber.fyL3
            if ranCompare <= fyL3:
                return 'Fluo'
            else:
                return 'Non-rad'

    def fluo(self, absorber, shell):  # create a new photon
        vacancy = False
        if shell == 'K':
            # kA1 = 0.5024
            # kA2 = 0.2915
            # kB1 = 0.1106
            # kB2 = 0.0402
            # kB3 = 0.0553
            # alpha1 causes a vacancy in L3
            kA1 = absorber.kA1
            # alpha2 causes a vacancy in L2
            kA2 = absorber.kA2
            kB1 = absorber.kB1
            kB2 = absorber.kB2
            kB3 = absorber.kB3
            tran = Comparison.compareKtransition(kA1, kA2, kB1, kB3, kB2)
        elif shell == 'L1':
            l1G3 = absorber.l1G3
            l1B3 = absorber.l1B3
            l1B4 = absorber.l1B4
            tran = Comparison.compareL1transition(l1G3, l1B3, l1B4)
        elif shell == 'L2':
            l2B1 = absorber.l2B1
            l2G1 = absorber.l2G1
            tran = Comparison.compareL2transition(l2B1, l2G1)
        elif shell == 'L3':
            l3L = absorber.l3L
            l3A2 = absorber.l3A2
            l3A1 = absorber.l3A1
            l3B2 = absorber.l3B2
            tran = Comparison.compareL3transition(l3L, l3A2, l3A1, l3B2)
        energyLine = absorber.fluoTranEnergies[tran]
        # if type(absorber) == Absorber_flex.Detector or issubclass(type(absorber), Absorber_flex.Detector):
        #     depenergies.append(self.energy - energyLine)
        if energyLine >= 1:
            fluPhoton = Photon(self.x, self.y, self.z, energyLine)
        else:
            fluPhoton = None
        if tran == 'Kalpha1':
            vacancy = 'L3'
        if tran == 'Kalpha2':
            vacancy = 'L2'
        return fluPhoton, vacancy

    def nonRad(self, absorber, shell):
        if shell == 'K':
            # kL1L1 = 0.1
            # kL1L2 = 0.148
            # kL1L3 = 0.102
            # kL2L2 = 0.01
            # kL2L3 = 0.188
            # kL3L3 = 0.086
            # kL1X = 0.123
            # kL2X = 0.089
            # kL3X = 0.12
            # kXX = 0.034
            kL1L1 = absorber.kL1L1
            kL1L2 = absorber.kL1L2
            kL1L3 = absorber.kL1L3
            kL2L2 = absorber.kL2L2
            kL2L3 = absorber.kL2L3
            kL3L3 = absorber.kL3L3
            kL1X = absorber.kL1X
            kL2X = absorber.kL2X
            kL3X = absorber.kL3X
            kXX = absorber.kXX
            tran = Comparison.compareNon_RaKtransition(kL1L1, kL1L2, kL1L3, kL2L2, kL2L3,
                                                       kL3L3, kL1X, kL2X, kL3X, kXX)
        elif shell == 'L1':
            # l1L2X = 0.2
            # l1L3X = 0.382
            # l1XX = 0.472
            l1L2X = absorber.l1L2X
            l1L3X = absorber.l1L3X
            l1XX = absorber.l1XX
            rateSum = l1L2X + l1L3X + l1XX
            l1L2X = l1L2X / rateSum
            l1L3X = l1L3X / rateSum
            l1XX = l1XX / rateSum
            tran = Comparison.compareNon_RaL1transition(l1L2X, l1L3X, l1XX)
        elif shell == 'L2':
            # l2L3X = 0.182
            # l2XX = 0.818
            l2L3X = absorber.l2L3X
            l2XX = absorber.l2XX
            tran = Comparison.compareNon_RaL2transition(l2L3X, l2XX)
        elif shell == 'L3':
            tran = Comparison.compareNon_RaL3transition()
        # else:
        #     print(shell)
        return tran

    def secondartEffects(self, absorber, shell, pixelArray, collimatorCount, RSCount, eventSequence, depenergies):  # w, wDisc
        # global peCount, totalInteractionCount, fluoPhotonCount, augerCount
        # secondartEffects(absorber, shell, pinhole, disc)
        # if type(absorber) is Absorber_flex.Detector_sub:
        #     print(self.energy)
        secondaryEffect = self.predictSecondEffect(absorber, shell)
        if secondaryEffect == 'Fluo':
            # print('self', self.energy)
            photon, flouVacancy = self.fluo(absorber, shell)
            # print('fluPhoton', photon.energy)

            # if self.currentAbsorber == scintillator:
            #     print(self.energy-photon.energy, 'PE Flu Deposit', photon.x, photon.y, photon.z)

            # fluoPhotonCount += 1
            if photon is not None:
              self.energy = self.energy - photon.energy
            if flouVacancy != False:
                # print(peCount, totalInteractionCount)
                self.secondartEffects(absorber, flouVacancy, pixelArray, collimatorCount, RSCount, eventSequence, depenergies)
                # print(peCount, totalInteractionCount)
            if photon != None and photon.energy > 1:  # unit Kev
                if type(absorber) is Absorber_flex.Detector_sub:
                    # if absorber == Absorber.CS or absorber == Absorber.I:
                    # absorber = Absorber.CsI
                    photon.currentAbsorber = absorber
                    # photon.currentAbsorber = scintillator
                # if absorber != air:
                #     ifInAir = False
                #     # photon.getCoefficients(ifInAir)
                # else:
                #     ifInAir = True
                # photon.getCoefficients(ifInAir)
                # fluoPhotonCount +=1
                # print(self.theta, self.phi)
                # print('Fluroscence Photon')
                # print(photon.theta, photon.phi)
                pathPts, absorbers, traLengths = [], [], []
                # print(photon.energy)
                return photon.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                                   eventSequence, depenergies)  # , density=absorber.density)
            # print(secondaryEffect)
            else:
                return collimatorCount, RSCount
        else:
            # print(peCount, totalInteractionCount)
            nonRadiTran = self.nonRad(absorber, shell)
            oriShell = nonRadiTran[0]
            if oriShell == 'K':
                if nonRadiTran[1:2] != 'X':
                    shell1 = nonRadiTran[1:3]
                else:
                    shell1 = 'X'
                if shell1 != 'X':
                    shell2 = nonRadiTran[3:]
                else:
                    shell2 = nonRadiTran[2]
            elif oriShell != 'M':
                if nonRadiTran[2] != 'X':
                    shell1 = nonRadiTran[2:4]
                else:
                    shell1 = 'X'
                shell2 = nonRadiTran[-1]
            else:
                shell1 = 'X'
                shell2 = 'X'
            if shell1 != 'X' and shell2 != 'X':
                # peC1, totalC1 = self.secondartEffects(absorber, shell1, pinhole, disc)
                # peC2, totalC2 = self.secondartEffects(absorber, shell2, pinhole, disc)
                # peCount = peC1+peC2
                # if shell1 == 'XX' or shell2 == 'XX':
                #     print(shell1, shell2)
                # print(totalInteractionCount)
                self.secondartEffects(absorber, shell1, pixelArray, collimatorCount, RSCount, eventSequence, depenergies)
                # print(totalInteractionCount)
                return self.secondartEffects(absorber, shell2, pixelArray, collimatorCount, RSCount, eventSequence, depenergies)
                # return tuple(map(operator.add, self.secondartEffects(absorber, shell1, pinhole, disc), self.secondartEffects(absorber, shell2, pinhole, disc)))
            elif shell1 != 'X':
                return self.secondartEffects(absorber, shell1, pixelArray, collimatorCount, RSCount, eventSequence, depenergies)
            elif shell2 != 'X':
                return self.secondartEffects(absorber, shell2, pixelArray, collimatorCount, RSCount, eventSequence, depenergies)
            # elif shell1 == 'XX' or shell2 == 'XX':
            #     print(shell1, shell2)
            else:
                # peCount += 1
                # totalInteractionCount += 1

                # if self.currentAbsorber == scintillator:
                #     print(shell1, shell2, self.energy, 'PE NonRa Deposit', self.x, self.y, self.z)

                return collimatorCount, RSCount  # peCount, totalInteractionCount

    def samplingAngles(self):
        phi = random.uniform(0, 2) * np.pi
        r = random.random()
        theta = np.arccos(1 - 2 * r)  # theta is a full range for a sphere because 2r, r for hemisphere
        # print(phi, theta)
        return theta, phi

    def csInteractTest(self, pixelArray, collimatorCount, RSCount, eventSequence):  # w, wDisc
        # global peCount, totalInteractionCount
        # print(self.phi, self.theta)
        # angleRad = 290
        # scattered_rad_pairs = []
        # while angleRad <=290:
        psiOri = self.csAngleDist()  # angleRad/180*np.pi#self.csAngleDist() # require a coordinate system rotation to original photon direction, angle between y' and x' axes.
        # print('CS', psiOri, self.energy, absorber)
        # if psiOri > np.pi :
        #     psi = 2*np.pi-psiOri
        # if psiOri <0:
        #     psi = -psiOri

        r = np.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

        # just for validating
        # rx =  r* np.sin(self.theta) * np.cos(self.phi)
        # ry =  r* np.sin(self.theta) * np.sin(self.phi)
        # rz =  r* np.cos(self.theta)
        # currentDirection = np.matrix([[rx], [ry], [rz]])

        rotation_phi = -self.phi
        rotation_z_vector = Rz(rotation_phi)

        rotation_theta = -self.theta
        rotation_y_vector2 = Ry(rotation_theta)

        rotation = rotation_y_vector2 * rotation_z_vector

        # rotated_vector = np.round(rotation*currentDirection, 6) # for rotation test

        # for new azimuth angle phi in the new coordinates
        # scattered_rad = 0
        # scattered_rad_pairs = []
        # while scattered_rad <=360:
        scattered_phi = random.uniform(0, 2) * np.pi  # -45/180*np.pi#random.uniform(0, 2) * np.pi
        rrx = r * np.sin(psiOri) * np.cos(scattered_phi)
        rry = r * np.sin(psiOri) * np.sin(scattered_phi)
        rrz = r * np.cos(psiOri)

        scatteredDirection = np.matrix([[rrx], [rry], [rrz]])

        rotation_inv = rotation.I
        # scattered_vector_test = np.round(rotation_inv * rotated_vector, 6)
        scattered_vector = np.round(rotation_inv * scatteredDirection, 6)
        sx, sy, sz = scattered_vector[0][0], scattered_vector[1][0], scattered_vector[2][0]
        # scattered_vector = np.matrix([[sx], [sy], [sz]])
        # scattered_angle_phi = np.arctan(sy/sx)
        # scattered_angle_theta = np.arctan(np.sqrt(sx**2+sy**2)/sz)
        scattered_angle_phi = np.arctan2(sy, sx)
        scattered_angle_theta = np.arctan2(np.sqrt(sx ** 2 + sy ** 2), sz)

        # m0c2 = (9.1e-31) * (3.0e+8) ** 2 #Joules -> 511 Kev
        energy = self.energy / (1 + self.energy / 511 * (1 - np.cos(psiOri)))
        scatteredPhoton = Photon(self.x, self.y, self.z, energy)
        scatteredPhoton.setDirection(scattered_angle_theta, scattered_angle_phi)
        # scatteredPhoton.currentAbsorber = absorber
        # print('CS', absorber, self.x, self.y, self.z)
        scatteredPhoton.Cscattered = True
        # print('CS')
        # print(scattered_angle_theta, scattered_angle_phi)
        pathPts, absorbers, traLengths = [], [], []
        return scatteredPhoton.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                                    eventSequence)  #

    def csInteract(self, traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                   eventSequence, depenergies, absorber):  # w, wDisc
        # global peCount, totalInteractionCount
        # print(self.phi, self.theta)
        # angleRad = 290
        # scattered_rad_pairs = []
        # while angleRad <=290:
        psiOri = self.csAngleDist()  # angleRad/180*np.pi#self.csAngleDist() # require a coordinate system rotation to original photon direction, angle between y' and x' axes.
        # print('CS', psiOri, self.energy, absorber)
        # if psiOri > np.pi :
        #     psi = 2*np.pi-psiOri
        # if psiOri <0:
        #     psi = -psiOri

        r = np.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

        # just for validating
        # rx =  r* np.sin(self.theta) * np.cos(self.phi)
        # ry =  r* np.sin(self.theta) * np.sin(self.phi)
        # rz =  r* np.cos(self.theta)
        # currentDirection = np.matrix([[rx], [ry], [rz]])

        rotation_phi = -self.phi
        rotation_z_vector = Rz(rotation_phi)

        rotation_theta = -self.theta
        rotation_y_vector2 = Ry(rotation_theta)

        rotation = rotation_y_vector2 * rotation_z_vector

        # rotated_vector = np.round(rotation*currentDirection, 6) # for rotation test

        # for new azimuth angle phi in the new coordinates
        # scattered_rad = 0
        # scattered_rad_pairs = []
        # while scattered_rad <=360:
        scattered_phi = random.uniform(0, 2) * np.pi  # -45/180*np.pi#random.uniform(0, 2) * np.pi
        rrx = r * np.sin(psiOri) * np.cos(scattered_phi)
        rry = r * np.sin(psiOri) * np.sin(scattered_phi)
        rrz = r * np.cos(psiOri)

        scatteredDirection = np.matrix([[rrx], [rry], [rrz]])

        rotation_inv = rotation.I
        # scattered_vector_test = np.round(rotation_inv * rotated_vector, 6)
        scattered_vector = np.round(rotation_inv * scatteredDirection, 6)
        sx, sy, sz = scattered_vector[0][0], scattered_vector[1][0], scattered_vector[2][0]
        # scattered_vector = np.matrix([[sx], [sy], [sz]])
        # scattered_angle_phi = np.arctan(sy/sx)
        # scattered_angle_theta = np.arctan(np.sqrt(sx**2+sy**2)/sz)
        scattered_angle_phi = np.arctan2(sy, sx)
        scattered_angle_theta = np.arctan2(np.sqrt(sx ** 2 + sy ** 2), sz)

        # m0c2 = (9.1e-31) * (3.0e+8) ** 2 #Joules -> 511 Kev
        energy = self.energy / (1 + self.energy / 511 * (1 - np.cos(psiOri)))
        if type(absorber) == Absorber_flex.Detector or issubclass(type(absorber), Absorber_flex.Detector):
            depenergies.append(self.energy - energy)
        # scatteredPhoton = Photon(self.x, self.y, self.z, energy)
        self.setDirection(scattered_angle_theta, scattered_angle_phi)

        # if self.currentAbsorber == scintillator:
        #     print(self.energy - energy,'CS Deposit', self.x, self.y, self.z)

        self.energy = energy
        self.Cscattered = True
        self.initialCheck = False
        # pathPts, absorbers, traLengths = [], [], []
        return self.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount, eventSequence, depenergies)  #

    def rsInteract(self, traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                   eventSequence, depenergies):  # w, wDisc
        # global peCount, totalInteractionCount
        # print(self.phi, self.theta)
        # angleRad = 290
        # scattered_rad_pairs = []
        # while angleRad <=290:
        absorber = self.currentAbsorber
        if self.currentAbsorber == scintillator:
            if self.currentAbsorber.energies is None or self.energy < 1:
                CSpeCo, CScsCo, CSrsCo, CStotalCo, IpeCo, IcsCo, IrsCo, ItotalCo = self.getCoefficients(
                    self.currentAbsorber)
            else:
                CSpeCo, CScsCo, CSrsCo, CStotalCo, IpeCo, IcsCo, IrsCo, ItotalCo = self.getCoefficients_file(
                    self.currentAbsorber)
            Pcs = (detector_sub2.mass * detector.atom2_weight / scintillator.mass) * (CStotalCo) / self.totalCo
            r = random.uniform(0, 1)
            if r <= Pcs:
                absorber = detector_sub2
                if absorber.energies is None or self.energy < 1:
                    self.getCoefficients(absorber)
                else:
                    self.getCoefficients_file(absorber)
            else:
                absorber = detector_sub1
                if absorber.energies is None or self.energy < 1:
                    self.getCoefficients(absorber)
                else:
                    self.getCoefficients_file(absorber)
        if absorber.ffs is not None:
            psiOri = self.rsAngleDist_file(absorber)
        else:
            psiOri = self.rsAngleDist(
                absorber)  # angleRad/180*np.pi#self.csAngleDist() # require a coordinate system rotation to original photon direction, angle between y' and x' axes.
        # print('RS', psiOri, self.energy, absorber)
        # if psiOri > np.pi :
        #     psi = 2*np.pi-psiOri
        # if psiOri <0:
        #     psi = -psiOri

        r = np.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

        # just for validating
        # rx =  r* np.sin(self.theta) * np.cos(self.phi)
        # ry =  r* np.sin(self.theta) * np.sin(self.phi)
        # rz =  r* np.cos(self.theta)
        # currentDirection = np.matrix([[rx], [ry], [rz]])

        rotation_phi = -self.phi
        rotation_z_vector = Rz(rotation_phi)

        rotation_theta = -self.theta
        rotation_y_vector2 = Ry(rotation_theta)

        rotation = rotation_y_vector2 * rotation_z_vector

        # rotated_vector = np.round(rotation*currentDirection, 6) # for rotation test

        # for new azimuth angle phi in the new coordinates
        # scattered_rad = 0
        # scattered_rad_pairs = []
        # while scattered_rad <=360:
        scattered_phi = random.uniform(0, 2) * np.pi  # -45/180*np.pi#random.uniform(0, 2) * np.pi

        rrx = r * np.sin(psiOri) * np.cos(scattered_phi)
        rry = r * np.sin(psiOri) * np.sin(scattered_phi)
        rrz = r * np.cos(psiOri)

        scatteredDirection = np.matrix([[rrx], [rry], [rrz]])

        rotation_inv = rotation.I
        # scattered_vector_test = np.round(rotation_inv * rotated_vector, 6)
        scattered_vector = np.round(rotation_inv * scatteredDirection, 6)
        sx, sy, sz = scattered_vector[0][0], scattered_vector[1][0], scattered_vector[2][0]
        # scattered_vector = np.matrix([[sx], [sy], [sz]])
        # scattered_angle_phi = np.arctan(sy/sx)
        # scattered_angle_theta = np.arctan(np.sqrt(sx**2+sy**2)/sz)
        scattered_angle_phi = np.arctan2(sy, sx)
        scattered_angle_theta = np.arctan2(np.sqrt(sx ** 2 + sy ** 2), sz)

        self.setDirection(scattered_angle_theta, scattered_angle_phi)
        # self.currentAbsorber = absorber
        # print('CS', absorber, self.x, self.y, self.z)
        self.Rscattered = True
        self.initialCheck = False
        # print('CS')
        # print(scattered_angle_theta, scattered_angle_phi)
        return self.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                         eventSequence, depenergies)  # , self.currentAbsorber.density)

    # def intersectWithScintillators(self, pt, cylinderDict, currentPosition, nextPosition,  pinhole, disc, density):
    #     # scintillatorsIntersectionDict = {}
    #     distanceToScintillator = {}
    #     pts =list(cylinderDict.keys())
    #     columnCentres = KDTree(pts)
    #     distance, index = columnCentres.query(pt, 1)
    #     columnKey = pts[index]
    #     column= cylinderDict.get(columnKey)
    #     # for cylinder in columnCentres:
    #     pts= column.intersectWithLine(currentPosition, nextPosition, tol=1e-04)
    #     if len(pts) != 0:
    #         pts = [tuple(e for e in pt) for pt in pts if tuple(e for e in pt) != currentPosition]
    #         for pt in pts:
    #             line = vedo.Line([pt, currentPosition])
    #             distanceToScintillator[pt] = line.length()
    #         # scintillatorsIntersectionDict[cylinder] = pts
    #     # for pts in scintillatorsIntersectionDict.values():
    #     #     # pt = tuple(e for e in pt)
    #     #     for pt in pts:
    #     #         line = vedo.Line([pt, currentPosition])
    #     #         distanceToScintillator[pt] = line.length()
    #     if len(distanceToScintillator.keys()) != 0:
    #         # newpoints = vedo.Points(list(distanceBetPts.keys()), r=10).color('red')
    #         # vedo.show(tcone, bcone, newpoints, axes=True)
    #         #     print(pinhole.r)
    #         # distanceBetPts = {k: v for k, v in sorted(distanceBetPts.items(), key=lambda  item:item[1])} # Python 3.7+ version, sort the dict by its values
    #         # possibleTravelLenLst = []
    #         pt = min(distanceToScintillator, key=distanceToScintillator.get)
    #         self.x = pt[0]
    #         self.y = pt[1]
    #         self.z = pt[2]
    #         self.materialChange = True
    #         if density != air.getDensity():
    #             density = air.getDensity()
    #             # ifInAir = True
    #         else:
    #             density = scintillator.getDensity()
    #             # ifInAir = False
    #             self.currentAbsorber = scintillator
    #         return self.main(traLengths, pinhole, disc,  density)

    # for pt in distanceBetPts.keys():

    def csAngleDist(self):
        r0 = 2.82e-13
        # m0c2 = (9.1e-31) * (3.0e+8) ** 2 #Joules
        m0c2 = 511  # Kev
        distls = []
        angles = []
        ka = self.energy / m0c2
        # psi = -180 #np.pi
        psiStart = 0
        stepSize = 0.1  # 0.1
        psiend = 180 + stepSize
        for psi in np.arange(psiStart, psiend, stepSize):
            # while (psi <= 180):
            psiRa = psi / 180 * np.pi
            c = ka * (1 - np.cos(psiRa))
            a = 1 + c
            b = 1 + (np.cos(psiRa)) ** 2.0
            first = (1 / a) ** 2
            second = b / 2
            third = 1 + c ** 2 / (b * a)
            d = first * second * third
            res = r0 ** 2 * d
            angles.append(psiRa)

            # if len(distls) >0:
            #     distls.append(np.sum(distls))
            # else:
            distls.append(res)
            # psi += 1/180*np.pi

            # psi += 0.1
        # max, min, sum = np.argmax(distls), np.argmin(distls), np.sum(distls)
        # distls = distls/np.max(distls)
        # ax = plt.subplot(111, projection='polar')
        # name = str(141) + " keV"
        # ax.plot(angles, distls, label=name)
        # plt.savefig("CSTest1.png", format="png", bbox_inches="tight")

        # max, min = np.max(distls), np.min(distls)
        # print(sum(distls))
        # plt.polar(angles, distls)
        # plt.show()
        # r= random.random()
        csAngle = Comparison.comparecsAngle(angles, distls)
        if csAngle == None:
            print('error')
        return csAngle

    # def ffInterpolat(self,x0, atomicFFs):
    #     for i in range(len(XS) - 1):
    #         if x0 == XS[i]:
    #             return atomicFFs[i]
    #         if x0 > XS[i] and x0 <= XS[i + 1]:
    #             return (atomicFFs[i] + (x0 - XS[i]) * (atomicFFs[i + 1] - atomicFFs[i]) / (XS[i + 1] - XS[i]))

    def ffInterpolat(self, x0, atomicFFs):
        for i in range(len(XS) - 1):
            if x0 == XS[i]:
                return atomicFFs[i]
            if x0 > XS[i] and x0 <= XS[i + 1]:
                return (atomicFFs[i] + (x0 - XS[i]) * (atomicFFs[i + 1] - atomicFFs[i]) / (XS[i + 1] - XS[i]))

    def rsAngleDist_file(self, absorber):
        r0 = (2.82e-15)
        # angles = np.arange(0, 180.1, 0.1) /180 *np.pi
        angles = np.arange(0, 180.1, 0.1) / 180 * np.pi
        a2 = r0 ** 2 * ((1 + np.cos(angles) ** 2) / 2)
        wavelen = 12398 / (self.energy * 1000)  # unit 0.1 nm
        q = np.sin(angles) / wavelen
        # f0 = absorber.offset
        # for s, e in zip(absorber.scale,
        #                 absorber.exponents):
        #     f0 += s * np.exp(-e * q * q)

        # sign = np.sign(q)
        q_abs = np.abs(q)
        # ffInterpolat_vc = np.vectorize(self.ffInterpolat)

        # if type(absorber) is Absorber.Cs:
        atomicFFs = absorber.ffs
        # elif type(absorber) is Absorber.I:
        #     atomicFFs = IATOMICFFs
        # elif type(absorber) is Absorber.W:
        #     atomicFFs = WATOMICFFs
        # elif type(absorber) is Absorber.Al:
        #     atomicFFs = AlATOMICFFs
        # else:
        #     print(absorber, 'extra absorber type needs to be dealt with')

        # f0 = ffInterpolat_vc(q_abs, atomicFFs)
        f0 = [self.ffInterpolat(q, atomicFFs) for q in q_abs]
        # ff0 = np.multiply(sign, f0)
        f0 = np.array(f0)
        dis = np.multiply(a2, f0 ** 2)  # / np.max(np.multiply(a2, f0))

        # dis = np.multiply(a2, f0) #/ np.max(np.multiply(a2, f0))
        # dis =dis/np.sum(dis)
        # r = random.random()
        rsAngle = Comparison.comparecsAngle(angles, dis)
        return rsAngle

    def rsAngleDist(self, absorber):
        r0 = (2.82e-15)
        # angles = np.arange(0, 180.1, 0.1) /180 *np.pi
        angles = np.arange(0, 180.1, 0.1) / 180 * np.pi
        a2 = r0 ** 2 * ((1 + np.cos(angles) ** 2) / 2)
        wavelen = 12398 / (self.energy * 1000)  # unit 0.1 nm
        q = np.sin(angles) / wavelen
        # f0 = absorber.offset
        # for s, e in zip(absorber.scale,
        #                 absorber.exponents):
        #     f0 += s * np.exp(-e * q * q)

        # sign = np.sign(q)
        q_abs = np.abs(q)
        # ffInterpolat_vc = np.vectorize(self.ffInterpolat)
        #
        # if type(absorber) is Absorber.Cs:
        #     atomicFFs = CsATOMICFFs
        # elif type(absorber) is Absorber.I:
        #     atomicFFs = IATOMICFFs
        # elif type(absorber) is Absorber.W:
        #     atomicFFs = WATOMICFFs
        # elif type(absorber) is Absorber.Al:
        #     atomicFFs = AlATOMICFFs
        # else:
        #     print(absorber, 'extra absorber type needs to be dealt with')

        # f0 = ffInterpolat_vc(q_abs, atomicFFs)
        # try:
        # f0 = ffInterpolat_vc(q_abs, atomicFFs)
        # try:
        #     for qq in q_abs:
        #         xraylib.FF_Rayl(absorber.atomNum, qq)
        # except Exception as e:
        #     print(e)
        f0 = [
            xraylib.FF_Rayl(absorber.atomNum, qq) if qq > 10 ** (-7) else xraylib.FF_Rayl(absorber.atomNum, 10 ** (-7))
            for qq in q_abs]
        # except Exception as e:
        #     print(e)

        # f0 = [self.ffInterpolat(q, atomicFFs) for q in q_abs]
        # ff0 = np.multiply(sign, f0)
        f0 = np.array(f0)
        dis = np.multiply(a2, f0 ** 2)  # / np.max(np.multiply(a2, f0))

        # dis = np.multiply(a2, f0) #/ np.max(np.multiply(a2, f0))
        # dis =dis/np.sum(dis)
        # r = random.random()
        rsAngle = Comparison.comparecsAngle(angles, dis)
        return rsAngle

    def sampleProbability(self):
        P = random.random()  # testing extreme case by setting random range original: random.random()
        return P

    def getTravellingLength(self, density, absorber):
        P = self.sampleProbability()
        # distance = np.log(1/(P-1))/(self.totalCo*W.density)
        # self.getCoefficients()
        # if self.ifPositionInsideAbsorber(W):
        if self.currentAbsorber != air and absorber.energies is None:
            self.getCoefficients(absorber)
        else:
            self.getCoefficients_file(absorber)
        # self.getCoefficients(absorber)
        # print(P)
        distance = np.log(1 / P) / (self.totalCo * density)  # cm^2/g *  g/cm^3 -> cm, P is probability passing through
        # distance = np.log(1 / P) / (self.totalCo * Air.density) # cm^2/g *  g/cm^3 -> cm
        return distance

    # def getTravellingLength(self):
    #     P = self.sampleProbability()
    #     # distance = np.log(1/(P-1))/(self.totalCo*W.density)
    #     # self.getCoefficients()
    #     if self.ifPositionInsideAbsorber(W):
    #         self.getCoefficients()
    #         distance = np.log(1 / P) / (self.totalCo * W.density)
    #     else:
    #         distance = np.log(1 / P) / (self.totalCo * Air.density)  # cm^2/g *  g/cm^3 -> cm
    #     # distance = np.log(1 / P) / (self.totalCo * Air.density) # cm^2/g *  g/cm^3 -> cm
    #     return distance
    # def getNextPostionDebug(self, absorber, debug):
    #     # dr = 0.01
    #     distance = self.getTravellingLength(absorber.density, absorber)
    #     # traLengths.append(distance)
    #     # self.theta = np.pi # fixed the direction of the generated photons
    #     # self.phi = 0
    #     dx = distance * np.sin(self.theta) * np.cos(self.phi)
    #     dy = distance * np.sin(self.theta) * np.sin(self.phi)
    #     dz = distance * np.cos(self.theta)
    #     if debug:
    #         self.x = self.x + dx
    #         self.y = self.y + dy
    #         self.z = self.z + dz
    #     # print(distance)
    #     # print(absorber)
    #     else:
    #         # for testing intersection points
    #         self.x = 0
    #         self.y = wDisc.r
    #         self.z = 0
    #         self.debug = True
    #
    #     return distance

    def getNextPostion(self, absorber):
        # dr = 0.01
        distance = self.getTravellingLength(absorber.density, absorber)
        # traLengths.append(distance)
        # self.theta = np.pi # fixed the direction of the generated photons
        # self.phi = 0
        dx = distance * np.sin(self.theta) * np.cos(self.phi)
        dy = distance * np.sin(self.theta) * np.sin(self.phi)
        dz = distance * np.cos(self.theta)
        # if debug:
        self.x = self.x + dx
        self.y = self.y + dy
        self.z = self.z + dz
        # print(distance)
        # print(absorber)
        # else:
        # # for testing intersection points
        #     self.x = 0
        #     self.y = wDisc.r
        #     self.z = wDisc.z
        #     self.debug = True

        return distance

    def insideCurrentAbsorberUpdated(self, absorber):  # check if the current photon is inside the material
        # if self.z >= W.z and self.z <= (W.z + W.t) and  np.abs(self.y) <= 150 and np.abs(self.x)<=150:
        # if self.z <= W.z and self.z >= (W.z - W.t) and np.abs(self.y) <= 150 and np.abs(self.x) <= 150:

        if self.z >= absorber.z and self.z <= (absorber.z + absorber.t):  # infinite long absorber with a thickness t
            if absorber.getSepration() == 0:
                return True
            else:
                sep = absorber.getSepration()
                d = np.sqrt(self.x * self.x + self.y * self.y)  #### UPDATE self.x and self.y #########
                if d > (sep / 2):
                    return True
                else:
                    return False

        else:
            return False

    def ifPositionInsideAbsorber(self, absorber):  # check if the current photon is inside the material
        # if self.z >= W.z and self.z <= (W.z + W.t) and  np.abs(self.y) <= 150 and np.abs(self.x)<=150:
        # if self.z <= W.z and self.z >= (W.z - W.t) and np.abs(self.y) <= 150 and np.abs(self.x) <= 150:

        if self.z >= absorber.z and self.z <= (absorber.z + absorber.t):  # infinite long absorber with a thickness t
            if absorber.getSepration() == 0:
                return True
            else:
                sep = absorber.getSepration()
                if np.abs(self.x) >= (sep / 2) or np.abs(self.y) >= (sep / 2):
                    return True
                else:
                    return False

        else:
            return False

    def passThroughAbsorber(self, absorber):
        if self.z > (absorber.z + absorber.t):
            return True
        else:
            return False

    # def outsideAbsorberUpdatedBelowOrOutXY(self):  # check if the current photon is inside the material
    #     # if self.z >= (W.z + W.t) or  np.abs(self.y) > 150 or np.abs(self.x) > 150:
    #     if self.z < w.z - w.t or np.abs(self.y) > 150 or np.abs(self.x) > 150:
    #         return True
    #     else:
    #         return False

    def checkTowards(self):
        if self.theta < np.pi and self.theta > 0:
            return True
        else:
            return False


def generateCylindricalSource(h, r, regionC):
    s = random.uniform(0, 1)
    theta = random.uniform(0, 2 * np.pi)
    regionX, regionY, regionZ = regionC
    regionTop = regionZ - h / 2
    regionBottom = regionZ + h / 2
    z = random.uniform(regionBottom, regionTop)
    r = np.sqrt(s) * r
    x = r * np.cos(theta) + regionX
    y = r * np.sin(theta) + regionY
    # z = z
    return x, y, z


def generateSphericalSource(radiusSphere, regionC):
    phi = random.uniform(0, 2) * np.pi
    ran = random.random()
    theta = np.arccos(1 - 2 * ran)  # theta is a full range for a sphere because 2r, r for hemisphere
    # print(phi, theta)
    # radiusRan = random.uniform()
    radiusRan = radiusSphere * np.cbrt(random.uniform(0, 1))
    xc, yc, zc = regionC
    x = radiusRan * np.sin(theta) * np.cos(phi)
    y = radiusRan * np.sin(theta) * np.sin(phi)
    z = radiusRan * np.cos(theta)
    x = x + xc
    y = y + yc
    z = z + zc

    return x, y, z


def visualisation(absorbers, scintillator):
    abs_list = []
    for ab in absorbers:
        if type(ab.shape) == vedo.shapes.Sphere or type(ab.shape) == vedo.shapes.Cone:
            if ab.r > 0:
                abs_list.append(ab.shape)
        elif ab.t > 0.0001:
            abs_list.append(ab.shape)
        if ab.hole is not None and ab.t > 0.0001:
            if np.array(ab.holeShape).size == 2:
                abs_list.append(ab.holeShape[0])
                abs_list.append(ab.holeShape[1])
            elif np.array(ab.holeShape).size == 1:
                if type(ab.holeShape) is vedo.shapes.Cylinder:
                    abs_list.append(
                        vedo.Cylinder(pos=(ab.x, ab.y, ab.z), r=ab.hole_side_dia,
                                      height=ab.t,
                                      axis=(0, 0, 1), res=120, cap=True, c='yellow'))
                else:
                    abs_list.append(ab.holeShape)
        # else:
        #     if ab.t > 0.001:
        #         abs_list.append(ab.shape)
    # ScintillatorBlock = scintillator.getScitillator

    abs_list.append(scintillator.getScitillator())
    abs_list.append(scintillator.getPlane())

    vedo.show(*abs_list, axes=True)
    print('vis')
    # sys.stdout.write('vis')
    # if np.array(ab.holeShape).size == 2:
    #     # tcone, bcone = w.holeShape
    #     vedo.show(ab.shape, ab.holeShape, window, scintillator.getScitillator, scintillator.getPlane(), axes=True)
    # elif np.array(w.holeShape).size == 1:
    #     # tcone, bcone = None, None
    #     if type(w.holeShape) is vedo.shapes.Cylinder:
    #         vedo.show(discCylinder,
    #                   vedo.Cylinder(pos=[(w.x, w.y, w.z), (w.x, w.y, w.z + w.t)], r=w.hole_side_dia, height=w.t,
    #                                 axis=(0, 0, 1), res=120, cap=True, c='yellow'), ScintillatorBlock, CCDPlane,
    #                   axes=True)
    #     else:
    #         vedo.show(discCylinder, w.holeShape, ScintillatorBlock, CCDPlane, axes=True)


def generatePlaneSource(regionL, regionW, regionC):
    regionX, regionY, regionZ = regionC
    regionL = np.abs(regionL)
    regionW = np.abs(regionW)
    regionFront = regionY - regionW / 2
    regionEnd = regionY + regionW / 2
    regionRight = regionX + regionL / 2
    regionLeft = regionX - regionL / 2
    regionSy = random.uniform(regionFront, regionEnd)
    regionSx = random.uniform(regionRight, regionLeft)
    return regionSx, regionSy, regionZ


# samplesNum =100
distances = []
threads = []
# energy = 141#141 #kEv 68 and 70 for testing W K-Edge
energy = source_energy
energies = []
# totalPeCounts = []
means = []
stds = []
meansP = []
stdsP = []
# ds = []
# d = 0 #max 20
# traLengths = []
# source to collimator distance: 1cm, 5cm, 10cm
# distanceSA = -5 #0 #cm -1, -5, -10cm
# distanceSA = args.distance
# pixelArray = Detector.getCCDArray()
# eventSequence = []
# peCounts = []
# detectorCounts = []
# outsideCounts = []
# airPECounts = []
startTime = time.time()
# for i in range (1): #9
# TOTALPNUM = 1e5
TOTALPNUM = source_ppNum


# DEG = args.angle
# DEG = 90

def concurrencyTest(count):
    pixelArray = detector.getCCDArray()
    eventSequence = []
    depenergies = []
    collimatorCount = 0
    # for j in range(1):
    # distanceSA =- 0.5
    # distances.append(np.abs(distanceSA))
    # energys.append(energy)
    # ds.append(d)
    # for j in range(samplesNum):
    # traLengths =[]
    initialCount = 0
    RSCount = 0
    # file = 'events_test_consistent' + '_' + str(distanceSA) + '.txt'
    # f = open(file, 'w')
    # start =time.time()
    while initialCount < source_ppNum:  # 100000
        # p = Photon(0, 0, d, 141)
        # print(initialCount)
        # if DEG == 90:
        #   sideLen = 0
        # else:
        #   sideLen = (np.abs(distanceSA)+w.t/2)/(np.tan(DEG/180*np.pi))
        #   sideLen = np.round(sideLen, 6)
        # if source_pos == 'x':
        source_x, source_y, source_z = source_pos
        if source_type == 'Po':
            p = Photon(source_x, source_y, source_z, source_energy)  # p = Photon(0, 0, 0, energy)
        elif source_type == 'Pl':
            regionSx, regionSy, regionSz = generatePlaneSource(source_length, source_width, source_pos)
            p = Photon(regionSx, regionSy, regionSz, source_energy)
        elif source_type == 'Sp':
            regionSx, regionSy, regionSz = generateSphericalSource(source_radius, source_pos)
            p = Photon(regionSx, regionSy, regionSz, source_energy)
        elif source_type == 'V_Cy':
            regionSx, regionSy, regionSz = generateCylindricalSource(source_height, source_radius, source_pos)
            p = Photon(regionSx, regionSy, regionSz, source_energy)
        elif source_type == 'H_Cy':
            regionSx, regionSy, regionSz = generateCylindricalSource(source_height, source_radius, source_pos)
            translatedCylinderPoint = np.matrix(
                [[regionSx - source_pos[0]], [regionSy - source_pos[1]], [regionSz - source_pos[2]]])
            regionSx, regionSy, regionSz = np.matrix([[source_pos[0]], [source_pos[1]], [source_pos[2]]]) + Rx(
                np.radians(90)) * translatedCylinderPoint
            regionSx, regionSy, regionSz = regionSx.item(0, 0), regionSy.item(0, 0), regionSz.item(0, 0)
            p = Photon(regionSx, regionSy, regionSz, source_energy)
        else:
            # source_type = 'Po'
            p = Photon(source_x, source_y, source_z, source_energy)  # by defalut is point source

        traLengths = []
        pathPts = []
        absorbers = []
        collimatorCount, RSCount = p.main(traLengths, pathPts, absorbers, pixelArray, collimatorCount, RSCount,
                                          eventSequence, depenergies)
        initialCount += 1
        # print(initialCount)
    # f.close()
    # peCounts.append(peCount)
    # detectorCounts.append(detectorCount)
    # outsideCounts.append(outsideCount)
    # airPECounts.append(airPECount)

    # print(totalInteractionCount)

    # p.start()
    # for t in threads:
    #     t.join()
    #     print(j)
    # plt.scatter(CCDPlotsX, CCDPlotsY)
    # plt.imshow(pixelArray)
    # plt.show()
    print((time.time() - startTime) / 60)
    print(np.abs(source_z), np.sum(pixelArray), len(eventSequence), len(depenergies))
    print(collimatorCount, RSCount)
    print('CPU: ' + str(count))


    # scoop.logger.log(20,(str(time.time() - startTime) / 60))
    # scoop.logger.log(20,("\r%d, %d, %d\n" % (np.abs(source_z), np.sum(pixelArray), len(eventSequence))))
    # scoop.logger.log(20,("\r%d, %d\n" %(collimatorCount, RSCount)))
    # scoop.logger.log(20,("CPU "+str(count)+'\n'))

    # sys.stdout.write(str((time.time() - startTime) / 60)+'\n')
    # sys.stdout.write("\r%d, %d, %d\n" % (np.abs(source_z), np.sum(pixelArray), len(eventSequence)))
    # sys.stdout.write("\r%d, %d\n" %(collimatorCount, RSCount))
    # sys.stdout.write("CPU "+str(count)+'\n')
    # print(np.mean(peCounts), np.mean(detectorCounts), np.mean(otherCount),totalInteractionCount, fluoPhotonCount, csCount, np.abs(distanceSA)) # np.mean(outsideCounts),  np.mean(airPECounts)
    # print(str((start - time.time()) / 60) + ('seconds'))

    # sys.stdout.flush()

    return initialCount, collimatorCount, pixelArray, RSCount, eventSequence, depenergies


# print(concurrencyTest(0)[:-1])
totalCount = 0
totalPixelArray = np.zeros(pixelArray.shape)
totalEventSequence = []
totalDepositEnergies = []
# pool = mp.Pool(mp.cpu_count())
totalCollimatorCount = 0
totalRsCount = 0
from datetime import datetime

# datetime object containing current date and time
now = datetime.now()

# dd/mm/YY H:M:S
dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")

if __name__ == "__main__":
    startTimeTotal = time.time()
    import faulthandler

    faulthandler.enable()
    # with MPIPoolExecutor(mp.cpu_count()) as executor:
    if secondary != True:
        fname = str(source_pos) + '_' + str(source_energy) + '_' + str(source_type) + '_' + str(
            TOTALPNUM)  + dt_string + '_NoSecondEffects.txt'
    else:
        fname = str(source_pos) + '_' + str(source_energy) + '_' + str(source_type) + '_' + str(
            TOTALPNUM) + dt_string + '_WithSecondEffects.txt'
    # foutput = open('K-Energy Testing/Console_'+fname, 'w')
    # with redirect_stdout(foutput):
    # concurrencyTest(1)
    # pool = mp.Pool(3)
    # pool = mp.Pool(mp.cpu_count())
    # result= [pool.map(concurrencyTest,list(range(mp.cpu_count())))]
    vedo.Plotter(offscreen=not vis)
    # concurrencyTest()
    # pool = mp.Pool(3)
    if vis:
        visualisation(def_absorbers, scintillator)
    # import multiprocessing as mp
    # pool = mp.Pool(process_num)
    # result= [pool.map(concurrencyTest,list(range(process_num)))][0]
    result = list(futures.map(concurrencyTest, list(range(process_num))))
    # # fs = [futures.submit(concurrencyTest, c) for c in list(range(16))]
    # # for f in fs:
    # #   re = f.result()
    #   # pool.close()
    # pool.close()
    # pool.join()
    # print(len(result[0]))
    for re in result:
        individualCount, individualCollimatorCount, individualPixelArray, RSCount, eventSequence, depenergies = re
        totalRsCount += RSCount
        totalCount += individualCount
        totalPixelArray += individualPixelArray
        totalCollimatorCount += individualCollimatorCount
        totalEventSequence += eventSequence
        totalDepositEnergies += depenergies
    print(param_f, source_pos, source_type, str(secondary), energy, totalCount,
          totalCollimatorCount, totalRsCount, np.max(totalPixelArray), np.sum(totalPixelArray), len(totalEventSequence), len(totalDepositEnergies))
    # foutput.close()
    with open('Readout/Readout_' + fname, 'wb') as totalPixelArrayFile:
      np.savetxt(totalPixelArrayFile, totalPixelArray, fmt='%i')
    with open('Sequences/Sequences_' + fname, 'wb') as totalEventSequenceFile:
      np.savetxt(totalEventSequenceFile, np.array(totalEventSequence), fmt='%i')
    with open('DepEnergies/DepEnergies_' + fname, 'wb') as totalDepositEnergyFile:
      np.savetxt(totalDepositEnergyFile, np.array(totalDepositEnergies), fmt='%s')
    print('totalTime:', (time.time() - startTimeTotal) / 60)
    #plt.hist(totalDepositEnergies)
    #plt.savefig('DepEnergies/energyHist_'+fname+'.png')
    # result = list(futures.map(concurrencyTest, list(range(process_num))))
    # # # fs = [futures.submit(concurrencyTest, c) for c in list(range(16))]
    # # # for f in fs:
    # # #   re = f.result()
    # #   # pool.close()
    # # pool.close()
    # # pool.join()
    # # print(len(result[0]))
    # for re in result:
    #     individualCount, individualCollimatorCount, individualPixelArray, RSCount, eventSequence = re
    #     totalRsCount += RSCount
    #     totalCount += individualCount
    #     totalPixelArray += individualPixelArray
    #     totalCollimatorCount += individualCollimatorCount
    #     totalEventSequence += eventSequence
    # print(param_f, source_pos, source_type, str(secondary), str(w.getSepration()), energy, totalCount,
    #       totalCollimatorCount, totalRsCount, np.max(totalPixelArray), np.sum(totalPixelArray), len(totalEventSequence))
    # # foutput.close()
    # with open('Readout/Readout_' + fname, 'wb') as totalPixelArrayFile:
    #   np.savetxt(totalPixelArrayFile, totalPixelArray, fmt='%i')
    # with open('Sequences/Sequences_' + fname, 'wb') as totalEventSequenceFile:
    #   np.savetxt(totalEventSequenceFile, np.array(totalEventSequence), fmt='%i')
    # print('totalTime:', (time.time() - startTimeTotal) / 60)
#     np.savetxt('K-Energy Testing/Readout_'+fname, totalPixelArray)
# else:
#     np.savetxt('Readout_' + str(energy) + '_' + str(abs(distanceSA) + 0.3) + '_' + str(TOTALPNUM) + '_' + str(
#         DEG) + 'deg_WithSecondEffects.txt', totalPixelArray)
# plt.imshow(totalPixelArray)
# plt.show()
#########Uncooment for different energies plots

#         print(i)
#         peCounts.append(peCount)
#         # energy += 100
#         # d -= 0.5
#         # print(i)
#     peCounts = np.array(peCounts)
#
#     # peCountsP = peCounts/samplesNum
#     # yerr = np.sqrt(samplesNum)
#     # peCountsP = [x / totalInteractionCount for x in peCounts]
#     peCountsP = peCounts/totalInteractionCount
#     std = np.std(peCounts)
#     mean =np.mean(peCounts)
#     stdP = np.std(peCountsP)
#     meanP =np.mean(peCountsP)
#     stds.append(std)
#     means.append(mean)
#     stdsP.append(stdP)
#     meansP.append(meanP)
#     # energy += 100
#     # distanceSA -= 0.2
#     # energies.append(energy)
#     # yerr = np.sqrt(totalInteractionCount)
# fig, axs = plt.subplots(1, 2)
# # axs[0].errorbar(energies, means, yerr=stds, fmt='o')
# # axs[0].set_xlabel('photon energy / KeV')
#
# axs[0].errorbar(distances, means, yerr=stds, fmt='o')
# axs[0].set_xlabel('photon distance / cm for 141 kEv')
# # plt.plot(ds, peCounts)
# # plt.xlabel('distance of photon to matter')
# axs[0].set_ylabel('absorption photon counts for '+str(totalInteractionCount)+" photons")
# # axs[1].scatter(energys, peCountsP)
# # axs[1].set_xlabel('photon energy  / KeV')
# # axs[1].set_ylabel("The portion of absorption")
# axs[1].errorbar(energies, meansP, yerr=stdsP, fmt='s')
# axs[1].set_xlabel('photon energy  / KeV')
# # plt.plot(ds, peCounts)
# # plt.xlabel('distance of photon to matter')
# axs[1].set_ylabel('The portion of absorption')
# fig.tight_layout()
# # fig2 = plt.figure(2)
# # plt.plot(traLengths)
#
# # plt.hist(traLengths, bins = len(traLengths))
# # plt.xlabel('cm')
# # plt.ylabel('count')
#
#
#
# plt.show()

###################
