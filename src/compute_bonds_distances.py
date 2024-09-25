# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 21:09:15 2022

@author: jdkem
"""

####################
# Import Libraries #
####################
import math
 
# Function to find distance
def compute_distance(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return math.sqrt(dx*dx + dy*dy + dz*dz)
    
# Function to update bonds and bond coeffs based on hybridization
def reset_bondtypes_based_on_hybridization(mm):
    # Find unique bond types
    unique_bondtypes = set([]) # tuple(element1, element2)
    for i in mm.bonds:
        id1, id2 = mm.bonds[i].atomids
        element1 = mm.atoms[id1].element
        element2 = mm.atoms[id2].element
        hybrid1 = mm.atoms[id1].hybridization
        hybrid2 = mm.atoms[id2].hybridization
        atom1 = '{}_{}'.format(hybrid1, element1)
        atom2 = '{}_{}'.format(hybrid2, element2)
        bondtype = tuple(sorted([atom1, atom2]))
        unique_bondtypes.add(bondtype)
    
    # Set bond type ID
    unique_bondtypes = sorted(unique_bondtypes)
    bondtypes_forward = {} # { bondID : tuple(element1, element2) }
    bondtypes_reverse = {} # { tuple(element1, element2) : bondID }
    for n, i in enumerate(unique_bondtypes):
        bondtypes_forward[n+1] = i
        bondtypes_reverse[i] = n+1
    
    # Add bonds to m.bonds and update m.nbonds and m.nbondtypes
    new_bonds = {}
    class Bonds: pass # .type .atomids .symbols
    for i in mm.bonds:
        bond = mm.bonds[i]
        id1, id2 = bond.atomids
        element1 = mm.atoms[id1].element
        element2 = mm.atoms[id2].element
        hybrid1 = mm.atoms[id1].hybridization
        hybrid2 = mm.atoms[id2].hybridization
        atom1 = '{}_{}'.format(hybrid1, element1)
        atom2 = '{}_{}'.format(hybrid2, element2)
        bondtype = tuple(sorted([atom1, atom2]))
        b = Bonds()
        b.type = bondtypes_reverse[bondtype]
        b.atomids = bond.atomids
        b.symbols = bondtype
        new_bonds[i] = b

    # Update m.bond_coeffs
    new_bond_coeffs = {}
    class Coeff_class: pass  # .type .coeffs .symbols
    for i in bondtypes_forward:
        C = Coeff_class()
        C.type = i
        C.coeffs = []
        C.symbols = bondtypes_forward[i]
        new_bond_coeffs[i] = C
    return new_bonds, new_bond_coeffs

################################################################
# Statistic's functions to use for analyzing bond stats. *NOTE #
# not using numpy to make this code have zero dependancies*    #
################################################################
def compute_mean(data):
  return sum(data)/len(data)
 
def compute_variance(data):
  mean = compute_mean(data)
  deviations = [(x - mean)**2 for x in data]
  variance = sum(deviations)/len(data)
  return variance
 
def compute_standard_deviation(data):
  variance = compute_variance(data)
  return math.sqrt(variance)

####################################################
# Class to compute bond distances and analyze them #
####################################################
class Types: pass # .symbols .count .avg .min .max .std
class compute:
    def __init__(self, m, compute_bond_dists, struct):
        self.bonds = {} # { bondid : bond distance }
        self.types = {} # { bond type id : type object }
        
        # Reset bond types based on hybridization if the structure is kept and analyzed
        if struct == 'kept' and compute_bond_dists['hybridization']:
            bonds, bond_coeffs = reset_bondtypes_based_on_hybridization(m)
        else: bonds = m.bonds; bond_coeffs = m.bond_coeffs;
        
        # Find box dimensions to remove periodic boundary conditions
        xline = m.xbox_line.split()
        yline = m.ybox_line.split()
        zline = m.zbox_line.split()
        lx = float(xline[1])-float(xline[0])
        ly = float(yline[1])-float(yline[0])
        lz = float(zline[1])-float(zline[0])
        
        # create array to store bond types in
        bond_type_holder = {} # { bond type id : [list of all lengths] }
        for i in bond_coeffs:
            bond_type_holder[i] = []
            
        # Initialize self.types with zeros
        for i in bond_coeffs:
            t = Types()
            t.symbols = bond_coeffs[i].symbols
            t.count = 0
            t.avg = 0
            t.min = 0
            t.max = 0
            t.std = 0
            self.types[i] = t
       
        
        # Loop through bonds
        for i in bonds:
            bond = bonds[i]
            id1, id2 = bond.atomids
            
            # Find id1 and id2 x, y, and z
            id1x = m.atoms[id1].x
            id1y = m.atoms[id1].y
            id1z = m.atoms[id1].z
            
            id2x = m.atoms[id2].x
            id2y = m.atoms[id2].y
            id2z = m.atoms[id2].z
            
            # Shift atoms using minimum image convention if bond
            # is periodic (always shifting id2 atom if needed)
            bondlx = abs(id1x - id2x)
            if bondlx > 0.5*lx:
                diff = id1x - id2x
                if diff < 0:
                    id2x = id2x - lx
                elif diff > 0:
                    id2x = id2x + lx
                else:
                    id2x = id2x
                    
            bondly = abs(id1y - id2y)
            if bondly > 0.5*ly:
                diff = id1y - id2y
                if diff < 0:
                    id2y = id2y - ly
                elif diff > 0:
                    id2y = id2y + ly
                else:
                    id2y = id2y

            bondlz = abs(id1z - id2z)
            if bondlz > 0.5*lz:
                diff = id1z - id2z
                if diff < 0:
                    id2z = id2z - lz
                elif diff > 0:
                    id2z = id2z + lz
                else:
                    id2z = id2z

                
            # compute bond distance
            dist = compute_distance(id1x, id1y, id1z, id2x, id2y, id2z)
            
            # Add bond distance to self.bonds and append to holder based on type-id
            self.bonds[i] = dist
            bond_type_holder[bond.type].append(dist)
            
        ###########################
        # Find bonding statistics #
        ###########################
        for i in bond_type_holder:
            lst = bond_type_holder[i]            
            self.types[i].count = len(lst)
            self.types[i].avg = '{:.4f}'.format( compute_mean(lst) )
            self.types[i].min = '{:.4f}'.format( min(lst) )
            self.types[i].max = '{:.4f}'.format( max(lst) )
            self.types[i].std = '{:.4f}'.format( compute_standard_deviation(lst) )
