# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 16th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""


##############################
# Import Necessary Libraries #
##############################
from collections import OrderedDict
from tqdm import tqdm
import operator
import math
import sys

# Try importing math.dist and set flag. not all python versions have math.dist, but
# math.dist is quicker then native distance functions. So default is to try using
# math.dist and the user can switch out for compute_distance if needed with a performance
# trade off.
try:
    test = math.dist([0, 0, 0], [1, 1, 1])
except: 
    print('WARNING currently installed python math module does not contatin math.dist.')
    print('To use bonds_via_distanc.py please search through file and replace: ')
    print('    distance_fff = math.dist( [x1, y1, z1], [x2, y2, z2] ) w/')
    print('    distance_fff = compute_distance(x1, y1, z1, x2, y2, z2)     about line: 235')
    print('                           and')
    print('   distance_ppp = math.dist( [x1, y1, z1], [x2i, y2i, z2i] w/')
    print('   distance_ppp = compute_distance(x1, y1, z1, x2i, y2i, z2i)   about line: 255')
    print('Commenting and uncommenting will be the best method. math.dist is much quicker')
    print('the a native python distance function.')





##################################
# vdw radius to search for bonds #
##################################
# Set elemental vdw radius for bond search (User can adjust this as needed)
# Ref1 = Batsanov, Stepan S. "Van der Waals radii of elements." Inorganic materials 37.9 (2001): 871-885.
# Ref2 = https://en.wikipedia.org/wiki/Van_der_Waals_radius
# The goal is to set via Ref1, but if not found in Ref1 use wikipedia Ref2
vdw_radius = {'C':  1.70, # Ref1
              'H':  1.20, # Ref1
              'O':  1.55, # Ref1
              'N':  1.60, # Ref1
              'S':  1.80, # Ref1
              'F':  1.50, # Ref1
              'Si': 2.10, # Ref1
              'Xe': 2.16, # Ref2
              'Ne': 1.54, # Ref2
              'Kr': 2.02, # Ref2
              'He': 1.40, # Ref2
              'D':  2.40, # Ref1 (2*H vdw radius)
              'Cl': 1.80, # Ref1
              'Ca': 2.40, # Ref1
              'Br': 1.90, # Ref1
              'Ar': 1.88, # Ref2
              'P':  1.90, # Ref1
              'Al': 2.10, # Ref1
              'Mg': 2.20, # Ref1
              'Li': 2.20, # Ref1
              'Fe': 2.05, # Ref1
              'Na': 2.40, # Ref1
              'K':  2.80, # Ref1
              'Cs': 3.00, # Ref1
              'Ba': 2.70, # Ref1
              'Sr': 2.55, # Ref1
              'Pb': 2.30, # Ref1
              }

# Function to get vdw radii
def get_vdw_radii(element, vdw_radius):
    if element in vdw_radius:
        return vdw_radius[element]
    else: 
        print(f'\nERROR {element} element not in vdw_radius dictionary')
        print(f'inside bonds_via_distance.py PLEASE ADD {element}')
        sys.exit()
        
        
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
    
    
##############################################################
# Class to generate bonds via interatomic distance searching #
##############################################################
class Stats:
    pass # .count .avg .min .max .std .cutoff
    
class Element_nb:
    pass # 
    
class generate:
    def __init__(self, m, boundary, vdw_radius_scale, maxbonded):
        self.bonds = [] # lst of viable bonds
        self.flagged_bonds = [] # lst of flagged bonds
        self.statistics = {} # { bondtype : Stats object}
        self.nb_stats = {} # { element : Element_nb object }
        self.images = [] # lst of images searched
        self.nb_count = {} # { element : nb_dict }
    
    
        ################################################
        # Generate images via minimum image convention #
        ################################################
        # Boundary conditions
        pflags = boundary.split() # split boundary
        count = pflags.count('f') + pflags.count('p')
        if len(pflags) != 3 and count != 3:
            print('ERROR boundary does not contain 3-dimensions or boundary flags are not f or p'); sys.exit()
        if 'p' in pflags and m.xy != 0 or m.xz != 0 or m.yz != 0:
            print('ERROR simulaiton cell is NOT orthogonal and a periodic boundary is trying')
            print('to be applied. Currently PBCs are limited to orthogonal boxes only ...'); sys.exit()
        
        # Image information
        nimages = 1 # only use minimum image convention
        images = set([]) # set to hold unique images
        # Loop through nimages to generate image flags
        for ix in range(-nimages, nimages+1):
            for iy in range(-nimages, nimages+1):
                for iz in range(-nimages, nimages+1):
                    
                    # Update image based on boundary conditions
                    if ix != 0 and pflags[0] == 'f': ix = 0
                    if iy != 0 and pflags[1] == 'f': iy = 0
                    if iz != 0 and pflags[2] == 'f': iz = 0
                
                    # Log image
                    images.add( (ix, iy, iz) )
                    
        # Sort images in ascending order of magnitude to reduce run time when 
        # finding bond distances, since most bonds will be found at (0, 0, 0)
        self.images = sorted(images, key=lambda x: sum([abs(i) for i in x]) )
    
    
        ####################################
        # Find set simulataion cell bounds #
        ####################################
        xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
        xlo = float(xline[0]); xhi = float(xline[1]);
        ylo = float(yline[0]); yhi = float(yline[1]);
        zlo = float(zline[0]); zhi = float(zline[1]);
        lx = xhi-xlo; ly = yhi-ylo; lz = zhi-zlo;
        
        # Build scaled images for quicker calculations later on (Li directions already multiplied)
        scaled_images = [(ix*lx, iy*ly, iz*lz) for (ix, iy, iz) in self.images]            

        
        ##################################################################################
        # Function to check if the atom is near a box edge. Then build near edge dict to #
        # use as a look up table to aviod multiple re-calculations as the loops progress #
        ##################################################################################
        def check_near_edge(x, y, z, max_from_edge, xlo, xhi, ylo, yhi, zlo, zhi):
            if abs(x - xlo) < max_from_edge or abs(x + xlo) < max_from_edge: pbcflag = True
            elif abs(x - xhi) < max_from_edge or abs(x + xhi) < max_from_edge: pbcflag = True        
            elif abs(y - ylo) < max_from_edge or abs(y + ylo) < max_from_edge: pbcflag = True
            elif abs(y - yhi) < max_from_edge or abs(y + yhi) < max_from_edge: pbcflag = True
            elif abs(z - zlo) < max_from_edge or abs(z + zlo) < max_from_edge: pbcflag = True
            elif abs(z - zhi) < max_from_edge or abs(z + zhi) < max_from_edge: pbcflag = True
            else: pbcflag = False
            return pbcflag
        
        # Determine "nearness" to edge. Only search periodic bonds if atom is this close to an edge.
        # Only search largest possible bond found in the system and nothing more to manage run time
        elements = sorted({m.atoms[i].element for i in m.atoms})
        max_from_edge = vdw_radius_scale*max([get_vdw_radii(element, vdw_radius) for element in elements]) 
        edgeflags = {} # { atomID : edge flag }
        for i in m.atoms:
            x = m.atoms[i].x; y = m.atoms[i].y; z = m.atoms[i].z;
            edgeflags[i] = check_near_edge(x, y, z, max_from_edge, xlo, xhi, ylo, yhi, zlo, zhi)


        ###############################################################################################
        # Generate lookup tables for distance cutoffs based on bondtypes for quicker runtime later on #
        ###############################################################################################
        maxdist_dict = {} # { tuple(element1, element2) : cutoff distance }
        for element1 in elements:
            radius1 = vdw_radius_scale*get_vdw_radii(element1, vdw_radius)
            for element2 in elements:
                radius2 = vdw_radius_scale*get_vdw_radii(element2, vdw_radius)
                
                # Compute distance and log in forward and reverse ordering
                maxdistance = 0.5*radius1 + 0.5*radius2
                maxdist_dict[(element1, element2)] = maxdistance
                maxdist_dict[(element2, element1)] = maxdistance
                
        
        ###################################################################################
        # Function to compute distance when math.dist is not available. math.dist will be #
        # quicker but for some reason not all python has math.dist .... compute_distance  #
        # will be time optimized but still will be slower then math.dist.                 #
        ###################################################################################        
        def compute_distance(x1, y1, z1, x2, y2, z2):
            dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
            return math.sqrt(dx*dx + dy*dy + dz*dz)
        
        
        ####################################
        # Find possible bonds via distance #
        ####################################
        possible_bonds = {} # { tuple(id1, id2): distance_fff or distance_ppp }
        self.bond_status = {'periodic': 0, 'non-periodic': 0} # to tally bond types
        id2_atoms = [i for i in m.atoms] # atoms to loop over in inner loop (will be reduced each outer loop iteration)
        
        # Start finding interatomic distances
        print('\n\nFinding interatomic distances for bond creation ....')
        for id1 in tqdm(m.atoms):      
            # Find atom1 to access x1, y1, and z1
            atom1 = m.atoms[id1];
            
            # Find x1, y1, and z1 for interatomic searching 
            x1 = atom1.x; y1 = atom1.y; z1 = atom1.z;
            
            # remove id1 from id2_atoms since it will be searched from outer most loop
            id2_atoms.remove(id1)
            for id2 in id2_atoms:  
                # Skip over if outer most loop and innner most loop have the same atomID
                if id1 == id2: continue        
                      
                # Find atom2 and max distance cutoff
                atom2 = m.atoms[id2]; maxdist = maxdist_dict[(atom1.element, atom2.element)];
                 
                # Find x2, y2, and z2 for interatomic searching
                x2 = atom2.x; y2 = atom2.y; z2 = atom2.z;

                # Compute non-periodic distance. math.dist is quicker then compute_distance function, 
                # however for some reason on MTU's HPC the math library does not have a .dist class
                # so comment/uncomment as needed ....
                distance_fff = math.dist( [x1, y1, z1], [x2, y2, z2] )
                #distance_fff = compute_distance(x1, y1, z1, x2, y2, z2)
                    
                
                # If non-periodic distance is acceptable log
                if distance_fff <= maxdist:
                    possible_bonds[tuple(sorted([id1, id2]))] = distance_fff
                    self.bond_status['non-periodic'] += 1
                
                
                # elif both id1 and id2 are found to be a maximum of the largest possible bond distance
                # away from any of the simulaiton cell walls, attempt finding the periodic distance by
                # shifting id2 throughout all possible images that have already been scaled by box dims
                elif edgeflags[id1] and edgeflags[id2]:
                    for (ixlx, iyly, izlz) in scaled_images:
                        x2i = x2+ixlx; y2i = y2+iyly; z2i = z2+izlz;
                        
                        # Compute distance. math.dist is quicker then compute_distance function, however
                        # for some reason on MTU's HPC the math library does not have a .dist class so
                        # comment/uncomment as needed ....
                        distance_ppp = math.dist( [x1, y1, z1], [x2i, y2i, z2i] )
                        #distance_ppp = compute_distance(x1, y1, z1, x2i, y2i, z2i)
                        
                        # If periodic distance is acceptable, log and break the images loop
                        if distance_ppp <= maxdist:                          
                            possible_bonds[tuple(sorted([id1, id2]))] = distance_ppp
                            self.bond_status['periodic'] += 1
                            break


                    
        ####################################
        # Find bonded atoms to be able to  #
        # set max nb cut-off per atom type #
        ####################################
        # Intialize atoms dictionary
        bonded = {i:[] for i in m.atoms}; # { atomid : [lst of bonded atoms]}
       
        # add in connect atoms to bonded dictionary
        for id1, id2 in possible_bonds:     
            bonded[id1].append(id2)
            bonded[id2].append(id1)
        # ordering dictionary    
        bonded = dict(OrderedDict(sorted(bonded.items())))

       
        ########################################
        # Start finding and flagging bonds     #
        # that meet cut-off values set by user #
        ########################################
        bonds = set([]); flagged_bonds = set([]);
        for id1 in bonded:
            connected = bonded[id1]
           
            # Creating dictionary to sort in decending order
            dist_dict = {}
            for id2 in connected:
                bond = tuple( sorted([id1, id2]) ) 
                dist_dict[bond] = possible_bonds[bond]
               
            # Ordering dictionary by value (bond length) in ascending order
            dist_dict = dict(sorted(dist_dict.items(), key=operator.itemgetter(1)))
           
            # Find max_nb based on element
            max_nb = maxbonded[m.atoms[id1].element]
           
            # Creating bonds based on inputs values
            count_nb = 0;
            for bond in dist_dict:
                count_nb += 1
                # If number of bonds is less than specified
                if count_nb <= max_nb:
                    bonds.add(bond)
                # If bond does not meet criteria flag the bond to not create the bond
                else:
                    flagged_bonds.add(bond)
                               
        # Removing any flagged bonds if they were created      
        bonds = list(set(bonds) - set(flagged_bonds)) 
       
        # Sorting bonds and flagged bonds                
        bonds = sorted(bonds); flagged_bonds = sorted(flagged_bonds) 
        self.bonds = bonds; self.flagged_bonds = flagged_bonds
       
        
        # generate bondtype_stats_holder
        bondtype_stats_holder = {} # { tuple(element1, element2) : [lst of distances]}
        for id1, id2 in bonds:
            distance = possible_bonds[ tuple(sorted([id1, id2])) ]
            element1 = m.atoms[id1].element
            element2 = m.atoms[id2].element
            bondtype = tuple(sorted([element1, element2]))
            
            # Log bondtype
            if bondtype in bondtype_stats_holder:
                bondtype_stats_holder[bondtype].append(distance)
            else:
                bondtype_stats_holder[bondtype] = [distance]
        
        # Find bond order stats
        for i in bondtype_stats_holder:
            lst = bondtype_stats_holder[i]
    
            # Find bond name
            bond = '{}-{}'.format(i[0], i[1])
            
            # find distance cutoff used
            cutoff = maxdist_dict[(i[0], i[1])]
            
            # Add to stats if count if larger then zero
            if len(lst) > 0:
                s = Stats()
                s.count = int(len(lst))
                s.avg = '{:.4f}'.format( compute_mean(lst) )
                s.min = '{:.4f}'.format( min(lst) )
                s.max = '{:.4f}'.format( max(lst) )
                s.std = '{:.4f}'.format( compute_standard_deviation(lst) )
                s.cutoff = cutoff
                self.statistics[bond] = s
                
        # Find number of bonded stats
        nb_dict = {i:[] for i in m.atoms}
        for id1, id2 in bonds:
            nb_dict[id1].append(id2)
            nb_dict[id2].append(id2)
            
        # Find unique elements and generate nb_count dict to add to
        elements = sorted({m.atoms[i].element for i in m.atoms})
        self.maxbond = max([len(nb_dict[i]) for i in nb_dict]) # find max bond from nb
        for i in elements:
            self.nb_count[i] = {j:0 for j in range(self.maxbond+1)}
            
        # Add to self.nb_count
        for i in nb_dict:
            element = m.atoms[i].element
            nb = len(nb_dict[i])
            self.nb_count[element][nb] += 1