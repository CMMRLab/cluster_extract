# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 2.1
November 3rd, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

read_REAXC 1.1 is exactly the same as read_REAXC 1.0, but has had
the create_bonds function optimized for time.

read_REAXC 1.2 has been modified from 1.1 and has added an atoms
dictionary to the read function to store each line of info based 
on the atom id for quicker acces to information like atom type.
The atom .type information can be found more quickly and is used
to find each bond type in the create_bonds function. Each bond 
type then is given a minimum bond order cut off based on the bond
type to decide if the time averaged bond order for the bond meets
the user input and then creates or flags to bond depending on time
averaged BO vs minimum bond order cut of for the specific bond type.

read_REAXC 1.3 adds some statiticis encaposolateds in the class avgs 

read_REACx 1.4 changed max bonded atoms from 6 to 9 for SMP10

read_REAXC 2.0 was a complete re-write of 1.4 to add cleaner and more
read able coding to the function. Both the read class and create bonds
class was completely overhauled for a more cleaner code to read and 
add onto if need. That max bonded atoms for 2.0 is now 10 and the create
bonds class also can handle 10. The create bonds class also is easier to
adjust then previous function for if more bonds then 10 are needed.

read_REAXC 2.1 added max bonded atoms to 20 like how many LAMMPS can handle. 
"""

####################
# Import Libraries #
####################
from collections import OrderedDict
import operator
import numpy as np


######################
# Class to Read file #
######################
class Atoms: # bond_order_file is build out to id10 and bo10
    pass # .type .maxnb .id1 .id2 .idN ... .mol .bo1 .bo2 .boN .nlp .q

class read_BO_file:
    def __init__(self, filename):
        self.steps = {} # {step value : atoms object}
        
        # Open and read file
        with open(filename, 'r') as f:
            
            # Itializatizing flags
            timestep_flag = False
            bonds_flag = False
            
            # Loop through file
            for line in f:
                
                # split and strip line
                line = line.strip()
                line = line.split()
                
                # Setting flags
                if len(line) == 0 or len(line) == 1:
                    timestep_flag = False
                    bonds_flag = False
                elif 'Timestep' in line:
                    timestep_flag = True
                elif line[0] != '#':
                    timestep_flag = False
                    bonds_flag = True
                
                # find timestep
                if timestep_flag:  
                    step = int(line[2])
                    
                    # create temporary atoms dictionary to reset
                    # every time a new timestep is iterated over
                    atoms = {} # { atom id : atoms object }
                        
                # find bond data
                elif bonds_flag:                        
                    # Find values that are always constant
                    atomid = int(line[0])
                    atomtype = int(line[1])
                    nb = int(line[2])
                    
                    # Find bonding parameters based on nb
                    # If 0 bonded atoms to atomid
                    if nb == 0:
                        a = Atoms()
                        
                        a.id = atomid
                        a.type = atomtype
                        a.nb =  nb
                        
                        a.id1 = int(0)
                        a.id2 = int(0)
                        a.id3 = int(0)
                        a.id4 = int(0)
                        a.id5 = int(0)
                        a.id6 = int(0)
                        a.id7 = int(0)
                        a.id8 = int(0)
                        a.id9 = int(0)
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[3])
                        
                        a.bo1 = int(0)
                        a.bo2 = int(0)
                        a.bo3 = int(0)
                        a.bo4 = int(0)
                        a.bo5 = int(0)
                        a.bo6 = int(0)
                        a.bo7 = int(0)
                        a.bo8 = int(0)
                        a.bo9 = int(0)
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = int(0)
                        a.nlp = int(0)
                        a.q = int(0)
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 1 bonded atoms to atomid
                    elif nb == 1:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(0)
                        a.id3 = int(0)
                        a.id4 = int(0)
                        a.id5 = int(0)
                        a.id6 = int(0)
                        a.id7 = int(0)
                        a.id8 = int(0)
                        a.id9 = int(0)
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[4])
                        
                        a.bo1 = float(line[5])
                        a.bo2 = int(0)
                        a.bo3 = int(0)
                        a.bo4 = int(0)
                        a.bo5 = int(0)
                        a.bo6 = int(0)
                        a.bo7 = int(0)
                        a.bo8 = int(0)
                        a.bo9 = int(0)
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[6])
                        a.nlp = float(line[7])
                        a.q = float(line[8])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 2 bonded atoms to atomid
                    elif nb == 2:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(0)
                        a.id4 = int(0)
                        a.id5 = int(0)
                        a.id6 = int(0)
                        a.id7 = int(0)
                        a.id8 = int(0)
                        a.id9 = int(0)
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[5])
                        
                        a.bo1 = float(line[6])
                        a.bo2 = float(line[7])
                        a.bo3 = int(0)
                        a.bo4 = int(0)
                        a.bo5 = int(0)
                        a.bo6 = int(0)
                        a.bo7 = int(0)
                        a.bo8 = int(0)
                        a.bo9 = int(0)
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[8])
                        a.nlp = float(line[9])
                        a.q = float(line[10])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 3 bonded atoms to atomid
                    elif nb == 3:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(0)
                        a.id5 = int(0)
                        a.id6 = int(0)
                        a.id7 = int(0)
                        a.id8 = int(0)
                        a.id9 = int(0)
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[6])
                        
                        a.bo1 = float(line[7])
                        a.bo2 = float(line[8])
                        a.bo3 = float(line[9])
                        a.bo4 = int(0)
                        a.bo5 = int(0)
                        a.bo6 = int(0)
                        a.bo7 = int(0)
                        a.bo8 = int(0)
                        a.bo9 = int(0)
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[10])
                        a.nlp = float(line[11])
                        a.q = float(line[12])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 4 bonded atoms to atomid
                    elif nb == 4:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(0)
                        a.id6 = int(0)
                        a.id7 = int(0)
                        a.id8 = int(0)
                        a.id9 = int(0)
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[7])
                        
                        a.bo1 = float(line[8])
                        a.bo2 = float(line[9])
                        a.bo3 = float(line[10])
                        a.bo4 = float(line[11])
                        a.bo5 = int(0)
                        a.bo6 = int(0)
                        a.bo7 = int(0)
                        a.bo8 = int(0)
                        a.bo9 = int(0)
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[12])
                        a.nlp = float(line[13])
                        a.q = float(line[14])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 5 bonded atoms to atomid
                    elif nb == 5:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(0)
                        a.id7 = int(0)
                        a.id8 = int(0)
                        a.id9 = int(0)
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[8])
                        
                        a.bo1 = float(line[9])
                        a.bo2 = float(line[10])
                        a.bo3 = float(line[11])
                        a.bo4 = float(line[12])
                        a.bo5 = float(line[13])
                        a.bo6 = int(0)
                        a.bo7 = int(0)
                        a.bo8 = int(0)
                        a.bo9 = int(0)
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[14])
                        a.nlp = float(line[15])
                        a.q = float(line[16])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 6 bonded atoms to atomid
                    elif nb == 6:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(0)
                        a.id8 = int(0)
                        a.id9 = int(0)
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[9])
                        
                        a.bo1 = float(line[10])
                        a.bo2 = float(line[11])
                        a.bo3 = float(line[12])
                        a.bo4 = float(line[13])
                        a.bo5 = float(line[14])
                        a.bo6 = float(line[15])
                        a.bo7 = int(0)
                        a.bo8 = int(0)
                        a.bo9 = int(0)
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[16])
                        a.nlp = float(line[17])
                        a.q = float(line[18])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 7 bonded atoms to atomid
                    elif nb == 7:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(0)
                        a.id9 = int(0)
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[10])
                        
                        a.bo1 = float(line[11])
                        a.bo2 = float(line[12])
                        a.bo3 = float(line[13])
                        a.bo4 = float(line[14])
                        a.bo5 = float(line[15])
                        a.bo6 = float(line[16])
                        a.bo7 = float(line[17])
                        a.bo8 = int(0)
                        a.bo9 = int(0)
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[18])
                        a.nlp = float(line[19])
                        a.q = float(line[20])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 8 bonded atoms to atomid
                    elif nb == 8:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(0)
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[11])
                        
                        a.bo1 = float(line[12])
                        a.bo2 = float(line[13])
                        a.bo3 = float(line[14])
                        a.bo4 = float(line[15])
                        a.bo5 = float(line[16])
                        a.bo6 = float(line[17])
                        a.bo7 = float(line[18])
                        a.bo8 = float(line[19])
                        a.bo9 = int(0)
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[20])
                        a.nlp = float(line[21])
                        a.q = float(line[22])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 9 bonded atoms to atomid
                    elif nb == 9:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(0)
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[12])
                        
                        a.bo1 = float(line[13])
                        a.bo2 = float(line[14])
                        a.bo3 = float(line[15])
                        a.bo4 = float(line[16])
                        a.bo5 = float(line[17])
                        a.bo6 = float(line[18])
                        a.bo7 = float(line[19])
                        a.bo8 = float(line[20])
                        a.bo9 = float(line[21])
                        a.bo10 = int(0)
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[22])
                        a.nlp = float(line[23])
                        a.q = float(line[24])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 10 bonded atoms to atomid
                    elif nb == 10:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(0)
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[13])
                        
                        a.bo1 = float(line[14])
                        a.bo2 = float(line[15])
                        a.bo3 = float(line[16])
                        a.bo4 = float(line[17])
                        a.bo5 = float(line[18])
                        a.bo6 = float(line[19])
                        a.bo7 = float(line[20])
                        a.bo8 = float(line[21])
                        a.bo9 = float(line[22])
                        a.bo10 = float(line[23])
                        a.bo11 = int(0)
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[24])
                        a.nlp = float(line[25])
                        a.q = float(line[26])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 11 bonded atoms to atomid
                    elif nb == 11:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(0)
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[14])
                        
                        a.bo1 = float(line[15])
                        a.bo2 = float(line[16])
                        a.bo3 = float(line[17])
                        a.bo4 = float(line[18])
                        a.bo5 = float(line[19])
                        a.bo6 = float(line[20])
                        a.bo7 = float(line[21])
                        a.bo8 = float(line[22])
                        a.bo9 = float(line[23])
                        a.bo10 = float(line[24])
                        a.bo11 = float(line[25])
                        a.bo12 = int(0)
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[26])
                        a.nlp = float(line[27])
                        a.q = float(line[28])
                        
                    # Find bonding parameters based on nb
                    # If 12 bonded atoms to atomid
                    elif nb == 12:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(line[14])
                        a.id13 = int(0)
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[15])
                        
                        a.bo1 = float(line[16])
                        a.bo2 = float(line[17])
                        a.bo3 = float(line[18])
                        a.bo4 = float(line[19])
                        a.bo5 = float(line[20])
                        a.bo6 = float(line[21])
                        a.bo7 = float(line[22])
                        a.bo8 = float(line[23])
                        a.bo9 = float(line[24])
                        a.bo10 = float(line[25])
                        a.bo11 = float(line[26])
                        a.bo12 = float(line[27])
                        a.bo13 = int(0)
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[28])
                        a.nlp = float(line[29])
                        a.q = float(line[30])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 13 bonded atoms to atomid
                    elif nb == 13:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(line[14])
                        a.id13 = int(line[15])
                        a.id14 = int(0)
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[16])
                        
                        a.bo1 = float(line[17])
                        a.bo2 = float(line[18])
                        a.bo3 = float(line[19])
                        a.bo4 = float(line[20])
                        a.bo5 = float(line[21])
                        a.bo6 = float(line[22])
                        a.bo7 = float(line[23])
                        a.bo8 = float(line[24])
                        a.bo9 = float(line[25])
                        a.bo10 = float(line[26])
                        a.bo11 = float(line[27])
                        a.bo12 = float(line[28])
                        a.bo13 = float(line[29])
                        a.bo14 = int(0)
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[30])
                        a.nlp = float(line[31])
                        a.q = float(line[32])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                        
                    # Find bonding parameters based on nb
                    # If 14 bonded atoms to atomid
                    elif nb == 14:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(line[14])
                        a.id13 = int(line[15])
                        a.id14 = int(line[16])
                        a.id15 = int(0)
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[17])
                        
                        a.bo1 = float(line[18])
                        a.bo2 = float(line[19])
                        a.bo3 = float(line[20])
                        a.bo4 = float(line[21])
                        a.bo5 = float(line[22])
                        a.bo6 = float(line[23])
                        a.bo7 = float(line[24])
                        a.bo8 = float(line[25])
                        a.bo9 = float(line[26])
                        a.bo10 = float(line[27])
                        a.bo11 = float(line[28])
                        a.bo12 = float(line[29])
                        a.bo13 = float(line[30])
                        a.bo14 = float(line[31])
                        a.bo15 = int(0)
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[32])
                        a.nlp = float(line[33])
                        a.q = float(line[34])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 15 bonded atoms to atomid
                    elif nb == 15:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(line[14])
                        a.id13 = int(line[15])
                        a.id14 = int(line[16])
                        a.id15 = int(line[17])
                        a.id16 = int(0)
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[18])
                        
                        a.bo1 = float(line[19])
                        a.bo2 = float(line[20])
                        a.bo3 = float(line[21])
                        a.bo4 = float(line[22])
                        a.bo5 = float(line[23])
                        a.bo6 = float(line[24])
                        a.bo7 = float(line[25])
                        a.bo8 = float(line[26])
                        a.bo9 = float(line[27])
                        a.bo10 = float(line[28])
                        a.bo11 = float(line[29])
                        a.bo12 = float(line[30])
                        a.bo13 = float(line[31])
                        a.bo14 = float(line[32])
                        a.bo15 = float(line[33])
                        a.bo16 = int(0)
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[34])
                        a.nlp = float(line[35])
                        a.q = float(line[36])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 16 bonded atoms to atomid
                    elif nb == 16:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(line[14])
                        a.id13 = int(line[15])
                        a.id14 = int(line[16])
                        a.id15 = int(line[17])
                        a.id16 = int(line[18])
                        a.id17 = int(0)
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[19])
                        
                        a.bo1 = float(line[20])
                        a.bo2 = float(line[21])
                        a.bo3 = float(line[22])
                        a.bo4 = float(line[23])
                        a.bo5 = float(line[24])
                        a.bo6 = float(line[25])
                        a.bo7 = float(line[26])
                        a.bo8 = float(line[27])
                        a.bo9 = float(line[28])
                        a.bo10 = float(line[29])
                        a.bo11 = float(line[30])
                        a.bo12 = float(line[31])
                        a.bo13 = float(line[32])
                        a.bo14 = float(line[33])
                        a.bo15 = float(line[34])
                        a.bo16 = float(line[35])
                        a.bo17 = int(0)
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[36])
                        a.nlp = float(line[37])
                        a.q = float(line[38])
                        
                    # Find bonding parameters based on nb
                    # If 17 bonded atoms to atomid
                    elif nb == 17:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(line[14])
                        a.id13 = int(line[15])
                        a.id14 = int(line[16])
                        a.id15 = int(line[17])
                        a.id16 = int(line[18])
                        a.id17 = int(line[19])
                        a.id18 = int(0)
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[20])
                        
                        a.bo1 = float(line[21])
                        a.bo2 = float(line[22])
                        a.bo3 = float(line[23])
                        a.bo4 = float(line[24])
                        a.bo5 = float(line[25])
                        a.bo6 = float(line[26])
                        a.bo7 = float(line[27])
                        a.bo8 = float(line[28])
                        a.bo9 = float(line[29])
                        a.bo10 = float(line[30])
                        a.bo11 = float(line[31])
                        a.bo12 = float(line[32])
                        a.bo13 = float(line[33])
                        a.bo14 = float(line[34])
                        a.bo15 = float(line[35])
                        a.bo16 = float(line[36])
                        a.bo17 = float(line[37])
                        a.bo18 = int(0)
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[38])
                        a.nlp = float(line[39])
                        a.q = float(line[40])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                        
                    # Find bonding parameters based on nb
                    # If 18 bonded atoms to atomid
                    elif nb == 18:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(line[14])
                        a.id13 = int(line[15])
                        a.id14 = int(line[16])
                        a.id15 = int(line[17])
                        a.id16 = int(line[18])
                        a.id17 = int(line[19])
                        a.id18 = int(line[20])
                        a.id19 = int(0)
                        a.id20 = int(0)

                        a.mol = int(line[21])
                        
                        a.bo1 = float(line[22])
                        a.bo2 = float(line[23])
                        a.bo3 = float(line[24])
                        a.bo4 = float(line[25])
                        a.bo5 = float(line[26])
                        a.bo6 = float(line[27])
                        a.bo7 = float(line[28])
                        a.bo8 = float(line[29])
                        a.bo9 = float(line[30])
                        a.bo10 = float(line[31])
                        a.bo11 = float(line[32])
                        a.bo12 = float(line[33])
                        a.bo13 = float(line[34])
                        a.bo14 = float(line[35])
                        a.bo15 = float(line[36])
                        a.bo16 = float(line[37])
                        a.bo17 = float(line[38])
                        a.bo18 = float(line[39])
                        a.bo19 = int(0)
                        a.bo20 = int(0)
                        
                        a.abo = float(line[40])
                        a.nlp = float(line[41])
                        a.q = float(line[42])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 19 bonded atoms to atomid
                    elif nb == 19:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(line[14])
                        a.id13 = int(line[15])
                        a.id14 = int(line[16])
                        a.id15 = int(line[17])
                        a.id16 = int(line[18])
                        a.id17 = int(line[19])
                        a.id18 = int(line[20])
                        a.id19 = int(line[21])
                        a.id20 = int(0)

                        a.mol = int(line[22])
                        
                        a.bo1 = float(line[23])
                        a.bo2 = float(line[24])
                        a.bo3 = float(line[25])
                        a.bo4 = float(line[26])
                        a.bo5 = float(line[27])
                        a.bo6 = float(line[28])
                        a.bo7 = float(line[29])
                        a.bo8 = float(line[30])
                        a.bo9 = float(line[31])
                        a.bo10 = float(line[32])
                        a.bo11 = float(line[33])
                        a.bo12 = float(line[34])
                        a.bo13 = float(line[35])
                        a.bo14 = float(line[36])
                        a.bo15 = float(line[37])
                        a.bo16 = float(line[38])
                        a.bo17 = float(line[39])
                        a.bo18 = float(line[40])
                        a.bo19 = float(line[41])
                        a.bo20 = int(0)
                        
                        a.abo = float(line[42])
                        a.nlp = float(line[43])
                        a.q = float(line[44])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                        
                    # Find bonding parameters based on nb
                    # If 20 bonded atoms to atomid
                    elif nb == 20:
                        a = Atoms()
                        
                        a.id = int(line[0])
                        a.type = int(line[1])
                        a.nb =  int(line[2])
                        
                        a.id1 = int(line[3])
                        a.id2 = int(line[4])
                        a.id3 = int(line[5])
                        a.id4 = int(line[6])
                        a.id5 = int(line[7])
                        a.id6 = int(line[8])
                        a.id7 = int(line[9])
                        a.id8 = int(line[10])
                        a.id9 = int(line[11])
                        a.id10 = int(line[12])
                        a.id11 = int(line[13])
                        a.id12 = int(line[14])
                        a.id13 = int(line[15])
                        a.id14 = int(line[16])
                        a.id15 = int(line[17])
                        a.id16 = int(line[18])
                        a.id17 = int(line[19])
                        a.id18 = int(line[20])
                        a.id19 = int(line[21])
                        a.id20 = int(line[22])

                        a.mol = int(line[23])
                        
                        a.bo1 = float(line[24])
                        a.bo2 = float(line[25])
                        a.bo3 = float(line[26])
                        a.bo4 = float(line[27])
                        a.bo5 = float(line[28])
                        a.bo6 = float(line[29])
                        a.bo7 = float(line[30])
                        a.bo8 = float(line[31])
                        a.bo9 = float(line[32])
                        a.bo10 = float(line[33])
                        a.bo11 = float(line[34])
                        a.bo12 = float(line[35])
                        a.bo13 = float(line[36])
                        a.bo14 = float(line[37])
                        a.bo15 = float(line[38])
                        a.bo16 = float(line[39])
                        a.bo17 = float(line[40])
                        a.bo18 = float(line[41])
                        a.bo19 = float(line[42])
                        a.bo20 = float(line[43])
                        
                        a.abo = float(line[44])
                        a.nlp = float(line[45])
                        a.q = float(line[46])
                        
                        # Add to tmp atoms dictionary
                        atoms[atomid] = a
                    
                    # Raise exception of of 20 bonded atoms
                    else:
                        raise Exception('Max number of bonds per atom was larger then 20. read_REAXC needs to be added onto to handle that many bonded atoms')
                        
                    # Add atoms dict to step
                    self.steps[step] = atoms
                
#####################################################
# Function to create bonds by averaging bond orders #
# over timesteps and applying bonding constraints   #
#####################################################

class Stats:
    pass # .count .avg .min .max .std .cutoff
    
class avgs:
    pass # .nlps .abos
        
class create_bonds:
    def __init__(self, filename, bond_info_dict, bond_cut):
        self.timesteps = []      # list of timesteps found
        self.bonds = []          # list of found bonds
        self.flaggged_bonds = [] # list of flagged bonds
        self.statistics = {}     # { tuple bond type : stats object }
        self.abo_stats = {}      # { element symbol: stats object }
        
        
        #############
        # Read file #
        #############
        BO = read_BO_file(filename)
        self.timesteps = sorted(list(BO.steps.keys()))
        
        
        ##############################################################
        # Convert bond_cut into order tuple key of bonding elements  #
        # and build dictionaries to store found bond statistics      #
        ##############################################################
        bond_cut_ordered = {};      # { tuple of bonding elements : BO cutoff }
        bondtype_stats_holder = {}; # { tuple of bonding elements : [list of avgs] }
        for bond in bond_cut:
            # Dont modify key it is the unknown default key
            if bond != 'unknown':
                sorted_bond = tuple(sorted(bond))
            else:
                sorted_bond = bond
            # Build new bond_cut_ordered key dict and bond stats dict 
            bond_cut_ordered[sorted_bond] = bond_cut[bond]
            bondtype_stats_holder[sorted_bond] = []
        bondtype_stats_holder['override'] = []; bondtype_stats_holder['unknown'] = []
        
        
        #################################################
        # Set division factor since BO will appear      #
        # twice per each step set of atoms 2*ntimesteps #
        #################################################
        div = 2*len(self.timesteps)
               
        
        #################################################
        # Loop through read file to create intial bonds #
        # dictionary with list as value append BO's to  #
        #################################################
        bond_dict = {} # { tuple(id1, id2): [] }
        nlps = {} # { atom id: [list of nlps ]}
        abos = {} # { atom id: [list of abos ]}
        for i in BO.steps:
            step = BO.steps[i]
            
            # Loop through atoms info of step to create data structures
            for j in step:
                atom = step[j]
                
                # Find intial guess of bonds from the file
                bond1 = tuple(sorted([j, atom.id1])); BO1 = atom.bo1;
                bond2 = tuple(sorted([j, atom.id2])); BO2 = atom.bo2;
                bond3 = tuple(sorted([j, atom.id3])); BO3 = atom.bo3;
                bond4 = tuple(sorted([j, atom.id4])); BO4 = atom.bo4;
                bond5 = tuple(sorted([j, atom.id5])); BO5 = atom.bo5;
                bond6 = tuple(sorted([j, atom.id6])); BO6 = atom.bo6;
                bond7 = tuple(sorted([j, atom.id7])); BO7 = atom.bo7;
                bond8 = tuple(sorted([j, atom.id8])); BO8 = atom.bo8;
                bond9 = tuple(sorted([j, atom.id9])); BO9 = atom.bo9;
                bond10 = tuple(sorted([j, atom.id10])); BO10 = atom.bo10;
                
                bond11 = tuple(sorted([j, atom.id11])); BO11 = atom.bo11;
                bond12 = tuple(sorted([j, atom.id12])); BO12 = atom.bo12;
                bond13 = tuple(sorted([j, atom.id13])); BO13 = atom.bo13;
                bond14 = tuple(sorted([j, atom.id14])); BO14 = atom.bo14;
                bond15 = tuple(sorted([j, atom.id15])); BO15 = atom.bo15;
                bond16 = tuple(sorted([j, atom.id16])); BO16 = atom.bo16;
                bond17 = tuple(sorted([j, atom.id17])); BO17 = atom.bo17;
                bond18 = tuple(sorted([j, atom.id18])); BO18 = atom.bo18;
                bond19 = tuple(sorted([j, atom.id19])); BO19 = atom.bo19;
                bond20 = tuple(sorted([j, atom.id20])); BO20 = atom.bo20;

                
                # Add intial guess of bonds bond_dict if BON != 0:
                # since 0's are used as place holders for indexing
                if BO1 != 0: bond_dict[bond1] = []
                if BO2 != 0: bond_dict[bond2] = []
                if BO3 != 0: bond_dict[bond3] = []
                if BO4 != 0: bond_dict[bond4] = []
                if BO5 != 0: bond_dict[bond5] = []
                if BO6 != 0: bond_dict[bond6] = []
                if BO7 != 0: bond_dict[bond7] = []
                if BO8 != 0: bond_dict[bond8] = []
                if BO9 != 0: bond_dict[bond9] = []
                if BO10 != 0: bond_dict[bond10] = []
                
                if BO11 != 0: bond_dict[bond11] = []
                if BO12 != 0: bond_dict[bond12] = []
                if BO13 != 0: bond_dict[bond13] = []
                if BO14 != 0: bond_dict[bond14] = []
                if BO15 != 0: bond_dict[bond15] = []
                if BO16 != 0: bond_dict[bond16] = []
                if BO17 != 0: bond_dict[bond17] = []
                if BO18 != 0: bond_dict[bond18] = []
                if BO19 != 0: bond_dict[bond19] = []
                if BO20 != 0: bond_dict[bond20] = []
                
                # create nlp and abo to dictionaries
                nlps[j] = []; abos[j] = [];
        
        
        ###################################################
        # Loop through read file a second time add start  #
        # appeding bond orders to list in dict value to   #
        # perfrom BO averaging of bond order data         #
        ###################################################
        bond_BO_avg = {} # { tuple(id1, id2): bond order average }
        for i in BO.steps:
            step = BO.steps[i]
            
            # Loop through atoms info of step
            for j in step:
                atom = step[j]
                
                # Find intial guess of bonds from the file
                bond1 = tuple(sorted([j, atom.id1])); BO1 = atom.bo1;
                bond2 = tuple(sorted([j, atom.id2])); BO2 = atom.bo2;
                bond3 = tuple(sorted([j, atom.id3])); BO3 = atom.bo3;
                bond4 = tuple(sorted([j, atom.id4])); BO4 = atom.bo4;
                bond5 = tuple(sorted([j, atom.id5])); BO5 = atom.bo5;
                bond6 = tuple(sorted([j, atom.id6])); BO6 = atom.bo6;
                bond7 = tuple(sorted([j, atom.id7])); BO7 = atom.bo7;
                bond8 = tuple(sorted([j, atom.id8])); BO8 = atom.bo8;
                bond9 = tuple(sorted([j, atom.id9])); BO9 = atom.bo9;
                bond10 = tuple(sorted([j, atom.id10])); BO10 = atom.bo10;
                
                bond11 = tuple(sorted([j, atom.id11])); BO11 = atom.bo11;
                bond12 = tuple(sorted([j, atom.id12])); BO12 = atom.bo12;
                bond13 = tuple(sorted([j, atom.id13])); BO13 = atom.bo13;
                bond14 = tuple(sorted([j, atom.id14])); BO14 = atom.bo14;
                bond15 = tuple(sorted([j, atom.id15])); BO15 = atom.bo15;
                bond16 = tuple(sorted([j, atom.id16])); BO16 = atom.bo16;
                bond17 = tuple(sorted([j, atom.id17])); BO17 = atom.bo17;
                bond18 = tuple(sorted([j, atom.id18])); BO18 = atom.bo18;
                bond19 = tuple(sorted([j, atom.id19])); BO19 = atom.bo19;
                bond20 = tuple(sorted([j, atom.id20])); BO20 = atom.bo20;
                
                # Add intial guess of bonds bond_dict if BON != 0:
                # since 0's are used as place holders for indexing
                if BO1 != 0: bond_dict[bond1].append(BO1)
                if BO2 != 0: bond_dict[bond2].append(BO2)
                if BO3 != 0: bond_dict[bond3].append(BO3)
                if BO4 != 0: bond_dict[bond4].append(BO4)
                if BO5 != 0: bond_dict[bond5].append(BO5)
                if BO6 != 0: bond_dict[bond6].append(BO6)
                if BO7 != 0: bond_dict[bond7].append(BO7)
                if BO8 != 0: bond_dict[bond8].append(BO8)
                if BO9 != 0: bond_dict[bond9].append(BO9)
                if BO10 != 0: bond_dict[bond10].append(BO10)
                
                if BO11 != 0: bond_dict[bond11].append(BO11)
                if BO12 != 0: bond_dict[bond12].append(BO12)
                if BO13 != 0: bond_dict[bond13].append(BO13)
                if BO14 != 0: bond_dict[bond14].append(BO14)
                if BO15 != 0: bond_dict[bond15].append(BO15)
                if BO16 != 0: bond_dict[bond16].append(BO16)
                if BO17 != 0: bond_dict[bond17].append(BO17)
                if BO18 != 0: bond_dict[bond18].append(BO18)
                if BO19 != 0: bond_dict[bond19].append(BO19)
                if BO20 != 0: bond_dict[bond20].append(BO20)
                
                # Add nlp and abo to dictionaries
                nlps[j].append(atom.nlp); abos[j].append(atom.abo);
            
            
        ##################################################
        # Average together Bond orders from all read in  #
        # timesteps of Bond order data bond atom id pair #
        ##################################################
        bond_BO_avg = {} # { tuple(id1, id2): bond order average }
        for i in bond_dict:
            # averaging set by div factor based on 2*ntimesteps
            bond_BO_avg[i] = sum(bond_dict[i])/div 

        
        
        ####################################
        # Find bonded atoms to be able to  #
        # set max nb cut-off per atom type #
        ####################################
        # Intialize atoms dictionary
        bonded = {}; # { atomid : set([of bonded atoms])}
        atomtype = {}; # { atomid: atomtype numeric id }
        atom_abo = {} # { atomid : sum of bond orders}
        flag_abo = {} # { atomid : [list of counted atomids]}
        for i in BO.steps:
            step = BO.steps[i]
            # Loop through atoms info of step
            for j in step:
                atom = step[j]
                bonded[j] = set([])
                atomtype[j] = atom.type
                atom_abo[j] = 0
                flag_abo[j] = []
        
        # add in connect atoms to bonded dictionary
        for id1, id2 in bond_BO_avg:     
            bonded[id1].add(id2)
            bonded[id2].add(id1)
        # ordering dictionary    
        bonded = dict(OrderedDict(sorted(bonded.items())))


        ########################################
        # Start finding and flagging bonds     #
        # that meet cut-off values set by user #
        ########################################
        bonds = set([]); flagged_bonds = set([]);
        atom_elements_abo = {} # {element: [list of abo's]}
        for i in bonded:
            connected = bonded[i]
            
            # Creating dictionary to sort in decending order
            BO_dict = {}
            for j in connected:
                bond = tuple(sorted([i,j])) 
                BO_dict[bond] = bond_BO_avg[bond]

            # Ordering dictionary by value (averaged bond order number) in descending order
            BO_dict = dict(sorted(BO_dict.items(), key=operator.itemgetter(1),reverse=True))
            
            # Find element and max_nb based on atom type
            element, max_nb = bond_info_dict[atomtype[i]]
            
            # Creating bonds based on inputs values
            count_nb = 0;
            for id1, id2 in BO_dict:
                count_nb += 1
                avg_BO = BO_dict[(id1, id2)]
                bond = tuple(sorted([id1, id2]))
                
                # Find elements 1 and 2 to determine bond type
                element_1, max_nb_1 = bond_info_dict[atomtype[id1]]
                element_2, max_nb_2 = bond_info_dict[atomtype[id2]]
                
                # Add elements to atom_elements_abo for later tallying
                atom_elements_abo[element_1] = []
                atom_elements_abo[element_2] = []
                
                # Find minimum bond order cut-off (try forward and reverse orders
                # then try for 'unknown' key. If all fails set as 0.3 and warn)
                bondtype = tuple(sorted([element_1, element_2]))
                if bondtype in bond_cut:
                    min_BO = bond_cut_ordered[bondtype]
                    storage_type = bondtype
                    #bondtype_stats_holder[bondtype].append(avg_BO)
                elif 'unknown' in bond_cut:
                    min_BO = bond_cut['unknown']
                    storage_type = 'unknown'
                    #bondtype_stats_holder['unknown'].append(avg_BO)
                else:
                    min_BO = 0.3
                    storage_type = 'override'
                    #bondtype_stats_holder['override'].append(avg_BO)
                    print(f'WARNING: bonding elements {element_1}-{element_2} are not in bond_cut dictionary.')
                    print('User defined BO cut-off was overrode for a cut-off of 0.3 (standard cut-off value)')
                    
                
                # If number of bonds is less than specified and BO is higher then minimum create the bond
                if count_nb <= max_nb and avg_BO >= min_BO:
                    bonds.add(bond)
                    bondtype_stats_holder[storage_type].append(avg_BO)
                    
                    # Set atom_abo's if bond not already accounted for
                    if id2 not in flag_abo[id1]:
                        atom_abo[id1] += avg_BO
                        flag_abo[id1].append(id2)
                    if id1 not in flag_abo[id2]:
                        atom_abo[id2] += avg_BO
                        flag_abo[id2].append(id1)
                        
                # If bond does not meet criteria flag the bond to not create the bond
                else:
                    flagged_bonds.add(bond)
                            
        # Removing any flagged bonds if they were created      
        bonds = list(set(bonds) - set(flagged_bonds)) 
             
        
        # Sorting bonds and flagged bonds                
        bonds = sorted(bonds); flagged_bonds = sorted(flagged_bonds) 
        self.bonds = bonds; self.flagged_bonds = flagged_bonds;
        
        
        ###########################
        # Find bonding statistics #
        ###########################
        # Functions for finding bond stats
        def count(lst):
            return int(len(lst)/2)
        
        def avg(lst):
            if lst: a = sum(lst)/len(lst)
            else: a = 0
            return a
        
        def lo(lst):
            if lst: m = min(lst)
            else: m = 0
            return m
        
        def hi(lst):
            if lst: m = max(lst)
            else: m = 0
            return m
        
        def std(lst):
            if lst: s = np.std(np.array(lst))
            else: s = 0
            return s
            
        # Find bond order stats
        for i in bondtype_stats_holder:
            avg_BO_lst = bondtype_stats_holder[i]
    
            # Find bond name
            if bondtype != 'unknown' or bondtype != 'override':
                bond = '{}{}{}'.format(i[0], '-', i[1])
            else: 
                bond = i
                
            # Find info to save
            nbonds = count(avg_BO_lst)
            avgbonds = '{:.4f}'.format(avg(avg_BO_lst))
            lobonds = '{:.4f}'.format(lo(avg_BO_lst))
            hibonds = '{:.4f}'.format(hi(avg_BO_lst))
            stdbonds = '{:.4f}'.format(std(avg_BO_lst))
            
            # find cutoff used
            if i in bond_cut_ordered:
                cutoff = bond_cut_ordered[i]
            else:
                cutoff = 0.3
            
            # Add to stats if count if larger then zero
            if nbonds > 0:
                s = Stats()
                s.count = nbonds
                s.avg = avgbonds
                s.min = lobonds
                s.max = hibonds
                s.std = stdbonds
                s.cutoff = cutoff
                self.statistics[bond] = s
                #self.statistics[i] = s
                
        ########################
        # Find atom statistics #
        ########################
        max_nb_elemet_cutoff = {} # { element symbol: cutoff}
        for i in atom_abo:
            # Find element and max_nb based on atom type
            # and append to atom_elements_abo dictionary
            element, max_nb = bond_info_dict[atomtype[i]]            
            atom_elements_abo[element].append(atom_abo[i])
            max_nb_elemet_cutoff[element] = max_nb

        # Loop through atom_elements_abo and find statistics
        for i in atom_elements_abo:
            abo_lst = atom_elements_abo[i]
            
            # Find abo stats
            ss = Stats()
            ss.count = int(len(abo_lst))
            ss.avg = '{:.4f}'.format(avg(abo_lst))
            ss.min = '{:.4f}'.format(lo(abo_lst))
            ss.max = '{:.4f}'.format(hi(abo_lst))
            ss.std = '{:.4f}'.format(std(abo_lst))
            ss.cutoff = max_nb_elemet_cutoff[i]
            self.abo_stats[i] = ss
        
                
        ########################
        # Create avgs instance #
        ########################
        # average nlps and abos
        nlps_avg = {} # { atomid : avg nlp }
        abos_avg = {} # { atomid : avg abo }
        for i in nlps:
            nlps_avg[i] = sum(nlps[i])/len(nlps[i])
        for i in abos:
            abos_avg[i] = sum(abos[i])/len(abos[i])
            
        # Add to avg class
        a = avgs()
        a.nlps = nlps_avg
        a.abos = abos_avg
        self.avgs = a
        


###############################################    
### Testing reading reaxc and bond creation ###    
###############################################
if __name__ == '__main__':

    # Set reaxc file name
    reaxc_file = 'post_mix_bonds.reaxc'
        
    
    # Using create bonds function. Adjust these as needed for your case.
    #   {atomtype: ('element symbol', max number of bonded atoms)}  
    bond_info_dict = {1: ('C', 4), 2: ('H', 1), 3: ('O', 3)}
    
    # C H O N S possbile bonding configurations bond order cut offs. This bond order cut off will be used for the
    # time averaged bond order for each bond type. This gives the ability to adjust the minimum bond order cut off
    # for each bond type. If the bond does not exists in this bond_cut dictionary, the 'unknown' key bond order cut
    # will be used. Add onto the dictnary as needed, the code will adjust add scale based on the info provided in 
    # the bond_info list and bond_cut dictionary. The code will reoder the dictionary bond type keys before using 
    # so only one ordering of a bond pair needs to be given the the code will worry about forward and reverse 
    # ordering (IE C-H bond only needs either be listed as ('C','H') or ('H','C') and the code will look for both
    # forward and reverse orderings).
    #          like-paired      pairing of all other configurations
    bond_cut = {('C','C'): 0.5, ('C','H'): 0.3, ('C','O'): 0.3, ('C','N'): 0.3, ('C','S'): 0.3, # C-other bonds
                ('H','H'): 0.3, ('H','O'): 0.3, ('H','N'): 0.3, ('H','S'): 0.3,                 # H-other bonds 
                ('O','O'): 0.3, ('O','N'): 0.3, ('O','S'): 0.3,                                 # O-other bonds
                ('N','N'): 0.3, ('N','S'): 0.3,                                                 # N-other bonds
                ('S','S'): 0.3,                                                                 # S-other bonds
                'unknown': 0.3}                                                                 # Default
    
    # Find bonds with create_bonds function that calls the read_BO_file class
    bonds = create_bonds(reaxc_file, bond_info_dict, bond_cut)
    
    # Print info on bonds found
    print('\n\n--------------------Bonding info found--------------------')
    print('timsteps average over:', bonds.timesteps)
    print('total bonds found: ', len(bonds.bonds))
    print('total flagged bonds:', len(bonds.flagged_bonds))
    
    # Print bonding statistics
    print('\n\n--------------------------Bonding statistics found--------------------------')
    print('{:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format('Bond', 'avg', 'count', 'min BO', 'max BO', 'std', 'cutoff'))
    print('-------------------------------------------------------------------------')
    for i in bonds.statistics:
        stats = bonds.statistics[i]
        print('{:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(str(i), stats.count, stats.avg, stats.min, stats.max, stats.std, stats.cutoff))
    
    