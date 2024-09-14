# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2 Splits lines everhy 3 chars (including white space to handle triple digit atomIDs)
May 19th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This script allows for a .mol file to be read in
and an Molecule class 'm' to be created on the info
in the .mol file, such that the class is exactly the
same that read_lmp creates, when reading in a LAMMPS
data file. This allows for a .mol file to have the exact
same class type objects that exists when reading a LAMMPS
data file and allows the set of codes to read a .mol file
as if it was a LAMMPS data file.

The atom cooridinates will be centered around the geometrical
center of the read in molecule and then the box will be centered
around the newly re-centered moleucle. The box will be 0.5 angstroms
larger in the x, y, and z directions and all image flags will be set
to zero. The atom charge will be set to zero and updated later on and 
the atom molid will be set as 1 and will stay as 1 throughout the
rest of the codes that use the information found by this code.

This function has been testing on .mol file created by chemdraw
and .mol files that have been download from the NIST Chemistry
webook site. If you have a .mol file that causes an error with
this script, two option exist:
    - address the issue and re-work the read_mol class
    - open .mol file in chemdraw and save from chemdraw to create
      new file that this script should be able to read in.
      
This script could serve as a starting point to create other
_ _ _2lmp scripts to have the entire set of scripts be able
to read in other molecule file types and be able to convert
a large number of file types to a lammps datafile.  


Format of .mol files from (character spacing):
    http://biotech.fyicenter.com/1000024_SDF_File_Format_Specification.html
"""


########################################    
### Class for reading chemdraw files ###
########################################  
class Atom_mol:
    pass # .element .x .y .z

class Bond_mol:
    pass  # .type .atomids = [atom1id, atom2id]
            
# for reading .mol file
class read_mol:    
    def __init__(self, inmolfile):
        self.natoms = 0 # total atoms in .mol file
        self.nbonds = 0 # total bonds in .mol file
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        
        ####################################
        # Testing if integer function for  #
        # finding bonds line in .mol file. #
        ####################################
        def is_integer(n):
            try:
                float(n)
            except ValueError:
                return False
            else:
                return float(n).is_integer()
        
        #########################################
        ### Opening and reading the .mol file ###
        #########################################
        with open(inmolfile, 'r') as f:
            
            # Initializing flags
            info_flag = False; atoms_flag = False; bonds_flag = False;
            end_flag = False; atom_id = 0; bond_id = 0;
            
            # Looping through each line of the .mol file
            for n, line in enumerate(f):
                
                # Split every 3 chars to get natoms
                char_lst = [line[i:i+3] for i in range(0, len(line), 3)]
                
                # Setting flags
                if char_lst == ['\n']:
                    info_flag = False; atoms_flag = False;
                    bonds_flag = False; end_flag = False; 
                elif 'V2000' in line:
                    info_flag = True
                    natoms = int(char_lst[0].strip())
                elif len(char_lst) > 11 and n >= 4:
                    info_flag = False
                    atoms_flag = True
                elif line != '' and is_integer(char_lst[0].strip()) and n > 3+natoms:
                    atoms_flag = False
                    bonds_flag = True
                elif 'END' in line or 'M' in line:
                    atoms_flag = False
                    bonds_flag = False
                    end_flag = True
    
                
                # Finding atom and bond info
                if info_flag:
                    # Find char_lst be splitting every 3-characters and then split whitespace from needed info
                    char_lst = [line[i:i+3] for i in range(0, len(line), 3)] # Split every 3 chars to get natoms and nbonds
                    self.natoms = int(char_lst[0].strip())
                    self.nbonds = int(char_lst[1].strip())
                elif atoms_flag:
                    atom_id += 1
                    a = Atom_mol()
                    char_coords = [line[i:i+10] for i in range(0, len(line), 10)] # Split every 10 chars to get coords
                    char_element = line[30:33] # Split line at 30-33 chars to get element
                    a.x = float(char_coords[0].strip())
                    a.y = float(char_coords[1].strip())
                    a.z = float(char_coords[2].strip())
                    a.element = str(char_element.strip())
                    self.atoms[atom_id] = a
                elif bonds_flag:
                    char_bonds = [line[i:i+3] for i in range(0, len(line), 3)] # Split every 3 chars to get bonds
                    bond_id += 1; b = Bond_mol();
                    b.atomids = [int(char_bonds[0].strip()), int(char_bonds[1].strip())]
                    b.type = int(char_bonds[2].strip())
                    self.bonds[bond_id] = b
                elif end_flag:
                    pass



##################################################
### Class for converting all info in chemdraw  ### 
### .mol into a class that is exactly the same ###
### as it comes from read_lmp to be able to    ###
### use .mol files in remaining of the codes   ###
##################################################
class Atom:
    pass # .type .molid .charge .x .y .z .ix. iy .iz
    
class Bond:
    pass  # .type .atomids = [atom1id, atom2id]
    
import os
    
class Molecule_File:
    def __init__(self, inmolfile):

        # Read chemdraw file
        molecule = read_mol(inmolfile) 
        
        # Set filename
        self.filename = inmolfile
        
        # Set system information
        self.natoms = molecule.natoms # total atoms in .mol file
        self.nbonds = molecule.nbonds # total bonds in .mol file
        self.natomtypes = 0
        self.nbondtypes = 0
        self.atoms = {}  # {atom number : atom object}
        self.bonds = molecule.bonds  # {bond number : bond object}
        self.xbox_line = ''
        self.ybox_line = ''
        self.zbox_line = ''
        self.xy = 0;
        self.xz = 0;
        self.yz = 0;
        self.header = '{} {} {}'.format('HEADER, ', os.path.basename(inmolfile), ' read w/ mol2lmp')
        self.velocities = {} # { atomid : tuple(velx, vely, velz)}
        self.bond_coeffs = {}
        
        
        ########################
        # Find center of atoms #
        # and shift needed to  #
        # center atoms         #
        ########################
        x = []; y = []; z = [];
        for i in molecule.atoms:
            x.append(molecule.atoms[i].x)
            y.append(molecule.atoms[i].y)
            z.append(molecule.atoms[i].z)
        
        # Finding geometric center
        x_center = sum(x)/len(x); y_center = sum(y)/len(y); z_center = sum(z)/len(z)
        
        # Finding needed shift to center around geometric center
        x_shift = (0 - x_center); y_shift = (0 - y_center); z_shift = (0 - z_center);

        
        #####################################
        # Shifting atoms and building atoms #
        # object exactly as read_PCFF does  #
        #####################################
        for i in molecule.atoms:
            # adding in new atom types and shifting atoms to be centered around geometric center
            a = Atom()
            a.type = 1 # Set as 1 this information is not used by all2lmp
            a.molid = 1  # Set as 1 (option in main code to update later on)
            a.charge = 0 # Set as 0 (option in main code to update later on)
            a.element = molecule.atoms[i].element
            a.x = molecule.atoms[i].x + x_shift
            a.y = molecule.atoms[i].y + y_shift
            a.z = molecule.atoms[i].z + z_shift
            a.ix = 0 # box size will be found such that image flags from .mol files will be zero
            a.iy = 0 # box size will be found such that image flags from .mol files will be zero
            a.iz = 0 # box size will be found such that image flags from .mol files will be zero
            self.atoms[i] = a
            self.velocities[i] = (0, 0, 0)
            
        ###############################################
        # Find new atoms x, y, z position after shift #  
        # to build lammps box size to produce image   #
        # flags in x, y, and z direction as zeros     #
        ###############################################
        x = []; y = []; z = [];
        for i in self.atoms:
            x.append(self.atoms[i].x)
            y.append(self.atoms[i].y)
            z.append(self.atoms[i].z)
        
        ###################################################
        # Find x, y, z box dims (search for min/max and   #
        # then oversize slightly in each direction also   #
        # if certain dimensions are zero set default dim) #
        ###################################################
        oversize = 0.5 # default oversize of 0.25 angtroms in each value (total over size = 0.25*2 = 0.5 angstroms)
        zero_dim_default = 0.5 # default +- value of box dimension is zero
        xlo = min(x)-oversize; xhi = max(x)+oversize;
        ylo = min(y)-oversize; yhi = max(y)+oversize;
        zlo = min(z)-oversize; zhi = max(z)+oversize;
        
        # if xlo and xhi == 0 reset to +- zero_dim_default value
        if xlo == 0 and xhi == 0 or xlo == -oversize and xhi == oversize:
            xlo = -zero_dim_default
            xhi = zero_dim_default
            
        # if ylo and yhi == 0 reset to +- zero_dim_default value
        if ylo == 0 and yhi == 0 or ylo == -oversize and yhi == oversize:
            ylo = -zero_dim_default
            yhi = zero_dim_default
            
        # if zlo and zhi == 0 reset to +- zero_dim_default value
        if zlo == 0 and zhi == 0 or zlo == -oversize and zhi == oversize:
            zlo = -zero_dim_default
            zhi = zero_dim_default
        
        # Set box dimensions string
        self.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
        self.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
        self.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')
        