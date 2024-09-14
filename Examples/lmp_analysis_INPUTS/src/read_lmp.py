# -*- coding: utf-8 -*-
"""
@author: Benjamin Jensen, with modifications by Will Pisani, 
         with modification by Josh Kemppainen (image flags,
         type labels, headers, triclinic support, morse/
         harmonic support, moved mass instance to be called
         from Coeff_class, added ability to read style hints)
Revision 3.4
March 1st, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49913
"""

class Atom:
    pass # .type .molid .charge .x .y .z .ix. iy .iz .comment

class Bond:
    pass  # .type .atomids = [atom1id, atom2id]

class Angle:
    pass  # .type .atomids = [atom1id, atom2id, atom3id]

class Dihedral:
    pass  # .type .atomids = [atom1id, atom2id, atom3id, atom4id]

class Improper:
    pass  # .type .atomids = [atom1,atom2,atom3,atom4]

class Coeff_class:
    pass  # .type .coeffs = []
    
# Function to strip comments
def strip_comment(line):
    # remove comments
    end = line.find('#')
    if end >= 0:
        #comment = line[end:]
        line = line[:end]
    return line


# Class to read type labels (will be executed first so mapping can occur)
class Type_Labels:
    def __init__(self, inmolfile, method='forward'):
        self.atom_type_labels_forward = {}  # {atom type : atom type label}
        self.atom_type_labels_reverse = {}  # {atom type label: atom type}
        
        self.bond_type_labels_forward = {}  # {bond type: bond type label}
        self.bond_type_labels_reverse = {}  # {bond type label: bond type}
        
        self.angle_type_labels_forward = {}  # {angle type : angle type label}
        self.angle_type_labels_reverse = {}  # {angle type label : angle type}
        
        self.dihedral_type_labels_forward = {}  # {dihedral type : dihedral type label}
        self.dihedral_type_labels_reverse = {}  # {dihedral type label : dihedral type}
        
        self.improper_type_labels_forward = {}  # {improper type : improper type label}
        self.improper_type_labels_reverse = {}  # {improper type label: improper type}
        
        self.type_labels_flag = False # update if type labels exists
        
        # Open and read file
        with open(inmolfile, 'r') as f:
            
            # Initialize flags
            atom_flag = False
            bond_flag = False
            angle_flag = False
            dihedral_flag = False
            improper_flag = False
            skip = 0
            
            # Loop through lines of file
            for n, whole_line in enumerate(f):
                # skip comment lines
                skip -= 1
                if skip >= 0:
                    continue
    
                # remove comments
                line = strip_comment(whole_line)
                line = line.strip()
    
                # begining of a section, flag the start and skip one line
                if line == '':
                    atom_flag = False
                    bond_flag = False
                    angle_flag = False
                    dihedral_flag = False
                    improper_flag = False
                elif 'Atom Type Labels' in line:
                    atom_flag = True
                    self.type_labels_flag = True
                    skip = 1; continue;
                elif 'Bond Type Labels' in line:
                    bond_flag = True
                    self.type_labels_flag = True
                    skip = 1; continue;
                elif 'Angle Type Labels' in line:
                    angle_flag = True
                    self.type_labels_flag = True
                    skip = 1; continue;
                elif 'Dihedral Type Labels' in line:
                    dihedral_flag = True
                    self.type_labels_flag = True
                    skip = 1; continue ; 
                elif 'Improper Type Labels' in line:
                    improper_flag = True
                    self.type_labels_flag = True
                    skip = 1; continue;   

                # Find atom type labels
                if atom_flag:
                    line = line.split()
                    ID = int(line[0])
                    Type = line[1]
                    self.atom_type_labels_forward[Type] = ID  
                    self.atom_type_labels_reverse[ID] = Type     
                # Find bond type labels
                elif bond_flag:
                    line = line.split()
                    ID = int(line[0])
                    Type = line[1]
                    self.bond_type_labels_forward[Type] = ID
                    self.bond_type_labels_reverse[ID] = Type
                # Find angle type labels
                elif angle_flag:
                    line = line.split()
                    ID = int(line[0])
                    Type = line[1]
                    self.angle_type_labels_forward[Type] = ID
                    self.angle_type_labels_reverse[ID] = Type
                # Find dihedral type labels
                elif dihedral_flag:
                    line = line.split()
                    ID = int(line[0])
                    Type = line[1]
                    self.dihedral_type_labels_forward[Type] = ID
                    self.dihedral_type_labels_reverse[ID] = Type
                # Find improper type labels
                elif improper_flag:
                    line = line.split()
                    ID = int(line[0])
                    Type = line[1]
                    self.improper_type_labels_forward[Type] = ID
                    self.improper_type_labels_reverse[ID] = Type
                    
        # Check that forward and reverse dicts are the same length. If they are not the same length the label
        # types are not unique since keys were written to reverse dict more then once and raise exception
        if len(self.atom_type_labels_forward) != len(self.atom_type_labels_reverse):
            raise Exception('Atom labels are not Unique!')
        if len(self.bond_type_labels_forward) != len(self.bond_type_labels_reverse):
            raise Exception('Bond labels are not Unique!')
        if len(self.angle_type_labels_forward) != len(self.angle_type_labels_reverse):
            raise Exception('Angle labels are not Unique!')
        if len(self.dihedral_type_labels_forward) != len(self.dihedral_type_labels_reverse):
            raise Exception('Dihedral labels are not Unique!')   
        if len(self.improper_type_labels_forward) != len(self.improper_type_labels_reverse):
            raise Exception('Improper labels are not Unique!')
            
# Function to check if variable is an int
def check_integer(variable):
    try:
        int(variable)
        return_boolean = True
    except:
        return_boolean = False
    return return_boolean
                    
# Function for label type mapping
def label_type_mapping(Type, forward_map, reverse_map, method, section):
    # forward and reverse map will be used to map type to and from label type
    # method will set which mapping occurs (options: 'forward' or 'reverse')
    # 'forward' mapping will set all types as integers
    # 'reverse' mapping will set all types as sring type labels
    
    # section will determine what to do to each section to make it consistant (options: 'topology' or 'force-feild')
    # 'topology' section will be for atoms, bonds, angles, dihedrals, impropers mapping
    # 'force-feild' section will be for atoms, bonds, angles, dihedrals, impropers, and cross term coeffs
    
    # Topology section mapping (will convert topology paramters to be consistent with force-feild section str(type-label) -> int(type id) )
    if section == 'topology' and forward_map and reverse_map:
        # Foward mapping and check for integer, if it is already integer it does not need to be mapped!
        if  method == 'forward' and not check_integer(Type): Type = forward_map[Type]
        
        # Reverse mapping return string
        if  method == 'reverse': Type = str(Type)
        
    
    # Force-feild section mapping (will convert force-feild paramters to be consistent with topology section int(type id) -> str(type-label) )
    elif section == 'force-feild' and forward_map and reverse_map:
        # Foward mapping return int
        if  method == 'forward' and check_integer(Type): Type = int(Type) #reverse_map[int(type)]
        
        # Foward mapping and check for integer, if integer  set as as type label
        if  method == 'reverse' and check_integer(Type): Type = reverse_map[int(Type)]
        
    # else set type as integer if forward and revers map dicts are empty it must be a standard lammps datafile
    else:
        Type = int(Type)
    return Type           

# Function get get style hint
def get_style_hint(whole_line):
    if '#' in whole_line:
        try: style_hint = whole_line.split('#')[-1].strip()
        except: style_hint = 'N/A' 
    else: style_hint = 'N/A'
    return style_hint                  

# Class to read datafile
class Molecule_File:
    def __init__(self, inmolfile, method):
        # Find type labels if they are in file and give access to instance of each type label
        # (will loop through the entire 1st so mapping can be performed when needed, making
        # the ordering of the type files independent). Set method of mapping (forward converts
        # topologies to integere values to match coeffs and reverse converts coeffs to labels,
        # keeping topologies as label types). 'forward' mapping is recommended!
        type_labels = Type_Labels(inmolfile); self.type_labels_flag = type_labels.type_labels_flag
        self.atom_type_labels_forward = type_labels.atom_type_labels_forward
        self.atom_type_labels_reverse = type_labels.atom_type_labels_reverse
        self.bond_type_labels_forward = type_labels.bond_type_labels_forward
        self.bond_type_labels_reverse = type_labels.bond_type_labels_reverse
        self.angle_type_labels_forward = type_labels.angle_type_labels_forward
        self.angle_type_labels_reverse = type_labels.angle_type_labels_reverse
        self.dihedral_type_labels_forward = type_labels.dihedral_type_labels_forward
        self.dihedral_type_labels_reverse = type_labels.dihedral_type_labels_reverse
        self.improper_type_labels_forward = type_labels.improper_type_labels_forward 
        self.improper_type_labels_reverse = type_labels.improper_type_labels_reverse 
        
        # Filename
        self.filename = inmolfile
        
        # Structure info
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.angles = {}  # {angle number : angle object}
        self.dihedrals = {}  # {dihedral number : dihedral object}
        self.impropers = {}  # {improper number : improper object}
        self.velocities = {}  # {atom number : tuple of velocities}

        # Parameters
        self.masses = {}  # {atom type : list of coeffs}
        self.pair_coeffs = {}  # {atom type : list of coeffs}
        self.bond_coeffs = {}   # {bond type : list of coeffs}  {1: [340,1.5], 2: [450,1.2], ...}
        self.angle_coeffs = {}  # {angle type : list of coeffs}
        self.dihedral_coeffs = {}  # {dihedral type : list of coeffs}
        self.improper_coeffs = {}  # {improper type : list of coeffs}
        self.bondbond_coeffs = {} # {cross-term number : list of coeffs}, angles
        self.bondangle_coeffs = {} # {cross-term number : list of coeffs}, angles
        self.angleangle_coeffs = {} # {cross-term number : list of coeffs}, impropers
        self.angleangletorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.endbondtorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.middlebondtorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.bondbond13_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.angletorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        
        # Style hints
        self.mass_coeffs_style_hint = 'N/A'
        self.pair_coeffs_style_hint = 'N/A'
        self.bond_coeffs_style_hint = 'N/A'
        self.angle_coeffs_style_hint = 'N/A'
        self.dihedral_coeffs_style_hint = 'N/A'
        self.improper_coeffs_style_hint = 'N/A'
        self.bondbond_coeffs_style_hint = 'N/A'
        self.bondangle_coeffs_style_hint = 'N/A'
        self.angleangle_coeffs_style_hint = 'N/A'
        self.angleangletorsion_coeffs_style_hint = 'N/A'
        self.endbondtorsion_coeffs_style_hint = 'N/A'
        self.middlebondtorsion_coeffs_style_hint = 'N/A'
        self.bondbond13_coeffs_style_hint = 'N/A'
        self.angletorsion_coeffs_style_hint = 'N/A'


        self.total_line = ''
        self.ttype_line = ''
        self.xbox_line = ''
        self.ybox_line = ''
        self.zbox_line = ''
        self.xy = 0;
        self.xz = 0;
        self.yz = 0;
        self.extra_lines = ''
        self.header = ''

        self.total = 0
        self.natoms = 0
        self.natomtypes = 0
        self.nbonds = 0
        self.nbondtypes = 0
        self.nangles = 0
        self.nangletypes = 0
        self.ndihedrals = 0
        self.ndihedraltypes = 0
        self.nimpropers = 0
        self.nimpropertypes = 0
        self.nbondbond = 0
        self.nbondangle = 0
        self.nangleangle = 0
        self.nangleangletorsion = 0
        self.nendbondtorsion = 0
        self.nmiddlebondtorsion = 0
        self.nbondbond13 = 0
        self.nangletorsion = 0

        
        # Open and read file
        with open(inmolfile, 'r') as f:
            
            # Initialize flags
            coeff_flag = False
            atomflag = False
            bondflag = False
            angleflag = False
            dihedralflag = False
            improperflag = False
            velocityflag = False    
            skip = 0
            
            # Loop through lines of file
            for n, whole_line in enumerate(f):
                # skip comment lines
                skip -= 1
                if skip >= 0:
                    continue
    
                # remove comments
                line = strip_comment(whole_line)
                line = line.strip()
    
                # begining of a section, flag the start and skip one line
                if line == '':
                    coeff_flag = False
                    atomflag = False
                    bondflag = False
                    angleflag = False
                    dihedralflag = False
                    improperflag = False
                    velocityflag = False
                elif 'LAMMPS' in line or 'HEADER' in line or n == 0:
                    self.header = line
                    continue
                elif 'atoms' in line:
                    self.total_line = line
                    self.total = int(line.split()[0])
                    self.natoms = self.total
                    continue
                elif 'bonds' in line:
                    self.nbonds = int(line.split()[0])
                    continue
                elif 'angles' in line:
                    self.nangles = int(line.split()[0])
                    self.nbondbond = int(line.split()[0])
                    self.nbondangle = int(line.split()[0])
                    continue
                elif 'dihedrals' in line:
                    self.ndihedrals = int(line.split()[0])
                    self.nangleangletorsion = int(line.split()[0])
                    self.nendbondtorsion = int(line.split()[0])
                    self.nmiddlebondtorsion = int(line.split()[0])
                    self.nbondbond13 = int(line.split()[0])
                    self.nangletorsion = int(line.split()[0])
                    continue
                elif 'impropers' in line:
                    self.nimpropers = int(line.split()[0])
                    self.nangleangle = int(line.split()[0])
                    continue
                elif 'atom types' in line:
                    self.ttype_line = line
                    self.natomtypes = int(line.split()[0])
                    continue
                elif 'bond types' in line:
                    self.nbondtypes = int(line.split()[0])
                    continue
                elif 'angle types' in line:
                    self.nangletypes = int(line.split()[0])
                    continue
                elif 'dihedral types' in line:
                    self.ndihedraltypes = int(line.split()[0])
                    continue
                elif 'improper types' in line:
                    self.nimpropertypes = int(line.split()[0])
                    continue
                elif 'per atom' in line:
                    self.extra_lines += line + '\n'
                    continue
                elif 'xlo' in line:
                    self.xbox_line = line
                    continue
                elif 'ylo' in line:
                    self.ybox_line = line
                    continue
                elif 'zlo' in line:
                    self.zbox_line = line
                    continue
                elif 'xy' in line:
                    self.xy = float(line.split()[0])
                    continue   
                elif 'xz' in line:
                    self.xz = float(line.split()[1])
                    continue  
                elif 'yz' in line:
                    self.yz = float(line.split()[2])
                    continue                   
                elif line == 'Masses':
                    self.mass_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.atom_type_labels_forward
                    reverse_map = self.atom_type_labels_reverse
                    coeffs = self.masses
                    skip = 1
                    continue
                elif line == 'Pair Coeffs':
                    self.pair_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.atom_type_labels_forward
                    reverse_map = self.atom_type_labels_reverse
                    coeffs = self.pair_coeffs
                    skip = 1
                    continue
                elif line == 'Bond Coeffs':
                    self.bond_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.bond_type_labels_forward
                    reverse_map = self.bond_type_labels_reverse
                    coeffs = self.bond_coeffs
                    skip = 1
                    continue
                elif line == 'Angle Coeffs':
                    self.angle_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.angle_type_labels_forward
                    reverse_map = self.angle_type_labels_reverse
                    coeffs = self.angle_coeffs
                    skip = 1
                    continue
                elif line == 'Dihedral Coeffs':
                    self.dihedral_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.dihedral_type_labels_forward
                    reverse_map = self.dihedral_type_labels_reverse
                    coeffs = self.dihedral_coeffs
                    skip = 1
                    continue
                elif line == 'Improper Coeffs':
                    self.improper_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.improper_type_labels_forward
                    reverse_map = self.improper_type_labels_reverse
                    coeffs = self.improper_coeffs
                    skip = 1
                    continue
                elif line == 'BondBond Coeffs':
                    self.bondbond_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.angle_type_labels_forward
                    reverse_map = self.angle_type_labels_reverse
                    coeffs = self.bondbond_coeffs
                    skip = 1
                    continue
                elif line == 'BondAngle Coeffs':
                    self.bondangle_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.angle_type_labels_forward
                    reverse_map = self.angle_type_labels_reverse
                    coeffs = self.bondangle_coeffs
                    skip = 1
                    continue
                elif line == 'AngleAngle Coeffs':
                    self.angleangle_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.improper_type_labels_forward
                    reverse_map = self.improper_type_labels_reverse
                    coeffs = self.angleangle_coeffs
                    skip = 1
                    continue
                elif line == 'AngleAngleTorsion Coeffs':
                    self.angleangletorsion_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.dihedral_type_labels_forward
                    reverse_map = self.dihedral_type_labels_reverse
                    coeffs = self.angleangletorsion_coeffs
                    skip = 1
                    continue
                elif line == 'EndBondTorsion Coeffs':
                    self.endbondtorsion_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.dihedral_type_labels_forward
                    reverse_map = self.dihedral_type_labels_reverse
                    coeffs = self.endbondtorsion_coeffs
                    skip = 1
                    continue
                elif line == 'MiddleBondTorsion Coeffs':
                    self.middlebondtorsion_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.dihedral_type_labels_forward
                    reverse_map = self.dihedral_type_labels_reverse
                    coeffs = self.middlebondtorsion_coeffs
                    skip = 1
                    continue
                elif line == 'BondBond13 Coeffs':
                    self.bondbond13_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.dihedral_type_labels_forward
                    reverse_map = self.dihedral_type_labels_reverse
                    coeffs = self.bondbond13_coeffs
                    skip = 1
                    continue
                elif line == 'AngleTorsion Coeffs':
                    self.angletorsion_coeffs_style_hint = get_style_hint(whole_line)
                    coeff_flag = True
                    forward_map = self.dihedral_type_labels_forward
                    reverse_map = self.dihedral_type_labels_reverse
                    coeffs = self.angletorsion_coeffs
                    skip = 1
                    continue
                elif line == 'Atoms':
                    atomflag = True
                    atomstyle = whole_line.split()[-1]
                    skip = 1
                    continue
                elif line == 'Bonds':
                    bondflag = True
                    skip = 1
                    continue
                elif line == 'Angles':
                    angleflag = True
                    skip = 1
                    continue
                elif line == 'Dihedrals':
                    dihedralflag = True
                    skip = 1
                    continue
                elif line == 'Impropers':
                    improperflag = True
                    skip = 1
                    continue
                elif line == 'Velocities':
                    velocityflag = True
                    skip = 1
                    continue

                # Find coeffs    
                if coeff_flag:
                    line = line.split()
                    ID = label_type_mapping(line[0], forward_map, reverse_map, method, section='force-feild')
                    c = Coeff_class()
                    idcoeffs = []
                    for i in line[1:]:
                        if '.' not in i and 'e' not in i and 'class2' not in i and 'morse' not in i:
                            idcoeffs.append(int(i))
                        elif 'e' not in i and 'class2' not in i and 'morse' not in i:
                            idcoeffs.append(float(i))
                           
                        # for morse and class2 coeffs
                        elif 'class2' in i:
                            idcoeffs.append(i)
                            #print(line)
                        elif 'morse' in i:
                            idcoeffs.append(i)
                            #print(line)
                        else:
                            idcoeffs.append(float(i))
                            
                    if '#' not in whole_line:
                        c.type = 'N/A'
                    else:
                        c.type = whole_line.split('#')[-1].rstrip().lstrip()
                    c.coeffs = idcoeffs
                    coeffs[ID] = c
    
                # Find atoms info
                if atomflag:
                    line = line.split()
                    if atomstyle == "charge":
                        # try to get image flags if they exists
                        try:
                           ID = int(line[0])
                           Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                           charge = float(line[2])
                           x = float(line[3])
                           y = float(line[4])
                           z = float(line[5])
                           imageflag_x = int(line[6])
                           imageflag_y = int(line[7])
                           imageflag_z = int(line[8])
                           if '#' not in whole_line: comment = 'N/A'
                           else: comment= whole_line.split('#')[-1].rstrip().lstrip()
                        except:
                           ID = int(line[0])
                           Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                           charge = float(line[2])
                           x = float(line[3])
                           y = float(line[4])
                           z = float(line[5])
                           imageflag_x = 0
                           imageflag_y = 0
                           imageflag_z = 0
                           if '#' not in whole_line: comment = 'N/A'
                           else: comment= whole_line.split('#')[-1].rstrip().lstrip()
                        # Save atom info
                        a = Atom()
                        a.type = Type
                        a.charge = charge
                        a.x = x
                        a.y = y
                        a.z = z
                        a.ix = imageflag_x
                        a.iy = imageflag_y
                        a.iz = imageflag_z
                        a.comment = comment
                        self.atoms[ID] = a
                    elif atomstyle == "molecular":
                        # try to get image flags if they exists
                        try:
                           ID = int(line[0])
                           Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                           charge = float(line[2])
                           x = float(line[3])
                           y = float(line[4])
                           z = float(line[5])
                           imageflag_x = int(line[6])
                           imageflag_y = int(line[7])
                           imageflag_z = int(line[8])
                           a.comment = comment
                           if '#' not in whole_line: comment = 'N/A'
                           else: comment= whole_line.split('#')[-1].rstrip().lstrip()
                        except:
                           ID = int(line[0])
                           Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                           charge = float(line[2])
                           x = float(line[3])
                           y = float(line[4])
                           z = float(line[5])
                           imageflag_x = 0
                           imageflag_y = 0
                           imageflag_z = 0 
                           if '#' not in whole_line: comment = 'N/A'
                           else: comment= whole_line.split('#')[-1].rstrip().lstrip()
                        # Save atom info
                        a = Atom()
                        a.type = type
                        a.charge = charge
                        a.x = x
                        a.y = y
                        a.z = z
                        a.ix = imageflag_x
                        a.iy = imageflag_y
                        a.iz = imageflag_z
                        a.comment = comment
                        self.atoms[ID] = a
                    elif atomstyle == "full":
                        # try to get image flags if they exists
                        try:
                           ID = int(line[0])
                           molid = int(line[1])
                           Type = label_type_mapping(line[2], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                           charge = float(line[3])
                           x = float(line[4])
                           y = float(line[5])
                           z = float(line[6])
                           imageflag_x = int(line[7])
                           imageflag_y = int(line[8])
                           imageflag_z = int(line[9])
                           if '#' not in whole_line: comment = 'N/A'
                           else: comment= whole_line.split('#')[-1].rstrip().lstrip()
                        except:
                           ID = int(line[0])
                           molid = int(line[1])
                           Type = label_type_mapping(line[2], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                           charge = float(line[3])
                           x = float(line[4])
                           y = float(line[5])
                           z = float(line[6])
                           imageflag_x = 0
                           imageflag_y = 0
                           imageflag_z = 0 
                           if '#' not in whole_line: comment = 'N/A'
                           else: comment= whole_line.split('#')[-1].rstrip().lstrip()
                        # Save atom info
                        a = Atom()
                        a.molid = molid
                        a.type = Type
                        a.charge = charge
                        a.x = x
                        a.y = y
                        a.z = z
                        a.ix = imageflag_x
                        a.iy = imageflag_y
                        a.iz = imageflag_z
                        a.comment = comment
                        self.atoms[ID] = a
                    else:
                        print(atomstyle)
                        raise Exception('Atom Style not defined as a style hint in the LAMMPS datafile.')
    
                # Find bonds
                elif bondflag:
                    line = line.split()
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.bond_type_labels_forward, self.bond_type_labels_reverse, method, section='topology')
                    atom1id = int(line[2])
                    atom2id = int(line[3])
                    b = Bond()
                    b.type = Type
                    b.atomids = [atom1id, atom2id]
                    self.bonds[ID] = b
    
                # Find angles
                elif angleflag:
                    line = line.split()
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.angle_type_labels_forward, self.angle_type_labels_reverse, method, section='topology')
                    atom1id = int(line[2])
                    atom2id = int(line[3])
                    atom3id = int(line[4])
                    c = Angle()
                    c.type = Type
                    c.atomids = [atom1id, atom2id, atom3id]
                    self.angles[ID] = c
    
                # Find dihedrals
                elif dihedralflag:
                    line = line.split()
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.dihedral_type_labels_forward, self.dihedral_type_labels_reverse, method, section='topology')
                    atom1id = int(line[2])
                    atom2id = int(line[3])
                    atom3id = int(line[4])
                    atom4id = int(line[5])
                    d = Dihedral()
                    d.type = Type
                    d.atomids = [atom1id, atom2id, atom3id, atom4id]
                    self.dihedrals[ID] = d
    
                # Find impropers
                elif improperflag:
                    line = line.split()
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.improper_type_labels_forward, self.improper_type_labels_reverse, method, section='topology')
                    atom1id = int(line[2])
                    atom2id = int(line[3])
                    atom3id = int(line[4])
                    atom4id = int(line[5])
                    i = Improper()
                    i.type = Type
                    i.atomids = [atom1id, atom2id, atom3id, atom4id]
                    self.impropers[ID] = i
    
                # Find velocities
                elif velocityflag:
                    line = line.split()
                    ID = int(line[0])
                    vx = float(line[1])
                    vy = float(line[2])
                    vz = float(line[3])
                    self.velocities[ID] = (vx, vy, vz)
                    
                    
#############################
# TESTING READ FUNCTINALITY #
#############################
if __name__ == "__main__":
    standard_file = 'read_TypeLabels_TEST/detda_typed_IFF_GT_all2lmp_standard.data'
    type_labels_file = 'read_TypeLabels_TEST/detda_typed_IFF_GT_all2lmp_TypeLabels.data'
    type_labels_file = 'read_TypeLabels_TEST/detda_typed_IFF_GT_LAMMPS_TypeLabels.data'
    max_print = 5 # Set max id's to print
    
    
    #############################
    # Try reading standard file #
    #############################
    m = Molecule_File(standard_file, method='forward')
    
    # Print findings
    print('\n\n\n-----------------------------------standard datafile-----------------------------------')
    print('{:^5} {:^15} {:^15} {:^15} {:^15} {:^5} {:^5} {:^5}'.format('atomid', 'atomtype', 'x', 'y', 'z', 'ix', 'iy', 'iz'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.atoms):
        atom = m.atoms[i]
        if n < max_print:
            print('{:^5} {:^15} {:^15} {:^15} {:^15} {:^5} {:^5} {:^5}'.format(i, atom.type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz))
    print('----------------------------------------------------------------------------------------')
    print('{:^5} {:^15} {:^15} {:^15} {:^25}'.format('bondid', 'bondtype', 'id1', 'id2', 'Python Data Type'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.bonds):
        bond = m.bonds[i]
        id1, id2 = bond.atomids
        if n < max_print:
            print('{:^5} {:^15} {:^15} {:^15} {:^25}'.format(i, bond.type, id1, id2, str(type(bond.type))))
    print('----------------------------------------------------------------------------------------')
    print('{:^15} {:^50} {:^25}'.format('bondtype', 'coeffs', 'Python Data Type'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.bond_type_labels_reverse):
        if n < max_print:
            print('{:^15} {:^50} {:^25}'.format(i, m.bond_type_labels_reverse[i], str(type(i))))
    print('----------------------------------------------------------------------------------------')
    print('{:^15} {:^50} {:^25}'.format('bondtype', 'coeffs', 'Python Data Type'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.bond_coeffs):
        coeff = m.bond_coeffs[i].coeffs
        if n < max_print:
            print('{:^15} {:^50} {:^25}'.format(i, str(coeff), str(type(i))))
            
            
    ####################################################
    # Try reading type label file with forward mapping #
    ####################################################
    m = Molecule_File(type_labels_file, method='forward')
    print(m.atom_type_labels_forward)
    print(m.atom_type_labels_reverse)
    print(m.bond_type_labels_forward)
    print(m.bond_type_labels_reverse)
    
    # Print findings
    print('\n\n\n\n-------------------------datafile w/type labels forward mapping------------------------')
    print('{:^5} {:^15} {:^15} {:^15} {:^15} {:^5} {:^5} {:^5}'.format('atomid', 'atomtype', 'x', 'y', 'z', 'ix', 'iy', 'iz'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.atoms):
        atom = m.atoms[i]
        if n < max_print:
            print('{:^5} {:^15} {:^15} {:^15} {:^15} {:^5} {:^5} {:^5}'.format(i, atom.type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz))
    print('----------------------------------------------------------------------------------------')
    print('{:^5} {:^15} {:^15} {:^15} {:^25}'.format('bondid', 'bondtype', 'id1', 'id2', 'Python Data Type'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.bonds):
        bond = m.bonds[i]
        id1, id2 = bond.atomids
        if n < max_print:
            print('{:^5} {:^15} {:^15} {:^15} {:^25}'.format(i, bond.type, id1, id2, str(type(bond.type))))
    print('----------------------------------------------------------------------------------------')
    print('{:^15} {:^50} {:^25}'.format('bondtype', 'coeffs', 'Python Data Type'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.bond_type_labels_reverse):
        if n < max_print:
            print('{:^15} {:^50} {:^25}'.format(i, m.bond_type_labels_reverse[i], str(type(i))))
    print('----------------------------------------------------------------------------------------')
    print('{:^15} {:^50} {:^25}'.format('bondtype', 'coeffs', 'Python Data Type'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.bond_coeffs):
        coeff = m.bond_coeffs[i].coeffs
        if n < max_print:
            print('{:^15} {:^50} {:^25}'.format(i, str(coeff), str(type(i))))
            
            
    ####################################################
    # Try reading type label file with reverse mapping #
    ####################################################
    m = Molecule_File(type_labels_file, method='reverse')
    print(m.atom_type_labels_forward)
    print(m.atom_type_labels_reverse)
    print(m.bond_type_labels_forward)
    print(m.bond_type_labels_reverse)
    
    # Print findings
    print('\n\n\n\n-------------------------datafile w/type labels reverse mapping------------------------')
    print('{:^5} {:^15} {:^15} {:^15} {:^15} {:^5} {:^5} {:^5}'.format('atomid', 'atomtype', 'x', 'y', 'z', 'ix', 'iy', 'iz'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.atoms):
        atom = m.atoms[i]
        if n < max_print:
            print('{:^5} {:^15} {:^15} {:^15} {:^15} {:^5} {:^5} {:^5}'.format(i, atom.type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz))
    print('----------------------------------------------------------------------------------------')
    print('{:^5} {:^15} {:^15} {:^15} {:^25}'.format('bondid', 'bondtype', 'id1', 'id2', 'Python Data Type'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.bonds):
        bond = m.bonds[i]
        id1, id2 = bond.atomids
        if n < max_print:
            print('{:^5} {:^15} {:^15} {:^15} {:^25}'.format(i, bond.type, id1, id2, str(type(bond.type))))
    print('----------------------------------------------------------------------------------------')
    print('{:^15} {:^50} {:^25}'.format('bondtype', 'coeffs', 'Python Data Type'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.bond_type_labels_reverse):
        if n < max_print:
            print('{:^15} {:^50} {:^25}'.format(i, m.bond_type_labels_reverse[i], str(type(i))))
    print('----------------------------------------------------------------------------------------')
    print('{:^15} {:^50} {:^25}'.format('bondtype', 'coeffs', 'Python Data Type'))
    print('----------------------------------------------------------------------------------------')
    for n, i in enumerate(m.bond_coeffs):
        coeff = m.bond_coeffs[i].coeffs
        if n < max_print:
            print('{:^15} {:^50} {:^25}'.format(i, str(coeff), str(type(i))))

