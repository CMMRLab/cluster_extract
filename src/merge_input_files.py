# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
May 18th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""


##############################
# Import Necessary Libraries #
##############################
import src.bonds_via_distance as bonds_via_distance
import src.mol2SYBYL2lmp as mol2SYBYL2lmp
import src.read_reaxff as read_reaxff
import src.read_lmp as read_lmp
import src.mol2lmp as mol2lmp
import sys


# Function for adding info to bonds for reaxff
def add_bond_info_for_reaxff(m, bonds):    
    # Find unique bond types
    unique_bondtypes = set([]) # tuple(element1, element2)
    for id1, id2 in bonds:
        element1 = m.atoms[id1].element
        element2 = m.atoms[id2].element
        bondtype = tuple(sorted([element1, element2]))
        unique_bondtypes.add(bondtype)
    
    # Set bond type ID
    unique_bondtypes = sorted(unique_bondtypes)
    bondtypes_forward = {} # { bondID : tuple(element1, element2) }
    bondtypes_reverse = {} # { tuple(element1, element2) : bondID }
    for n, i in enumerate(unique_bondtypes):
        bondtypes_forward[n+1] = i
        bondtypes_reverse[i] = n+1
    
    # Add bonds to m.bonds and update m.nbonds and m.nbondtypes
    class Bonds: pass # .type .atomids .symbols
    m.bonds = {}
    for n, bond in enumerate(bonds):
        element1 = m.atoms[bond[0]].element
        element2 = m.atoms[bond[1]].element
        bondtype = tuple(sorted([element1, element2]))
        b = Bonds()
        b.type = bondtypes_reverse[bondtype]
        b.atomids = sorted(bond)
        b.symbols = bondtype
        m.bonds[n+1] = b

    # Update m.bond_coeffs
    class Coeff_class: pass  # .type .coeffs .symbols
    m.bonds_coeffs = {}
    for i in bondtypes_forward:
        C = Coeff_class()
        C.type = i
        C.coeffs = []
        C.symbols = bondtypes_forward[i]
        m.bond_coeffs[i] = C    
        
    # Update quantities
    m.nbonds = len(m.bonds); m.nbondtypes = len(unique_bondtypes);
    return m

# Function for adding new bond info into a LAMMPS datafile w/bonds
def add_bond_symbols2LAMMPS_datafile_w_bonds(m):
    for i in m.bond_coeffs:
        coeff = m.bond_coeffs[i]
        for j in m.bonds:
            bond = m.bonds[j]
            id1, id2 = bond.atomids
            if i == bond.type:
                element1 = m.atoms[id1].element
                element2 = m.atoms[id2].element
                coeff.symbols = tuple(sorted([element1, element2]))         
    return m
    


######################################################################################################
# Function to read -> merge -> add missing data from each file that atom_typing will be able to read #
######################################################################################################
def merge(topofile, bondfile, mass_map, bondorder, maxbonded, boundary, vdw_radius_scale):
    reaxff_flag = False # will update if files are for reaxff conversion
    bonddist_flag = False # will update if files have no bonds and bonds are found via distance
    
    #----------------------------------------------#
    # Read lammps data file and add what is needed #
    #----------------------------------------------#
    if topofile.endswith('data') or topofile.endswith('dat'):
        m = read_lmp.Molecule_File(topofile, method='forward', sections=['Atoms'])
            
        # If read in a LAMMPS .data file it will not have an element 
        # attribute in m.atoms so find element symbol from mass_map
        for i in m.atoms:
            # Find element from mass_map and add to atom instance
            atom = m.atoms[i]; mass = m.masses[atom.type].coeffs[0];
            try:
                atom.element = [i for i in mass_map if mass in mass_map[i]][0];
                atom.comment = [i for i in mass_map if mass in mass_map[i]][0];
            except:
                print(f'Not all masses in {topofile} are in the mass_map dictionary. Failed for mass: {mass}')
                sys.exit();
            
        # Add element symbol to m.masses[ID].type
        for i in m.masses:
            mass = m.masses[i]
            mass.type = [i for i in mass_map if mass.coeffs[0] in mass_map[i]][0]
            
        # Check if len(m.bonds) == 0, if so it is a reaxff file and read
        # the reaxff orginal bond order file and add to m.bonds instance
        if len(m.bonds) == 0 and bondfile != 'n.u.':
            reaxff_flag = True
            
            # update bondfile name if bondfile starts with topofile build bondfile actual name
            if bondfile.startswith('topofile'):            
                base = topofile[:topofile.rfind('.')]
                ext = bondfile.split('.')[-1]
                bondfile = '{}.{}'.format(base, ext)
                print(f'Using path and topfile name from topofile to set bondfile = {bondfile}')
            
            # Generate bond_info dictionary
            bond_info = {} # {atomtype: ('element symbol', max number of bonded atoms)} 
            for i in m.masses:
                mass = m.masses[i].coeffs[0]; element = [i for i in mass_map if mass in mass_map[i]][0];
                bond_info[i] = (element, maxbonded[element])

            # Find bonding information from reaxff file and abos/nlps avgs class
            reaxff = read_reaxff.create_bonds(bondfile, bond_info, bondorder)
            m = add_bond_info_for_reaxff(m, reaxff.bonds)
            
            
            # Add avgs, abo_stats, statistics to reaxff class and append to m
            class reaxff_stats(): pass
            R = reaxff_stats()
            R.abo_stats = reaxff.abo_stats
            R.bo_stats = reaxff.statistics
            R.ntimesteps = len(reaxff.timesteps)
            R.nbonds = len(reaxff.bonds)
            R.nflaggedbonds = len(reaxff.flagged_bonds)
            R.maxbonded = maxbonded
            m.reaxff = R
            
        # else was a lammps datafile with bonds so reset bond info
        else: m = add_bond_symbols2LAMMPS_datafile_w_bonds(m)
        

    #-----------------------------------------#
    # Read other files and add what is needed #
    #-----------------------------------------#
    # Read .mol file 
    elif topofile.endswith('mol') or topofile.endswith('sdf'):
        m = mol2lmp.Molecule_File(topofile)
               
    # Read .mol2 file (VMD MDL MOL2 file)
    elif topofile.endswith('mol2'):
        m = mol2SYBYL2lmp.Molecule_File(topofile)
       
    # If topofile was a .mol or .sdf or .mol2 file it will not have m.masses nor correct atom types so find elements and set types
    if topofile.endswith('mol') or topofile.endswith('sdf') or topofile.endswith('mol2'):
        # Find all elements in system to set atomtypes
        elements = sorted({m.atoms[i].element for i in m.atoms})
        atomtypes = {i:n+1 for n, i in enumerate(elements)}

        # Set mass from for each atomtype in atomtypes and
        # build m.masses instance as read_lmp would have
        class Coeff_class: pass  # .type .coeffs = []
        m.masses = {}
        for i in atomtypes:
            c = Coeff_class()
            c.type = i # set type as element
            c.coeffs = [mass_map[i][0]] # set mass mass_map
            m.masses[atomtypes[i]] = c
        m.natomtypes = len(m.masses)
        
        # Assign atom type in m.atoms based on element
        for i in m.atoms:
            atom = m.atoms[i];
            atom.type = atomtypes[atom.element];
            atom.comment = atom.element
            
        # Add bond coeffs to file if file is .mol or .sdf or .mol2
        m = add_bond_info_for_reaxff(m, [m.bonds[i].atomids for i in m.bonds])
            
        
    #---------------------------------------------------------------------------#
    # Check if len(m.bonds) == 0, if so generate bonds via vdw distance cutoffs #
    #---------------------------------------------------------------------------#
    if len(m.bonds) == 0:
        bonddist_flag = True
        
        # Finding bonds
        bond_creation = bonds_via_distance.generate(m, boundary, vdw_radius_scale, maxbonded)
        m = add_bond_info_for_reaxff(m, bond_creation.bonds)
         
        # Add avgs, abo_stats, statistics to reaxff class and append to m
        class bonddist_stats(): pass
        b = bonddist_stats()
        b.dist_stat = bond_creation.statistics
        b.nbonds = len(bond_creation.bonds)
        b.nflaggedbonds = len(bond_creation.flagged_bonds)
        b.maxbonded = maxbonded
        b.boundary = boundary
        b.images = bond_creation.images
        b.nb_count = bond_creation.nb_count
        b.maxbond = bond_creation.maxbond
        b.bond_status = bond_creation.bond_status
        b.vdw_radius_scale = vdw_radius_scale
        m.bonds_via_dist = b
            
    # create a new instance in m called reaff_flag to keep track of typing inputs
    m.reaxff_flag = reaxff_flag
    
    # create a new instance in m called bonddist_flag to keep track of input creation
    m.bonddist_flag = bonddist_flag
    
    # find all elements in file and add to m
    m.elements = sorted({m.masses[i].type for i in m.masses})
    return m