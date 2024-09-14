# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
May 19h, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
from collections import OrderedDict


##########################################################################
# Writing lammps datafile: https://docs.lammps.org/2001/data_format.html #
##########################################################################
def lmp(mm, filename, version, atom_style, bondflag=False):

    # Writing new file with bonds information
    with open(filename,'w') as f: 
        header = '{} > cluster_extract {}'.format(mm.header, version)
        f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters
        f.write(f'{mm.natoms} atoms\n')
        if bondflag: f.write(f'{mm.nbonds} bonds\n')
        f.write('\n')
        
        f.write(f'{mm.natomtypes} atom types\n')
        if bondflag: f.write(f'{mm.nbondtypes} bond types\n')
        f.write('\n')
        
        f.write(f'{mm.xbox_line}\n{mm.ybox_line}\n{mm.zbox_line}\n')
        if mm.xy != 0 or mm.xz != 0 or mm.yz != 0:
            f.write('{:>12.9f} {:^9.9f} {:^9.9f} {} {} {}\n'.format(mm.xy, mm.xz, mm.yz, 'xy', 'xz', 'yz'))
        
        # Write masses
        f.write('\nMasses\n\n')
        for i in mm.masses: 
            comment = '{:^2} {:5}'.format('#', mm.masses[i].type)
            mass = mm.masses[i].coeffs[0]
            f.write('{:^3} {:^10.5f} {:^2}\n'.format(i, mass, comment))
            
        # Write bond coeffs section
        if bondflag:
            f.write('\nBond Coeffs\n\n')
            for i in mm.bond_coeffs:
                bond = mm.bond_coeffs[i]
                comment = '{:^2} {}-{}'.format('#', bond.symbols[0], bond.symbols[1])
                f.write('{:^3} {:^2}\n'.format(i, comment))
            
        # Write atoms    
        f.write(f'\nAtoms # {atom_style}\n\n')
        for i in mm.atoms:
            atom = mm.atoms[i]
            
            # Try to get comment data
            try: comment = '{:^5} {:>3} {:>10}'.format('#', atom.comment, atom.hybridization)
            except: comment = ''
            
            # write charge atom style
            if atom_style == 'charge':
                f.write('{:^3} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.type, atom.charge, atom.x, atom.y,\
                                                                                                                atom.z, atom.ix, atom.iy, atom.iz, comment)) 
            # write molecular atom style
            elif atom_style == 'molecular':
                f.write('{:^3} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.type, atom.charge, atom.x, atom.y,\
                                                                                                                atom.z, atom.ix, atom.iy, atom.iz, comment))                          
            # write full atom style
            elif atom_style == 'full':
                f.write('{:^3} {:^2} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, atom.type, atom.charge, atom.x, atom.y,\
                                                                                                                atom.z, atom.ix, atom.iy, atom.iz, comment)) 
            # else raise and exception
            else: raise Exception('Atom Style Error - must be charge, molecular, or full')

        # Write bonds
        if bondflag:
            f.write('\nBonds\n\n')
            for i in mm.bonds:
                bond = mm.bonds[i]
                f.write('{:^2} {:^2} {:^2} {:^2}\n'.format(i, bond.type, bond.atomids[0], bond.atomids[1])) 
    return


##########################################################################################
# Writing new file mol2 with bonds information                                           #
# https://chemicbook.com/2021/02/20/mol2-file-format-explained-for-beginners-part-2.html #
##########################################################################################
def mol2(m, filename, version, restype='molIDs', remove_PBC_bonds=True, addbox=True):
    ################################################################
    # Check if bond is periodic using the minimum image convention #
    ################################################################
    # Find box dimensions to remove periodic boundary conditions
    x = m.xbox_line.split(); y = m.ybox_line.split(); z = m.zbox_line.split();
    xlo = x[0]; xhi = x[1]; ylo = y[0]; yhi = y[1]; zlo = z[0]; zhi = z[1];
    lx = float(x[1])-float(x[0]); ly = float(y[1])-float(y[0]); lz = float(z[1])-float(z[0]);
    
    # atom positons and bonds to define simulation cell in VMD
    boxatoms = [(xhi, ylo, zhi), (xlo, ylo, zhi), (xlo, yhi, zhi), (xhi, yhi, zhi),
                (xhi, ylo, zlo), (xlo, ylo, zlo), (xlo, yhi, zlo), (xhi, yhi, zlo)]
    boxbonds = [(1, 2), (2, 3), (3, 4), (4, 1), (5, 6), (6, 7),
                (7, 8), (8, 5), (4, 8), (5, 1), (3, 7), (6, 2)]
    
    # set max_x, max_y, max_z w/ minimum image convention
    max_x = lx/2; max_y = ly/2; max_z = lz/2;
    
    # Function to check bond periodicity status
    def check_bond_periodicity(m, id1, id2):
        pbc_flag = False # Intialize and update if bond is periodic
        x1 = m.atoms[id1].x; y1 = m.atoms[id1].y; z1 = m.atoms[id1].z
        x2 = m.atoms[id2].x; y2 = m.atoms[id2].y; z2 = m.atoms[id2].z
        
        # if bond fails minimum image convention it is periodic
        if abs(x2 - x1) > max_x: pbc_flag = True
        if abs(y2 - y1) > max_y: pbc_flag = True
        if abs(z2 - z1) > max_z: pbc_flag = True
        return pbc_flag
    
    # Find bonds that we want to write and find which bonding atoms are periodic to rm
    bondIDs2write = []
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        bond = tuple(sorted([id1, id2]))
        
        # Skip over if bond if for fused rings (restype='fused') and bond not in fused ring bonds
        if restype == 'fused' and bond not in m.rings.fused.bonds: continue
        
        # Find bond periodicity if attempting to remove them else set as False
        if remove_PBC_bonds: pbc_flag = check_bond_periodicity(m, id1, id2)
        else: pbc_flag = False
        
        # Log bondID to write if not pbc_flag
        if not pbc_flag: bondIDs2write.append(i)


    # Write file
    with open(filename,'w') as f: 
        natoms = m.natoms; nbonds = len(bondIDs2write);
        if addbox: natoms += len(boxatoms); nbonds += len(boxbonds)
        # Write molecule section
        f.write('@<TRIPOS>MOLECULE\n')
        f.write(f'{m.header} > lmp2SYBYLmol2 {version} w/remove_PBC_bonds={str(remove_PBC_bonds)}\n')
        f.write(f'  {natoms} {nbonds}    0    0    0\n')
        f.write('SMALL\n')
        f.write('NO_CHARGES\n')
        f.write('****\n')
        f.write('Energy = 0\n')
        
        # Write Atoms info
        f.write('\n@<TRIPOS>ATOM\n')
        m.atoms = dict(OrderedDict(sorted(m.atoms.items()))) # sort to keep IDs as close as possible to orginal
        id_map = {} # { orginal atomID : New atomID } to make IDs contiguous if not already
        molids = []
        for n, i in enumerate(m.atoms):
            atom = m.atoms[i]
            
            # Find atoms info
            element = atom.element
            x = '{:>17.4f}'.format(atom.x) # float point
            y = '{:>10.4f}'.format(atom.y) # float point
            z = '{:>10.4f}'.format(atom.z) # float point
            charge = '{:>10.4f}'.format(atom.charge)
            
            # Find ResName info
            if restype == 'molIDs':
                subst_name = '{:>4}'.format('m'+str(atom.molid)) # Will give access through [VMD ResName Coloring]
                subst_id = atom.molid # The ID number of the substructure containing the atom (int)  [VMD RESID Coloring]
            elif restype == 'ringed':
                subst_name = '{:>4}'.format('r'+str(atom.ring)) # Will give access through [VMD ResName Coloring]
                subst_id = atom.ring # The ID number of the substructure containing the atom (int)  [VMD RESID Coloring]
            elif restype == 'fused':
                subst_name = '{:>4}'.format('f'+str(atom.fusedringID)) # Will give access through [VMD ResName Coloring]
                subst_id = atom.fusedringID # The ID number of the substructure containing the atom (int)  [VMD RESID Coloring]
            else: 
                subst_name = '****' # The name of the substructure containing the atom (string)
                subst_id = atom.molid # The ID number of the substructure containing the atom (int)  [VMD RESID Coloring]

            # Add id to map and log subst_id
            id_map[i] = n+1; molids.append(subst_id);
            
            # Write atoms info
            f.write('{:>7} {} {} {} {} {} {:>7} {:>7} {}\n'.format(n+1, element, x, y, z, element, subst_id, subst_name, charge))
        # Add in box atoms if flag
        if addbox:
            boxID = max(molids) + 1; box_map = {} # { index of location : new atomID }
            for ID, box in enumerate(boxatoms):
                x, y, z = box; n+=1; element='B'; box_map[ID+1]=n+1;
                x = '{:>17.4f}'.format(float(x)); y = '{:>10.4f}'.format(float(y))
                z = '{:>10.4f}'.format(float(z)); charge = '{:>10.4f}'.format(0);
                subst_name = 'BOX'; subst_id = boxID;
                f.write('{:>7} {} {} {} {} {} {:>7} {:>7} {}\n'.format(n+1, element, x, y, z, element, subst_id, subst_name, charge))
                
            
        # Write Bonds info
        f.write('@<TRIPOS>BOND\n')
        bondIDs2write = sorted(bondIDs2write) # sort to keep IDs as close as possible to orginal
        for n, i in enumerate(bondIDs2write):
            id1, id2 = m.bonds[i].atomids
                        
            # Set bond_type as 1 for now...
            # Possible options:
            #    1 = single
            #    2 = double
            #    3 = triple
            #    am = amide
            #    ar = aromatic
            #    du = dummy
            #    un = unknown (cannot be determined from the parameter tables)
            #    nc = not connected
            bond_type = '1'
            
            # Write bonds info
            new_id1 = id_map[id1]; new_id2 = id_map[id2]
            f.write('{:>6} {:>6} {:>6} {:>6}\n'.format(n+1, new_id1, new_id2, bond_type)) 
        # Add in box bonds if flag
        if addbox:
            for id1, id2 in boxbonds:
                new_id1 = box_map[id1]; new_id2 = box_map[id2]; bond_type = 'du'; n += 1;
                f.write('{:>6} {:>6} {:>6} {:>6}\n'.format(n+1, new_id1, new_id2, bond_type)) 
    return