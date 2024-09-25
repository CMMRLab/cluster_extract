# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 16th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""


###########################
# Writing lammps datafile #
###########################
def file(mm, filename, version, atom_style):

    # Writing new file with bonds information
    with open(filename,'w') as f: 
        header = '{} > cluster_extract {}'.format(mm.header, version)
        f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters
        f.write(f'{mm.natoms} atoms\n')
        f.write(f'{mm.nbonds} bonds\n\n')
        f.write(f'{mm.natomtypes} atom types\n')
        f.write(f'{mm.nbondtypes} bond types\n\n')
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
        f.write('\nBonds\n\n')
        for i in mm.bonds:
            bond = mm.bonds[i]
            f.write('{:^2} {:^2} {:^2} {:^2}\n'.format(i, bond.type, bond.atomids[0], bond.atomids[1])) 

    return