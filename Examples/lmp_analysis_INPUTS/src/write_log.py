# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 16th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to write log file
def file(mm, bp, logname, version, compute_bond_dists, N0):
    
    # write file
    with open(logname, 'w') as f:
        ######################################
        # Write header info and file purpose #
        ######################################
        f.write(f'This is a saved .txt file for the print outs that appear when running cluster_extract {version}\n\n')
        
        ##########################################
        # Write the elements found in the system #
        ##########################################
        f.write('\n\nElements found in system:\n')
        for i in mm.elements:
            f.write('{} {}\n'.format('-', i))
    
        ####################################
        # If reaxff flag write reaxff info #
        ####################################
        if mm.reaxff_flag:
            # Write reaxff header
            f.write('\n\n\nReaxFF atom-typing specific section:\n')
            
            # Write intial reaxff parameters of nbonds, nflaggedbonds, and ntimesteps
            f.write('{} {} {}\n'.format('Each bonding atom ids have had', mm.reaxff.ntimesteps, 'Bond orders (BO) averaged over to find average BO for each bonding pair'))
            f.write('{} {} {}\n'.format('Total bonds created     :', mm.reaxff.nbonds, '(due to meeting specified criteria)'))
            f.write('{} {} {}\n\n\n'.format('Total bonds not created :', mm.reaxff.nflaggedbonds, '(due to not meeting specified criteria)'))
            
            # Write bond type stats table
            f.write('---------------------------Bond type bond order statistics and info---------------------------\n')
            f.write('{:^10} {:^10} {:^10} {:^15} {:^15} {:^15} {:^15}\n'.format('Bond', 'Bond', 'Bond Order', 'Bond Order', 'Bond Order', 'Bond Order','Cut-off'))
            f.write('{:^10} {:^10} {:^10} {:^15} {:^15} {:^15} {:^10}\n'.format('Type', 'Count', 'Average', 'Minimum', 'Maximum', 'Standard Deviation','used'))
            f.write('----------------------------------------------------------------------------------------------\n')
            for bond in mm.reaxff.bo_stats:
                i = mm.reaxff.bo_stats[bond]
                if i: f.write('{:^10} {:^10} {:^10} {:^15} {:^15} {:^15} {:^15}\n'.format(bond, i.count, i.avg, i.min, i.max, i.std, i.cutoff))
                     
            # Write abo type stats
            f.write('\n\n')
            f.write('-------------------------------Element abo statistics and info--------------------------------\n')
            f.write('{:<10} {:^8} {:^10} {:^15} {:^13} {:^15} {:^19}\n'.format('Element', 'Element', 'abo', 'abo', 'abo', 'abo','Max bonded'))
            f.write('{:<10} {:^8} {:^10} {:^15} {:^13} {:^15} {:^14}\n'.format('Type', 'Count', 'Average', 'Minimum', 'Maximum', 'Standard Deviation','Cut-off'))
            f.write('----------------------------------------------------------------------------------------------\n')
            for element in mm.reaxff.abo_stats:
                i = mm.reaxff.abo_stats[element]
                f.write('{:<10} {:^8} {:^10} {:^15} {:^13} {:^15} {:^19}\n'.format(element, i.count, i.avg, i.min, i.max, i.std, i.cutoff))
         
            
        ###################################################
        # If bond distance flag write bond distance stats #
        ###################################################
        if mm.bonddist_flag:
            # Write bond status and boundary used to create bonds
            f.write('\n\n{} {}\n'.format('vdw radius scale', mm.bonds_via_dist.vdw_radius_scale))
            f.write('{} {} {} {}\n'.format('boundary used:', mm.bonds_via_dist.boundary, '  nimages searched:', len(mm.bonds_via_dist.images)))
            f.write('{} {} {} {}\n'.format('non-periodic bonds found:', mm.bonds_via_dist.bond_status['non-periodic'], '    periodic bonds found:', mm.bonds_via_dist.bond_status['periodic']))

            
            # Write bond type stats table
            f.write('------------------------------------Bond type bond length statistics and info------------------------------------\n')
            f.write('{:^10} {:^10} {:^15} {:^20} {:^20} {:^20} {:^18}\n'.format('Bond', 'Bond', 'Bond Length', 'Bond Length', 'Bond Length', 'Bond Length','Cut-off'))
            f.write('{:^10} {:^10} {:^15} {:^20} {:^20} {:^20} {:^18}\n'.format('Type', 'Count', 'Average', 'Minimum', 'Maximum', 'Standard Deviation','used'))
            f.write('-----------------------------------------------------------------------------------------------------------------\n')
            for bond in mm.bonds_via_dist.dist_stat:
                i = mm.bonds_via_dist.dist_stat[bond]
                if i:
                    f.write('{:^10} {:^10} {:^15} {:^20} {:^20} {:^20} {:^18}\n'.format(bond, i.count, i.avg, i.min, i.max, i.std, i.cutoff))
                    
            # Write nb table
            columns = 'element'
            for i in range(mm.bonds_via_dist.maxbond+1):
                tmp = '{:>10} {:^5}'.format('nb', 'count')
                columns += '{:^16}: {:^3}'.format(tmp, i)
            f.write('\n\n{}\n'.format( ''.join(len(columns)*['-']) ))
            f.write('{}\n'.format(columns))
            f.write('{}\n'.format( ''.join(len(columns)*['-']) ))
            for i in mm.bonds_via_dist.nb_count:
                row = '{:<6}'.format(i)
                for count in mm.bonds_via_dist.nb_count[i]:
                    row += '{:^21}'.format(mm.bonds_via_dist.nb_count[i][count])
                f.write('{}\n'.format(row))
                
        ######################
        # Write bond lengths #
        ######################
        # Write bond type dists table for kept clusters if user wants
        if compute_bond_dists['kept']:
            f.write('\n\n\n')
            f.write('-----------------------Bond type bond length statistics for kept custers----------------------\n')
            f.write('{:<20} {:^10} {:^10} {:^15} {:^15} {:^15}\n'.format('Bond', 'Bond', 'Bond Length', 'Bond Length', 'Bond Length', 'Bond Length'))
            f.write('{:<20} {:^10} {:^10} {:^15} {:^15} {:^15}\n'.format('Type', 'Count', 'Average', 'Minimum', 'Maximum', 'Standard Deviation'))
            f.write('----------------------------------------------------------------------------------------------\n')
            for i in mm.bonds_lengths_kept.types:
                bond = mm.bonds_lengths_kept.types[i]
                symbols = '{}-{}'.format(bond.symbols[0], bond.symbols[1])
                f.write('{:<20} {:^10} {:^10} {:^15} {:^15} {:^15}\n'.format(symbols, bond.count, bond.avg, bond.min, bond.max, bond.std))
                
        # Write bond type dists table for volatiles if user wants
        if compute_bond_dists['removed']:
            f.write('\n\n\n')
            f.write('------------------------Bond type bond length statistics for volatiles------------------------\n')
            f.write('{:<20} {:^10} {:^10} {:^15} {:^15} {:^15}\n'.format('Bond', 'Bond', 'Bond Length', 'Bond Length', 'Bond Length', 'Bond Length'))
            f.write('{:<20} {:^10} {:^10} {:^15} {:^15} {:^15}\n'.format('Type', 'Count', 'Average', 'Minimum', 'Maximum', 'Standard Deviation'))
            f.write('----------------------------------------------------------------------------------------------\n')
            for i in mm.bonds_lengths_removed.types:
                bond = mm.bonds_lengths_removed.types[i]
                symbols = '{}-{}'.format(bond.symbols[0], bond.symbols[1])
                f.write('{:<20} {:^10} {:^10} {:^15} {:^15} {:^15}\n'.format(symbols, bond.count, bond.avg, bond.min, bond.max, bond.std))
         
        ####################################
        # Write molecules/cluster findings #
        ####################################
        # Write molecules table          
        f.write('\n\n\n--------------------------------------------Cluster Analysis-------------------------------------\n')
        f.write('{:^10} {:^15} {:^20} {:^15} {:^15} {:^15}\n'.format('molID', 'Molecule Size', 'Mass', '%Mass', '%Size', 'Molecule Formula'))
        f.write('-------------------------------------------------------------------------------------------------\n')  
        for i in mm.molecules.data:
            data = mm.molecules.data[i]
            size = '{: >6}'.format(data.size)
            mass = '{:.2f}'.format(data.mass)
            pmass = '{:.2f}'.format(data.pmass)
            psize = '{:.2f}'.format(data.psize)
            formula = '{:^10}'.format(data.formula)
            if data.size > 0:
                f.write('{:^10} {:^15} {:^20} {:^15} {:^15} {:^15}\n'.format(i, size, mass, pmass, psize, formula))
                
        # Write polymer chemistry findings
        f.write('\n\n*NOTE: Xn and p calculations are based on clusters that meet your minimum cluster size/mass. Mn is computed before removal of molecules*\n')
        f.write('{:^25} {}\n'.format('N0 used for calculations:', N0))
        f.write('{:^30} {}\n'.format('Number average molar mass (Mn) :', mm.molecules.Mn))
        f.write('{:^30} {}\n'.format('Degree of Polymerization (Xn)  :', mm.Xn))
        f.write('{:^30} {} {}\n'.format('Extent of reaction (conversion):', mm.p, '%' ))

            
        #####################
        # Print By products #
        #####################
        # Write delete atoms criteria
        f.write('\n\nBy-products criteria: {}\n'.format(str(bp.delete_atoms)))
        
        # Warn if delete_atoms was used and files extension was not .dat or .data
        if bp.delete_atoms['criteria'] != 0 and not mm.filename.endswith('data') or mm.filename.endswith('dat'):
            f.write('WARNING Using delete_atoms option for a Non-LAMMPS file!\n\n')
            
        # Write by-products table
        f.write('------By Products Tally------\n')
        f.write('{:<14} {:>14}\n'.format('Type', 'Count'))
        f.write('-----------------------------\n')
        for i in bp.kept_molecules:
            f.write('{:<14} {:>14}\n'.format(i, bp.kept_molecules[i]))
            
        ###########################
        # Write ringed atoms data #
        ###########################
        # Write all inputs of find_rings
        f.write('\n\n\n----Inputs used for find_rings----\n')
        f.write('{} {}\n'.format('Walked along elements  : ', mm.find_rings['elements2walk']))
        f.write('{} {}\n'.format('Checked for ring sizes : ', mm.find_rings['rings2check']))
        f.write('{} {}\n'.format('Fused ring clusters    : ', mm.find_rings['fused-rings']))
        f.write('{} {}\n'.format('Fused rings checked    : ', mm.find_rings['fused2check']))
        f.write('{} {}\n'.format('Total rings found      : ', mm.rings.total))
        
        f.write(f'{mm.rings.partitioned_count} atoms along a ring seam had to be partitioned amoung the ring\n') 
        f.write('types (To find accurate %Mass of ring type to entire system).\n')
        f.write('Giving preference of partionioning in this order:\n')
        f.write('- 6 member ring\n')
        f.write('- 5 member ring\n')
        f.write('- 7 member ring\n')
        f.write('- 4 member ring\n')
        f.write('- 3 member ring\n')
        f.write('- 8 member ring\n')
        f.write('- minimum ring size\n')
        f.write('*NOTE: If count of rings exists, but no atoms exist for that ring, this means the\n')
        f.write('atoms for that ring were partionted to other rings that the atoms belong to.*\n\n')                                                                                                     
        for i in mm.rings.data:
            data = mm.rings.data[i]
            count = '{:^5}'.format(data.count)
            pcount = '{:.2f}'.format(data.pcount)
             
            # If count is greater then zero write
            if data.count > 0:
                f.write('---------------------------------------------------------------------------------------------\n')
                f.write('|{:^25} | {:^25} | {:^35}|\n'.format('Ring', 'Count', '%Ring count'))
                f.write('|{:^25} | {:^25} | {:^35}|\n'.format('Type', 'of Rings', 'per all rings'))
                f.write('---------------------------------------------------------------------------------------------\n')
                f.write('|{:^25} | {:^25} | {:^35}|\n'.format(i, count, pcount))
                f.write('---------------------------------------------------------------------------------------------\n')
                f.write('|{:^16} | {:^15} | {:^16} | {:^16} | {:^16}|\n'.format('Element', 'natoms', 'Mass', '%Mass', '%natoms'))
                f.write('---------------------------------------------------------------------------------------------\n')
                for j in mm.rings.data[i].partitioned:
                    data1 = mm.rings.data[i].partitioned[j]
                    size = '{:^5}'.format(data1.size)
                    mass = '{:.2f}'.format(data1.mass)
                    pmass = '{:.2f}'.format(data1.pmass)
                    mass = '{:.2f}'.format(data1.mass)
                    psize = '{:.2f}'.format(data1.psize)
                    f.write('|{:^16} | {:^15} | {:^16} | {:^16} | {:^16}|\n'.format(j, size, mass, pmass, psize))
                f.write('---------------------------------------------------------------------------------------------\n\n\n')


        #################################
        # Write ringed cluster findings #
        #################################
        if len(mm.rings.clusters) > 0:
            maxID = 10 # to stop printing after maxID
            # Write ringed clusters table          
            f.write('-------------------------------------Ringed Clusters--------------------------------------\n')
            f.write('{:^10} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('ringID', 'Ring Size', 'Ring Mass', '%Mass', '%Size', 'Ring Formula'))
            f.write('------------------------------------------------------------------------------------------\n')  
            for i in mm.rings.clusters:
                data = mm.rings.clusters[i]
                size = '{: >6}'.format(data.size)
                mass = '{:.2f}'.format(data.mass)
                pmass = '{:.2f}'.format(data.pmass)
                psize = '{:.2f}'.format(data.psize)
                formula = '{:^10}'.format(data.formula)
                if i <= maxID and data.size >= 3:
                    f.write('{:^10} {:>15} {:>15} {:>15} {:>15} {:>15} {}\n'.format(i, size, mass, pmass, psize, formula, str(data.atoms)))
            # if len(mm.rings.clusters.data) > maxID print ...
            if len(mm.rings.clusters) > maxID:
                size = '{: >6}'.format('.'); mass = '{:.2}'.format('.'); pmass = '{:.2}'.format('.'); psize = '{:.2}'.format('.'); formula = '{:^10}'.format('.');
                f.write('{:^10} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('.', size, mass, pmass, psize, formula))
                f.write('{:^10} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('.', size, mass, pmass, psize, formula))
                f.write('{:^10} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('.', size, mass, pmass, psize, formula))
                
        #############################################
        # Write fused ring findings if user desired #
        #############################################
        if mm.find_rings['fused-rings']:
            # Write ringed clusters table          
            f.write('\n\n\n--------------------------------------------------Fused Ring Clusters-----------------------------------------------------\n')
            f.write('{:<6} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>35}\n'.format('FusedID', 'Size', 'Mass', '%Mass', '%Size', 'Nrings', '%Rings', 'FusedRing Formula'))
            f.write('--------------------------------------------------------------------------------------------------------------------------\n')  
            for i in mm.rings.fused.data:
                data = mm.rings.fused.data[i]
                size = '{: >6}'.format(data.size)
                mass = '{:.2f}'.format(data.mass)
                pmass = '{:.2f}'.format(data.pmass)
                psize = '{:.2f}'.format(data.psize)
                nrings = '{:>6}'.format(data.nrings)
                prings = '{:>6}'.format(data.prings)
                formula = '{:^10}'.format(data.formula)
                if data.nrings > 1:
                    f.write('{:<6} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>35}\n'.format(i, size, mass, pmass, psize, nrings, prings, formula))


        ###############################
        # Write hybridization results #
        ###############################
        f.write('\n\n\n-----------------------------Hybridization Information-------------------------------\n')
        f.write('{:^20} {:^16} {:^16} {:^16} {:^16}\n'.format('Atom-Type', 'natoms', 'Mass', '%Mass', '%natoms')) 
        f.write('-------------------------------------------------------------------------------------\n')
        for element in mm.hybridization:
            for hybridization in mm.hybridization[element]:
                atomtype = '{} {}'.format(hybridization, element)
                size = mm.hybridization[element][hybridization].size
                mass = '{:.2f}'.format(mm.hybridization[element][hybridization].mass)
                pmass = '{:.2f}'.format(mm.hybridization[element][hybridization].pmass)
                psize = '{:.2f}'.format(mm.hybridization[element][hybridization].psize)
                
                # Only write data if natoms > 0 and do not write all hybridization
                if size > 0:# and hybridization != 'all':
                    f.write('{:^20} {:^16} {:^16} {:^16} {:^16}\n'.format(atomtype, size, mass, pmass, psize))   
                    
        ######################################
        # Write carbon yield and carbon loss #
        ######################################
        f.write('\n\n----------------Char yield data----------------\n')
        f.write('{:^15} {:.2f}\n'.format('Char yield             : ', mm.char_yield))
        f.write('{:^15} {:.2f}\n'.format('Extracted Cluster Mass : ', mm.total_system_mass))
        f.write('{:^15} {:.2f}\n'.format('Entire System Mass     : ', mm.molecules.mass_total))
    return


# use this function to print the file to the screen
def print_to_screen(name):
    # open the file
    with open(name, 'r') as f:
        for n, line in enumerate(f):
            
            # Only print line if it is at least 2 lines in
            if n > 1:
                print(line.strip())