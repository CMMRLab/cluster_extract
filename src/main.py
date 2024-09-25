# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 3.4
Septmeber 25th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

 
    *********************************************************
    * Requirements:                                         *
    *   python 3.7+                                         *
    *                                                       *   
    * Dependencies:                                         *
    *   python tqdm module:                                 *
    *    - pip3 install tqdm (if pip manager is installed)  *
    *                                                       *
    *   python numpy module:                                *
    *    - pip3 install numpy (if pip manager is installed) *
    *                                                       *
    *********************************************************
"""


##############################
# Import Necessary Libraries #
##############################
import src.spatial_distributions as spatial_distributions
import src.distance_and_orientation as distance_and_orientation
import src.compute_bonds_distances as bond_distances
import src.cluster_analysis as cluster_analysis
import src.spatial_binning as spatial_binning
import src.ring_analysis as ring_analysis
import src.hybridization as hybridization
import src.coordination as coordination
import src.write_files as write_files
import src.write_log as write_log
import src.atom_info as atom_info
import src.plotting as plotting
import os

######################################################
### Main function to perform all atom-typing tasks ###
######################################################
def main(m, delete_atoms, compute_bond_dists, find_rings, atom_style, files2write, parent_directory,
         mass_initial, N0, dist_orient, topofile, plot_options, spatial, binning, pflag=True):
    
    
    ###########################################
    # Initialize some preliminary information #
    ###########################################
    # set version and print starting information to screen
    version = 'v3.4 / 25 September 2024'
    if pflag: print(f'\n\nRunning cluster_extract {version}')
    
    
    
    #---------------------------------------------------------------------------------#
    # Find clusters and add cluster info to m. Will add the following instances to m: #
    #   m.molecules -> .data .atoms .clusters .formula .molids                        # 
    #   m.atoms[ID] -> .molecule -> .mass .size .formula                              #
    #---------------------------------------------------------------------------------#
    m = cluster_analysis.add_molecule_data2m(m)
    
    
    #-------------------------------------------------------------------------------------------------#
    # Extract clusters and rebuild m class as mm with all necessary information and contigous atomids #
    # with for .atoms and .bonds, transfering box info, header, filename, reaxff info, molecule info  #
    # .... method='extract' to find large clusters, method='by-product' to find by-products           #
    #-------------------------------------------------------------------------------------------------#
    mm = cluster_analysis.cluster_removal(m, delete_atoms, method='extract') # mm=large kept molecules
    bp = cluster_analysis.cluster_removal(m, delete_atoms, method='by-product') # bp=removed by-products
    
    
    #---------------------------------------------------------------------------------------------------#
    # FROM HERE ON OUT ONLY ANALYZING mm CLASS SINCE IT CONTAINS THE ATOMS WERE CARE ABOUT              #
    # Find and add neighbor_ids information to mm.atoms as mm.atoms[ID].neighbor_ids; were neighbor_ids #
    # is a dictionary as:                                                                               #
    #    {1st-neighs,   2nd-neighs,   Nth-neighs,}                                                      #
    #    {1: [1, 2, 3], 2: [4, 5, 6], ....}                                                             #
    #    Use a modified BFS algorithm to find Nth neighbors away from each atom to find neighboring IDs #
    #    set the Nth_neigh_depth as 4 and find up to 4 neighs away from atomid (adjust as needed)       #
    #---------------------------------------------------------------------------------------------------#
    mm = coordination.neighbors(mm, Nth_neigh_depth=4)
    
    
    #------------------------------------------------------------------------------------------------------#
    # Find and add ring data to mm class. ring data will be found via information in find_rings dictionary #
    # below; where 'elements2walk' sets the element types that can belong in the ring and 'rings2check are #
    # sizes of rings to look for. Code provided by Jake Gissinger with modifications and additions by Josh #
    # Then add .rings .ring .ringID .ringformula to mm.atoms[ID]                                           # 
    # .rings = [lst or ringsizes atom belongs too]                                                         #
    # .ring = integer value of partioned ring                                                              #
    # .ringID = ringID (zero if not logical)                                                               #
    # .ringformula = ring formula (blank if not logical)                                                   #
    #------------------------------------------------------------------------------------------------------#
    mm = ring_analysis.add_ring_data2m(mm, find_rings)


    
    #--------------------------------------------------------------------------------------------------------#
    # Find final atom info needed to perform atom-typing and add tom m.atoms[ID] instance for each atom.     #
    # Full verbose set of instances that mm.atoms[ID] will contain for determing atom type:                  #
    #     .molecule.formula = formula atomID belongs to (sorted by element with '-' delimiter. I.E. C1-H4)   #
    #     .element = element type of atomID                                                                  #
    #     .neighbor_ids = {1: [1, 2, 3], 2: [4, 5, 6], ....} dict or neighbors at certain depths             #  
    #     .nb = number of bonded 1st neighbors                                                               #
    #     .info = [element, ringsize, nb]                                                                    #
    #     .neighbor_info = {1: [], 2: []}; lst = [['C', 6, 3], ['H', 0, 1]]; sorted by nb -> ring -> element #
    #--------------------------------------------------------------------------------------------------------#
    mm = atom_info.add(mm)


    
    #------------------------------------------------------------------------------------------------------------------#
    # Find hybrdization and add to mm.atoms[ID].hybridization. hybridization will be intialized as 'Unknown' and if    #
    # found will be set as 'Sp0' or 'Sp1' or 'Sp2' or 'Sp3'. Where the meanings are as follows:                        #
    #    'Sp1' is a linear geometry set by searching for angles around 180 degrees                                     # 
    #    'Sp2' is a Trigonal Planar geometry set in the order 1.) ring size 2.) search for angles around 120 degrees   #
    #    'Sp3' is a Tetrahedral geometry set in the order of 1.) C/N nb==4 2.) search for angles around 109.5          #
    #    'Terminal' is an atom that is terminating and does not fall into the Sp1/Sp2/Sp3 categories domain (H, F, ..) #
    #'   'unknown' is the intialized hybridization state and is the default if the hybridization could not be found    #
    #                                                                                                                  #
    # Also add mm.atoms[ID].avg_angle to atoms dict, since this maybe useful outside of the hybridization              #
    # characterization. If the atom is not in the center of an angle this will be set to zero.                         #
    #------------------------------------------------------------------------------------------------------------------#
    mm = hybridization.partially_topological_partially_geometric(mm) # Find per atom hybridization
    mm = hybridization.hybridization_data(mm) # Tally and find global hybridization data table
    
    
    #---------------------------------#
    # Find distances and orientraions #
    #---------------------------------#
    if dist_orient['perform']:
        mm.dist_orient = dist_orient
        mm = distance_and_orientation.analysis(mm)
        
    

    
    
    #------------------------------------------------------------------#
    # Find bond distances for kept cluster and volatiles if user wants #
    #------------------------------------------------------------------#
    if compute_bond_dists['kept']: mm.bonds_lengths_kept = bond_distances.compute(mm, compute_bond_dists, struct='kept')
    if compute_bond_dists['removed']: mm.bonds_lengths_removed = bond_distances.compute(bp, compute_bond_dists, struct='removed')
    
    
    #-------------------------------------------------------------------------------#
    # Compute some global quantities of the system after all analysis has been done #
    #-------------------------------------------------------------------------------#
    N = mm.nmolecules; mm.char_yield = 0; mm.p = 0; mm.Xn = 0; # Initialize w/zeros and find N
    if mass_initial > 0: mm.char_yield = 100*round(mm.total_system_mass/mass_initial, 4) # Char Yield
    if N > 0: mm.Xn = round(N0/N, 3) # Degree of Polyermization
    if N0 > 0: mm.p = round(100*((N0 - N)/N0), 3) # Conversion
    
    
    #----------------------------------------------------------------------#
    # Setting up directories and where to write final files and results to #
    #----------------------------------------------------------------------#
    pwd = os.getcwd() # Find present working directory
    path = os.path.join(pwd, parent_directory) # Find/create paths to store code results
    # If parent_directory == 'topofile' use topofile path as path
    if parent_directory == 'topofile':
        print('Using path from topofile to set parent_directory ...')
        path = os.path.dirname(os.path.abspath(topofile))
    if not os.path.isdir(path): os.mkdir(path) # Check if path1 exists. IF not create.
        
    # Change the current working directory to path
    # to write all new files to outputs directory
    os.chdir(path)

    
               
    #---------------------------------------------------------------------------#
    # Write output files based on files2write dictionary for user desired files #
    #---------------------------------------------------------------------------#
    # Find basename of file from m
    basename = os.path.basename(m.filename)
    basename = basename[:basename.rfind('.')]
        
    # Write logfile and print to screen if pflag
    logname  = '{}_clusters.txt'.format(basename)
    write_log.file(mm, bp, logname, version, compute_bond_dists, N0, dist_orient)
    if pflag: write_log.print_to_screen(logname)
    #if not files2write['cluster_extract_log-file']: os.remove(logname)
    
    # Write lmp file(s)
    if files2write['clusters_lmpfile_bonds=Y']:
        filename  = '{}_clusters_w_bonds.data'.format(basename)
        write_files.lmp(mm, filename, version, atom_style, bondflag=True)
    if files2write['clusters_lmpfile_bonds=N']:
        filename  = '{}_clusters_wo_bonds.data'.format(basename)
        write_files.lmp(mm, filename, version, atom_style, bondflag=False)
    
    # Write mol2 file(s)
    if files2write['clusters_mol2file_molIDs']:
        filename  = '{}_clusters_molIDs.mol2'.format(basename)
        write_files.mol2(mm, filename, version, restype='molIDs', remove_PBC_bonds=True, addbox=files2write['add_box_in_all_mol2files'])
    if files2write['clusters_mol2file_ringed']:
        filename  = '{}_clusters_ringed.mol2'.format(basename)
        write_files.mol2(mm, filename, version, restype='ringed', remove_PBC_bonds=True, addbox=files2write['add_box_in_all_mol2files'])
    if files2write['clusters_mol2file_fusedR'] and find_rings['fused-rings']:
        filename  = '{}_clusters_fusedR.mol2'.format(basename)
        write_files.mol2(mm, filename, version, restype='fused', remove_PBC_bonds=True, addbox=files2write['add_box_in_all_mol2files'])
    if files2write['clusters_mol2file_hybrid']:
        filename  = '{}_clusters_hybridization.mol2'.format(basename)
        write_files.mol2(mm, filename, version, restype='hybridization', remove_PBC_bonds=True, addbox=files2write['add_box_in_all_mol2files'])
    for size in files2write['custom_ring_partitioning']:
        filename  = '{}_clusters_ringed_{}.mol2'.format(basename, size)
        res = 'partitioning:{}'.format(size)
        write_files.mol2(mm, filename, version, restype=res, remove_PBC_bonds=True, addbox=files2write['add_box_in_all_mol2files'])
        
    #------------------------#
    # If dist_orient is used #
    #------------------------#
    if dist_orient['perform']:
        # Find spatial_distributions
        spatial_distributions.find(mm, basename, spatial)
        write_log.csv(mm, basename+'.csv', version) # Write csv
        
        # Plot results
        print('\nPlotting results ....')
        plotting.master(mm, basename, plot_options)
        
    #----------------------------------------------------------------------------------------------#
    # Perform spatial binning (will write files internally, so do this after changing directories) #
    #----------------------------------------------------------------------------------------------#
    if binning['perform']:
        mm = spatial_binning.analysis(mm, binning, basename)
        
    # Change back to the intial directory to keep directory free for deletion
    os.chdir(pwd)
    
    # Print file locations and completeion of code
    if pflag:
        print(f'\n\nAll outputs can be found in {path} directory')
        print('\n\nNormal program termination\n\n')
    return mm, bp