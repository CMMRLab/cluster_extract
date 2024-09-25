# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 3.1
May 26th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


    **********************************************************
    * Requirements:                                          *
    *   python 3.7+                                          *
    *                                                        *   
    * Dependencies:                                          *
    *   python tqdm module:                                  *
    *    - pip3 install tqdm (if pip manager is installed)   *
    *                                                        *
    *   python numpy module:                                 *
    *    - pip3 install numpy (if pip manager is installed)  *
    *                                                        *
    **********************************************************
"""
##############################
# Import Necessary Libraries #
##############################
import src.merge_input_files as merge_input_files
import src.auto_logging as auto_logging
import src.rm_files as rm_files
from src.main import main
import os

# Lammps functions to import
from lammps import lammps
from mpi4py import MPI

##############
### Inputs ###
##############
#--------------------------------------------------------------------------------------------------------------#
# Python dictionary to control cluster size that will be kept. Anything cluster less than this size will not   #
# show up in the written *.data file nor the written *.nta file (atomIDs will be sorted and reset to make them #
# contiguous). If you want to use this code to keep all atoms, set this variable to zero. This will allow the  #
# code to find all clusters and all clusters. Dictionary keys/values and their meanings:                       #
# delete_atoms = {'method': METHOD-OPTION ('mass' or 'size' string)                                            #
#                 'criteria': CRITERIA-OPTION (float or int value) }                                           #
#--------------------------------------------------------------------------------------------------------------#
delete_atoms = {'method': 'mass', # 'size' or 'mass'
                'criteria': 0 }   # if 'size' natoms: if 'mass' mass in amu


#--------------------------------------------------------------------------------------------------------------#
# Initial Mass of system for calculating Char Yeild of the system. Char Yield (Cy) is computed via             #     
#    Cy = 100*[current_mass/mass_intial]; where current_mass is the mass of the system after delete_atoms      #
#                                                                                                              #
# Initial number of molecules (N0) to compute the conversion (p) and degree of polymerization (Xn):            #
#    p = 100*[(N0 - N)/N0]; where N is the number of molecules of the system after delete_atoms                #  
#    p = N0/N; where N is the number of molecules of the system after delete_atoms                             #
#--------------------------------------------------------------------------------------------------------------#                                                                                                              
mass_initial = 40334.07600000192  
N0 = 492

#-------------------------------------------------------------------------------------------------------------#
# Option to compute bond distances for kept clusters or volatiles (True or False). True will turn on analysis #
# and False will turn off analysis. For ReaxFF files with no bonds in them or .mol or .mol2 files, the bond   #
# types will be set by the elements in the bond. For LAMMPS datafiles with bonds in them the bond types will  #
# be left as is. For the 'kept' option (which will find bond distances of all kept molecules) there is an     #
# additional 'hybridization' option that may be turned on to reset the bond type based on element and         #
# hybridization, which gives a more a more realistic understanding of bond lengths based on bond type.        #
#-------------------------------------------------------------------------------------------------------------#
compute_bond_dists = {'kept': True,            # Option to find bond length of kept clusters
                      'hybridization': True,   # Option to reset bond types based on element and hybridization
                      'removed': True}         # Option to find bond length of removed clusters

#-------------------------------------------------------------------------------------------------------------#
# Option to check for rings in system and fused ring connectivity. Each option has a comment explaining what  #
# the option controls.                                                                                        #   
#-------------------------------------------------------------------------------------------------------------#
find_rings = {'elements2walk': ['C', 'H', 'O'], # List of elements to walk (recommended to walk along all elements in system).
              'rings2check': [3, 4, 5, 6, 7],   # List of ring sizes to check for in main ring analysis code.
              'fused-rings': True,              # True or False to run ring connectivty analysis ('perform' key must be True).
              'fused2check': [5, 6]}            # List of ring sizes to check for in fused ring analysis code.

#------------------------------------------------------------------------------------------------------------#
# Python string variable to set atom style format of the written data file. Currently supported atom styles: #
# charge, molecular, or full.                                                                                #
#------------------------------------------------------------------------------------------------------------#
atom_style = 'charge'

#------------------------------------------------------------------------------------------------------------#
# Commands for creating Files that this code can write. True or False response, where True will write the    #
# file and False will skip writing of the file. Each file has a comment explaining the main contents of the  #  
# file, along with the name of the file to be written.                                                       #
#------------------------------------------------------------------------------------------------------------#
files2write = {'cluster_extract_log-file': True, # cluster extract log file of all print outs (name: *_clusters.txt)
               'clusters_lmpfile_bonds=Y': True, # datafile of extracted clusters w/ bonds (name: *_clusters_w_bonds.data)
               'clusters_lmpfile_bonds=N': True, # datafile of extracted clusters w/o bonds (name: *_clusters_wo_bonds.data)
               'clusters_mol2file_molIDs': True, # mol2file of extract clusters w/ResName as molIDs (name: *_clusters_molIDs.mol2)
               'clusters_mol2file_ringed': True, # mol2file of extract clusters w/ResName as ring-sizes (name: *_clusters_ringed.mol2)
               'clusters_mol2file_fusedR': True, # mol2file of extract clusters w/ResName as FusedRing-sizes (name: *_clusters_fusedR.mol2)
               'clusters_mol2file_hybrid': True, # mol2file of extract clusters w/ResName as FusedRing-sizes (name: *_clusters_hybridization.mol2)
               'custom_ring_partitioning': [5, 6, 7],   # [ring sizes for preferential ring partitioning]; if empty not used (name *_clusters_ringed_SIZE.mol2)
               'add_box_in_all_mol2files': True, # add in dummy atoms in the 8-positons to make box and connect w/bonds; ResID = max(ResID)+1
               }

#------------------------------------------------------------------------------------------------------------#
# Python dictionary to set element symbol from the mass of the atom if reading in a LAMMPS .data file or used#
# to set the mass of an atom when reading in files that already have elemental information such as .mol,     #
# .mol2, ... Each key in the dictionary signifies an element and must have a value of a list with at least   #
# one mass in the list. For .mol or .mol2 or  ...files, the 1st mass index will be used for mass calculation #
# in the cluster analysis. No examples are given since this is self-explanatory, however if a read file does #
# cannot find the information needed from this dictionary and error will be issued and the code will exit.   #
#------------------------------------------------------------------------------------------------------------#
mass_map = {'C':   [12.0000, 12.01115, 12.01100, 12.0112, 10.01115],
            'H':   [1.0080, 1.00797, 1.00782, 1.0000],
            'O':   [15.9990, 15.99940, 15.99491, 15.994],
            'N':   [14.0000, 14.00670],
            'S':   [32.06400, 32.06],
            'F':   [18.9984, 18.998400],
            'Si':  [28.086, 28.085000, 28.085500],
             }  



#-----------------------------------------------------------------------------------------------------------------------#
# ReaxFF: Possbile bonding configurations bond order cut-offs. The bond order cut-off will be used for the time averaged#
# bond order for each bond type. This gives the ability to adjust the minimum bond order cut-off for each bond type. If #
# the bond does not exist in this bondorder dictionary, the 'unknown' key bondorder cut-off will be enforce. Add onto   #
# the dictionary as needed, the code will adjust based on the info provided in the bondorder and maxbonded dictionary.  #
#  The code will re-order the dictionary bonding pair element keys before use so only one ordering of a bond pair needs #
# to be given (IE forward and reverse ordering of C-H bond only needs either be listed as ('C','H') or ('H','C') and the#
# code will look for both forward and reverse orderings). Update as needed.                                             #
#-----------------------------------------------------------------------------------------------------------------------#
bondorder = {('C','C'):   0.3, ('C','H'):  0.3, ('C','O'):  0.3, ('C','N'):  0.3, ('C','S'):  0.3, ('C','Si'): 0.3, # C-other bonds
             ('H','H'):   0.3, ('H','O'):  0.3, ('H','N'):  0.3, ('H','S'):  0.3, ('H','Si'): 0.3,                  # H-other bonds 
             ('O','O'):   0.3, ('O','N'):  0.3, ('O','S'):  0.3, ('O','Si'): 0.3,                                   # O-other bonds
             ('N','N'):   0.3, ('N','S'):  0.3, ('N','Si'): 0.3,                                                    # N-other bonds
             ('S','S'):   0.3, ('S','Si'): 0.3,                                                                     # S-other bonds
             ('Si','Si'): 0.3,                                                                                      # Si-other bonds     
             'unknown':   0.3}                                                                                      # Unknown

#-----------------------------------------------------------------------------------------------------------------------#
# ReaxFF/Distance: maxbonded cut-off based on element type to reduce the number of bonded atoms to each element type in #
# ReaxFF systems or systems where bonds must be inferred via interatomic distances. Update as needed.                   #
#-----------------------------------------------------------------------------------------------------------------------#
maxbonded = {'C':4, 'H':1, 'O':2, 'N':4, 'S':2, 'F':1, 'Si':4, 'Xe':0, 'Ne':0, 'Kr':0, 'He':0, 'D':1, 'Cl':1,
             'Ca':1, 'Br':1, 'Ar':0, 'P':5, 'Al': 6, 'Mg': 6, 'Li': 6, 'Fe': 6, 'Na': 1, 'K':0, 'Cs':0, 'Ba':0,
             'Sr':0, 'Pb':0}


#-----------------------------------------------------------------------------------------------------------------------#
# Option to find bonds via interatomic distances. If the bondfile = 'n.u' and the topofile does not have bonds in it,   #
# these settings will control how bonds are found. Bonds will be searched for via interatomic distances that are within #
# the vdw waals radius cutoff (cut-offs used can be found in /modules/bonds_via_distance.py/vdw_radius dictionary). The #
# "scale" of the vdw radius cut-off can be adjusted by the vdw_radius_scale variable (python float value to set the     #
# multiplication value of the vdw radius). The bonds will be searched via the boundary conditions supplied by the       #
# boundary variable (python string with 3-white space separated boundary flags). Boundary flags are set up like LAMMPS  #
# where 'p' is periodic and 'f' is not periodic.                                                                        # 
#-----------------------------------------------------------------------------------------------------------------------#
# Distance: vdw_radius_scale will set the float value multiplier to adjust the vdw radius for interatomic distance searching.
# A value of 1.0 sets the vdw as it is, a value of 1.1 adds 10% more distance onto the vdw radius of each element, ...
vdw_radius_scale = 1.1

# Distance: Boundary to set the type of boundary to search for bonds. Setup like LAMMPS with the following meaning:
#   p=periodic; f=non-periodic
boundary = 'p p p'

#------------------------------------------------------------------------------------------------------------#
# Auto logging / LAMMPS analysis / LAMMPS deleting specific settings. Such as log file name and intervals ...#
#------------------------------------------------------------------------------------------------------------#
logfile = 'cluster_extract_lmp_analysis_test' # cluster extracts automated logging name
time = 0        # will be used to tally time
interval = 10  # is the time interval that gets tallied
maxtime = 400 # Sets the maximum time of files


####################################################################################
### Below are special options that "extend" cluster_extract for various projects ###
####################################################################################
##################################
# Belongs to: dist_orient method #
###############################################################################################################
# Option to for performing distances and orientation analysis on a grouping of atoms.                         # 
#    'x'   # will create a "dummy" plane with the normal orientation along the X-axis centered in the box     #
#    'y'   # will create a "dummy" plane with the normal orientation along the Y-axis centered in the box     #
#    'z'   # will create a "dummy" plane with the normal orientation along the Z-axis centered in the box     #
#    'xlo' # will create a "dummy" plane with the normal orientation along the X-axis on the lo X-plane       #
#    'xhi' # will create a "dummy" plane with the normal orientation along the X-axis on the hi X-plane       #
#    'ylo' # will create a "dummy" plane with the normal orientation along the X-axis on the lo Y-plane       #
#    'yhi' # will create a "dummy" plane with the normal orientation along the X-axis on the hi Y-plane       #
#    'zlo' # will create a "dummy" plane with the normal orientation along the X-axis on the lo Z-plane       #
#    'zhi' # will create a "dummy" plane with the normal orientation along the X-axis on the hi Z-plane       #
###############################################################################################################
dist_orient = {'perform': False,             # True or False to use this analysis
               'method': 'ringIDs',          # Grouping method: 'molIDs' or 'ringIDs' or 'fusedringIDs'
               'reference': ['x'],           # Reference IDs
               'mincluster': 5,              # Minimum cluster size to analyze
              }

##################################
# Belongs to: dist_orient method #
###############################################################################################################
# Plotting options to adjust as desired. More options can be added as desired, but for now these are what are #
# directly available to the user from the run script. Manually update can be performed on plotting script:    #
# src/molecule_distributions/plotting.py                                                                      #
###############################################################################################################
colors = ['black', 'tab:blue', 'tab:green', 'tab:red', 'tab:grey', 'tab:cyan'] # Putting up here for cleaner plot_options dict
#colors = [] # comment/uncomment to empty list
formula2name = {'C21-H24-O2':'DABPA', 'C15-H10-N2-O4':'Tolyl', 'C21-H14-N2-O4':'BMPM', 'C1902-H3760':'CNT'} # { formula : name } Putting up here for cleaner plot_options dict
#formula2name = {} # comment/uncomment to empty dict
plot_options = {'min-molecule-count': 10,     # sets the minimum count of molecules w/ a specific formulas to plot (mainly for ReaxFF)
                'dist-bin-span': 10,          # sets the bin span for distance bin plot (MUST BE int)
                'angle-bin-span': 5,          # sets the bin span for angle bin plot (can be int or float)
                'max-angle-bins': 6,          # sets the number of angle bins on the x-axis (MUST BE int)
                'ddf-bin-span': 0.25,         # sets the bin span for ddf (d=distance;d=distribution;f=function) a psuedo RDF (can be int or float)
                'ddf-mv-avg': 25,             # sets moving averge window dimension for ddf. (MUST BE int) 1=no moving average; 2=average 2-data points together; ...
                'ddf-max': 'all',             # sets ddf max plotting distance. Can be float/int value (IE 10 or 20 or ...) or string value 'all', which will plot entire dataset
                'colors': colors,             # list of colors to specify. IF empty list, colors will be set by standard color cycle
                'formula2name': formula2name, # dict for mapping formula to formula molecule name. IF empty formula will be used in legends
                'dist-range4angles': [0, 12]  # list for setting distance range to be plotted in binned orientations. IF empty the entire distance range will be used
                }

##################################
# Belongs to: dist_orient method #
############################################################################################################
# SPATIAL DISTRIBUITONS OPTIONS. NOT COMMENT SUPER WELL YET                                                #
############################################################################################################
spatial = {'bin-span': 1,                # Sets bin size (MUST BE INT)
           'direction': 'z',             # Sets direction for spatial direction ('x', 'y' , or 'z')
           'colors': colors,             # list of colors to specify. IF empty list, colors will be set by standard color cycle
           'formula2name': formula2name, # dict for mapping formula to formula molecule name. IF empty formula will be used in legends
           'normalize-atom-count': {'C1902-H3760':3}, # To divide count of atoms for this formula by 'div' to take out count of pi-electons
           }


########################################################
# Stand alone method for Spatial binning of quantities #
###############################################################################################################
# Spatial binning options. Most are documented below except the 'logging' key/value pair. The 'logging' key/  #
# value pair sets up how to log spatially binned values. The 'logging' value is a list of strings. The strings#
# set the attributes to log and can be joined together based on the 'delimiter' value. The following strings  #
# are supported:                                                                                              #
#    'element'         element of atom                                                                        #
#    'ring'            partitioned ring size                                                                  #
#    'rings'           list of ring sizes atom belongs to (empty if not in a ring)                            #
#    'ringformula'     ring formula atom belongs to (uses partitioned ring sizes)                             #
#    'hybridization'   hybridization state of the atom                                                        #
#    'formula'         molecular formula atom belongs to                                                      #
#    'type'            LAMMPS atom TypeID                                                                     #
#    'charge'          per atom charge                                                                        #
#    'molid'           moleculeID atom belongs to                                                             #
#    and others. Ask Josh.                                                                                    #
# Attributes can be joined together by using a delimiter, set by the 'delimiter' key. For example say that    #
# binning['delimiter'] = '|', the joining delimiter is '|', where attributes can now be linked like:          #
#    'attr1|attr2|att3'                                                                                       #
# Examples of linking attributes:                                                                             #
#    'element|ring'                                                                                           #
#    'element|ring|hybridization'                                                                             #
#    'formula|rings|ringformula'                                                                              #
# Full example of binning['logging'] list:                                                                    #
#    ['element', 'element|ring', 'element|ring|hybridization']                                                #
###############################################################################################################
binning = {'perform'  : True,     # True or False to use this analysis
           'direction':'xyz',     # x or y or z or xy or yz or xyz or ...
           'increment': 3.5,      # bin span in angstroms
           'measure'  : 'pmass-local',  # 'size' or 'mass' or 'psize' or 'pmass' or 'psize-local' or 'pmass-local'
           'delimiter': '|',      # delimiter to join per-atom attributes
           'logging'  : ['element', 'element|ring', 'element|ring|hybridization', 'formula|rings|ringformula']
           }




#########################################################
### Import merge input files and main function to run ###
#########################################################
if __name__ == "__main__":      
    # Set Filenames with correct paths
    convergencefile = logfile+'.csv'
    byproductsfile = logfile+'.del'
        
    # Check for convergence file and byproducts file; if they exit remove before running code
    if os.path.isfile(convergencefile): os.remove(convergencefile)
    if os.path.isfile(byproductsfile): os.remove(byproductsfile)
    
    # Loop through files to merge and run 
    niterations = int(maxtime/interval) # Compute niterations
    lmp = lammps() # set lammps as lmp
    for n in range(niterations):
        
        # Set time in ps
        time += interval  
        
        # Names of topofile and bondfile to read and analyze (update as needed)
        myid = 'PFA_rep_4_reax_c_species_NPT_2300K_2_GPa_01ps_carb_cycle_time_' + str(time) + 'ps'
        
        # Standard LAMMPS script.
        #----------------------Variables----------------------
        lmp.command('variable        myid string {}'.format(myid))
        
        #----------------------Initialization----------------------
        lmp.command('units           real')
        lmp.command('dimension       3')
        lmp.command('boundary        p p p')
        lmp.command('atom_style      charge')
        
        #----------------------ForceField----------------------
        lmp.command('read_data       ${myid}.data')
        lmp.command('pair_style	     reax/c NULL  safezone 200 mincap 2000')
        lmp.command('pair_coeff      * * ffield_CHON_for_PAN_PBO.reax C H O')
        lmp.command('fix             charges all qeq/reax 1 0.0 10.0 1.0e-6 reax/c')
                   
        #----------------------Reset atomids-----------------
        lmp.command('reset_atom_ids  sort yes')
        
        #----------------------Settings----------------------
        lmp.command('timestep        0.1')
        lmp.command('thermo          100')
        lmp.command('thermo_style    custom step temp press etotal ke pe ebond eangle edihed evdwl density lx ly lz')
    
        #----------------------Get Bonding info----------------------
        lmp.command('fix             1 all nvt temp 2500 2500 100')
        lmp.command('fix             2 all reax/c/bonds 10 ${myid}_renumbered.reaxc') # Bond order will be averaged over 11 bond orders
        lmp.command('run             100 # 0.01 ps') 
        lmp.command('write_data      ${myid}_renumbered.data') # will overwrite original datafile with sorted atomids
        
        

        
        # Find number of processors and currently running processor
        running_processor = MPI.COMM_WORLD.Get_rank()
        total_processors = MPI.COMM_WORLD.Get_size()
        
        # Only run cluster_extract and write files on processor 0
        if running_processor == 0:
        
            # Find datafile and bondfile to give to cluster_extract
            topofile = '{}_renumbered.data'.format(myid)
            bondfile = '{}_renumbered.reaxc'.format(myid)

            # Read .data file or .mol or .mol2 or .sdf file based on extension and add necessary information
            m = merge_input_files.merge(topofile, bondfile, mass_map, bondorder, maxbonded, boundary, vdw_radius_scale)
            
            # Run main and get outputs class and byproducts class
            outputs, byproducts = main(m, delete_atoms, compute_bond_dists, find_rings, atom_style, files2write, '',
                                       mass_initial, N0, dist_orient, topofile, plot_options, spatial, binning, pflag=True)
            
            
            # Rules for logger dict. { key (column name - must be unique) : value (cluster_extract metric) }. This example using ''.format()
            # to set the decimal preceison in the writing of info. All % values in cluster_extract are rounded to two decimal places so .2f.
            logger = {
                'Iteration':                          '{}'.format(n+1),
                'Time (ps)':                          '{}'.format(time),
                'Char Yield':                         '{:.2f}'.format(outputs.char_yield),
                '%Count Carbon':                      '{:.2f}'.format(outputs.hybridization['C']['all'].psize),
                '%Count Unknown C':                   '{:.2f}'.format(outputs.hybridization['C']['unknown'].psize),
                '%Count Sp1 C':                       '{:.2f}'.format(outputs.hybridization['C']['Sp1'].psize),
                '%Count Sp2 C':                       '{:.2f}'.format(outputs.hybridization['C']['Sp2'].psize),
                '%Count Sp3 C':                       '{:.2f}'.format(outputs.hybridization['C']['Sp3'].psize),
                '%Count O':                           '{:.2f}'.format(outputs.hybridization['O']['all'].psize),
                '%Count H':                           '{:.2f}'.format(outputs.hybridization['H']['all'].psize),
                '%Count 5 Ring':                       '{:.2f}'.format(outputs.rings.data[5].partitioned['C'].psize),
                '%Count 6 Ring':                       '{:.2f}'.format(outputs.rings.data[6].partitioned['C'].psize),
                '%Count 7 Ring':                       '{:.2f}'.format(outputs.rings.data[7].partitioned['C'].psize),
                    }
        

            auto_logging.write_convergence(logger, convergencefile)
            auto_logging.write_by_products(time, byproducts, byproductsfile)
            
            # Try rm_files to delete files. rm_files.delete takes in a lst_of_filenames to loop through and try deleting.
            lst_of_filenames = [topofile, bondfile] # Will delete each specific renumbered .data and .reaxc files
            #lst_of_filenames = [] # Will delete nofiles
            rm_files.delete(lst_of_filenames)
        
        
        # Clear lammps for next iteration
        lmp.command('clear')
    # Finalize MPI rung    
    MPI.Finalize()
        

    
