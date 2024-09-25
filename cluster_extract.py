# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 3.4
Septmeber 25th, 2024
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


##############
### Inputs ###
##############
#------------------------------------------------------------------------------------------------------------#
# Files to be read in. Two system files are assigned to the following variables:                             #
#    topofile = the main file containing atom positions, element/type, and possibly bonds. Supported files:  #
#               LAMMPS .data/.dat or .mol or .sdf or .mol2 file formats.                                     #
#    bondfile = optional file that contains ReaxFF bonding information to create bonds for a ReaxFF system   #
#               to convert to a fix bond force field. If the topofile contains bonding info in it set        # 
#               bondfile = 'n.u' (short for 'not used'). Code WILL NOT attempt to read bondfile if topofile  #
#               has bonding info in it, but bondfile variable still needs to exist, hence 'n.u.' string      #
#------------------------------------------------------------------------------------------------------------#
topofile = 'Examples/auto_logging1/post_poly_rep1_105ps.data'
bondfile = 'Examples/auto_logging1/post_poly_rep1_105ps.reaxc'



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
find_rings = {'elements2walk': ['C', 'H', 'O', 'N'], # List of elements to walk (recommended to walk along all elements in system).
              'rings2check': [3, 4, 5, 6, 7],        # List of ring sizes to check for in main ring analysis code.
              'fused-rings': False,                   # True or False to run ring connectivty analysis ('perform' key must be True).
              'fused2check': [5, 6, 7]}                 # List of ring sizes to check for in fused ring analysis code.


#------------------------------------------------------------------------------------------------------------#
# Python string variable to set atom style format of the written data file. Currently supported atom styles: #
# charge, molecular, or full.                                                                                #
#------------------------------------------------------------------------------------------------------------#
atom_style = 'full'

#------------------------------------------------------------------------------------------------------------#
# Commands for creating Files that this code can write. True or False response, where True will write the    #
# file and False will skip writing of the file. Each file has a comment explaining the main contents of the  #  
# file, along with the name of the file to be written.                                                       #
#------------------------------------------------------------------------------------------------------------#
files2write = {'cluster_extract_log-file': False, # cluster extract log file of all print outs (name: *_clusters.txt)
               'clusters_lmpfile_bonds=Y': False, # datafile of extracted clusters w/ bonds (name: *_clusters_w_bonds.data)
               'clusters_lmpfile_bonds=N': False, # datafile of extracted clusters w/o bonds (name: *_clusters_wo_bonds.data)
               'clusters_mol2file_molIDs': False, # mol2file of extract clusters w/ResName as molIDs (name: *_clusters_molIDs.mol2)
               'clusters_mol2file_ringed': False, # mol2file of extract clusters w/ResName as ring-sizes (name: *_clusters_ringed.mol2)
               'clusters_mol2file_fusedR': False, # mol2file of extract clusters w/ResName as FusedRing-sizes (name: *_clusters_fusedR.mol2)
               'clusters_mol2file_hybrid': False, # mol2file of extract clusters w/ResName as FusedRing-sizes (name: *_clusters_hybridization.mol2)
               'custom_ring_partitioning': [],    # [ring sizes for preferential ring partitioning]; if empty not used (name *_clusters_ringed_SIZE.mol2)
               'add_box_in_all_mol2files': False, # add in dummy atoms in the 8-positons to make box and connect w/bonds; ResID = max(ResID)+1
               }

#------------------------------------------------------------------------------------------------------------#
# Python string variable type to set parent directory to store all new files this code will write. If the    #
# variable is left as an empty string, files will be written to the present working directory.               #
# Example to set parent_directory as with path from topofile:                                                #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile'                           #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.         #
#------------------------------------------------------------------------------------------------------------#
parent_directory = 'Examples/stand_alone'
parent_directory = 'topofile'


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
    import src.merge_input_files as merge_input_files
    from src.main import main


    # Read .data file or .mol or .mol2 or .sdf file based on extension and add necessary information
    m = merge_input_files.merge(topofile, bondfile, mass_map, bondorder, maxbonded, boundary, vdw_radius_scale)
    
    # Run main and get outputs class and byproducts class
    outputs, byproducts = main(m, delete_atoms, compute_bond_dists, find_rings, atom_style, files2write, parent_directory,
                               mass_initial, N0, dist_orient, topofile, plot_options, spatial, binning, pflag=True)
    
    
    # Examples of how to access outputs data
    def examples(find_rings, outputs):
        """
        Examples of how to access most data strcutures in cluster_extract. There will be data strcutures not
        shown in this example that you can dig around or ask Josh for help. You can also look at the write_log
        function to see how most data strcutures can be accessed as well as iterated through (mm is outputs).
        write_log function use a fair bit of data structures, but there are more deeper in the code. Most functions
        have an ending dictionary called .data were there maybe classes inside for further discretization of data.
        The .data dictionary contains most of the metrics the user is interested in.
        """
        ####################
        # neighbor example #
        ####################
        print('\n\n\n\n-----Example of accessing neighboring atoms-----')
        print('atomid 1: 1st neighboring atoms    : ', outputs.atoms[1].neighbor_ids[1])
        print('atomid 1: 1st neighboring info     : ', outputs.atoms[1].neighbor_info[1])
        
        print('\n n1=neigh1; n2=2nd neighs by n1  {n1: [2nd neighs], .... Nnumber of 1st neighs')
        print('atomid 1: 2nd neighboring atoms    : ', outputs.atoms[1].neighbor_ids)
        print('atomid 1: 2nd neighboring elements : ', outputs.atoms[1].neighbor_ids)
        
        ####################
        # molecule example #
        ####################
        # Showing how to access hybridization data.
        # DATA STRUCT: outputs.molecules.data[molID].OBJECT (key is an integer value staring from 1-to-Nmolecules)
        # OBJECTs:
        # .mass  .size  .pmass  .psize .formula
        print('\n\n-----Example of accessing molecule results -----')
        print('molID 1 mass    : ', outputs.molecules.data[1].mass)
        print('molID 1 size    : ', outputs.molecules.data[1].size)
        print('molID 2 %mass   : ', outputs.molecules.data[2].pmass)
        print('molID 2 %size   : ', outputs.molecules.data[2].psize)
        print('molID 3 formula : ', outputs.molecules.data[3].formula)
        print('Mn              : ', outputs.molecules.Mn)
        print('Mw              : ', outputs.molecules.Mw)
        
        ###############################################
        # post analysis qtys like char yield and Xn/p #
        ###############################################
        print('\n\n-----Example of accessing post analysis results -----')
        print('Char Yield               : ', outputs.char_yield)
        print('Conversion               : ', outputs.p)
        print('Degree of polymerization : ', outputs.Xn)
        
        #########################
        # hybridization example #
        #########################
        # Showing how to access hybridization data.
        # DATA STRUCT: outputs.hybridization.['elemenet']['hybridization'].OBJECT (both keys are strings
        # for atom of interest and 'element' is a string of capitol element symbol and 'hybridization' is 
        # 'Sp1' 'Sp2' 'Sp3' or 'all' for all elements).
        # OBJECTs:
        # .mass  .size  .pmass  .psize
        print('\n\n-----Example of accessing hybridization analysis-----')
        print('Mass of Sp2 C      : ', outputs.hybridization['C']['Sp2'].mass)
        print('Size of Sp2 C      : ', outputs.hybridization['C']['Sp2'].size)
        print('pmass of Sp2 C     : ', outputs.hybridization['C']['Sp2'].pmass)
        print('pmass of all Carbon: ', outputs.hybridization['C']['all'].pmass)
        
        #########################
        # Ring analysis Example #
        #########################
        # Showing how to access ringed count data.
        # DATA STRUCT: outputs.rings.data[ringsize].OBJECT (ringsize is an integer of ring of interest).
        # OBJECTs:
        # .count  .pcount
        print('\n\n-----Example of accessing ring data-----')
        print('5-member ring count  : ', outputs.rings.data[5].count)
        print('5-member ring pcount : ', outputs.rings.data[5].pcount)
        
        # Showing how to access ringed data.
        # DATA STRUCT: outputs.rings.data[ringsize].partitioned['element'].OBJECT (ringsize is an integer
        # of ring of interest and 'element' is a string of capitol element symbol or 'all' for all elements).
        # OBJECTs:
        # .mass  .size  .pmass  .psize
        print('\n\n-----Example of accessing partioned ring data-----')
        print('5-member Carbon ring atom count  : ', outputs.rings.data[5].partitioned['C'].size)
        print('5-member Carbon ring mass    : ', outputs.rings.data[5].partitioned['C'].mass)
        print('6-member Carbon ring mass    : ', outputs.rings.data[6].partitioned['C'].mass)
        print('5-member Carbon ring pmass   : ', outputs.rings.data[5].partitioned['C'].pmass)
        print('all 5-member ringed pmass    : ', outputs.rings.data[5].partitioned['all'].pmass)
        print('5-member Oxygen ring mass    : ', outputs.rings.data[5].partitioned['O'].pmass)
        print('6-member Oxygen ring mass    : ', outputs.rings.data[6].partitioned['O'].pmass)

        # Using ring analysis on clustered atom rings
        if find_rings['fused-rings']:
            # Showing how to access ringed cluster id data. 
            # DATA STRUCT: outputs.rings.clusters.data[clusterID].OBJECT (cluster ID is an interger of rank of
            # ringed cluster starting from 1). OBJECTs:
            # .mass  .size  .pmass  .psize  .nrings .rings  .formula
            print('\n\n-----Example of accessing ringed cluster data-----')
            print('Fused Ring ID 1 mass     : ', outputs.rings.fused.data[1].mass)
            print('Fused Ring ID 1 pmass    : ', outputs.rings.fused.data[1].pmass)
            print('Fused Ring ID 1 rings    : ', outputs.rings.fused.data[1].nrings)
            print('Fused Ring ID 1 %rings   : ', outputs.rings.fused.data[1].prings)
            print('Fused Ring ID 1 formula : ', outputs.rings.fused.data[1].formula)
        return
    
    # run examples function
    #examples(find_rings, outputs)