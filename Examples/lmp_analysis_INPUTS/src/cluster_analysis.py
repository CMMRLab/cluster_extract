# -*- coding: utf-8 -*-
"""
Evolution of molecular weight anlaysis
Original author  : Matt Radue
Revision author 1: Prathamesh Deshpande
Revision author 2: Josh Kemppainen modified this to work with reaxff codes 

Revision 1.5
May 5th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


Questions? jdkemppa@mtu.edu
"""
##############################
# Import Necessary Libraries #
##############################          
from collections import OrderedDict                                                                      
                                                                                                               

#############################
### cluster analysis code ###
#############################
class Data: pass # .size .mass .pmass .psize .formula
class cluster_analysis:                                                                                                           
    def __init__(self, m):
        self.data = {} # { cluster-id : Data object }
        self.atoms = {} # { cluster-id : {set of atoms in cluster id }}
        self.clusters = [] # { [list of atoms in cluster1], [nclusters]}
        self.formula = {} # { atom-id : formula atom belongs too }
        self.molids = {} # {atom-id : molid}
        self.nmolecules = 0 # count of molecules in system
        self.Mw = 0  # Weight-average molar mass
        self.Mn = 0  # Number-average molar mass
        self.Mz = 0  # Higher-average molar mass 
        self.Mz1 = 0 # Higher-average z+1 molar mass 
        self.RMW = 0 # Weight-averge reduced molecular weight
        nmol_initial_zeros = 10 # Set max molIDs to intialize w/zeros
        
        #----------------------------------#
        # Function to find cluster formula #
        #----------------------------------#
        def find_cluster_formula(cluster, m):
            elements = [m.atoms[atomid].element for atomid in cluster]
            base_elements = list(sorted({i for i in elements}))
            formula = ''
            for element in base_elements:
                formula += '{}{}-'.format(element, elements.count(element))
            return formula[:-1]        
        
        # Initialize clusters and molids
        for ID in m.atoms:
            self.clusters.append(set([ID])); self.molids[ID] = 1; # initialize as 1 and update later
    
        #------------------------------------------------------------------#
        # Find clusters and then sort in a manor to guarteen same ordering #
        #------------------------------------------------------------------#
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            for c1, cluster in enumerate(self.clusters):
                if id1 in cluster: break
            for c2, cluster in enumerate(self.clusters):
                if id2 in cluster: break
            if c1 != c2:
              self.clusters[c1].update(self.clusters[c2])
              del self.clusters[c2]
        self.clusters = [tuple(sorted(i)) for i in self.clusters] # Sort each individual cluster atomIDs in ascending order
        self.clusters = sorted(self.clusters, key=lambda x: x[0]) # Sort all clusters based on 1st atomID in cluster
        self.clusters = sorted(self.clusters, key=len, reverse=True) # Sort all clusters by number of atoms
        
        # Find molids after sorting clusters. molid will be count + 1 where the largest
        # cluster will be molid = 1 and smallet cluster will be molid = nclusters in system
        for molID, cluster in enumerate(self.clusters, 1):
            for ID in cluster:
                self.molids[ID] = molID
        
        # Analyze clusters
        self.natoms = [];  self.cmass = [];   
        for c in self.clusters:
            self.cmass.append(sum([m.masses[m.atoms[i].type].coeffs[0] for i in c]))   
            self.natoms.append(len(c))
                                                                                                                      
        # Summing mass and atoms                                                                                                                                                               
        self.mass_total = sum(self.cmass); self.size_total = sum(self.natoms);
                
        #---------------------------------#
        # Initialize and analyze clusters #
        #---------------------------------#
        for n in range(nmol_initial_zeros):  
            d = Data()
            d.atoms = []; d.size = 0;
            d.mass = 0; d.pmass = 0;
            d.psize = 0; d.formula = 'Initial';
            self.data[n+1] = d     

        # Save cluster info                                 
        for n, (size, mass) in enumerate(zip(self.natoms, self.cmass)):  
            # Find cluster formula and compute compute pmass and psize
            formula = find_cluster_formula(self.clusters[n], m)
            pmass = 0; psize = 0;
            if self.mass_total > 0: pmass = round(100*mass/self.mass_total, 2)
            if self.size_total > 0: psize = round(100*size/self.size_total, 2)
            
            # Save info for for log file
            d = Data()
            d.atoms = self.clusters[n]
            d.size = size
            d.mass = mass
            d.pmass = pmass
            d.psize = psize
            d.formula = formula
            self.data[n + 1] = d
            
            # atom atoms info to atoms object
            self.atoms[n + 1] = self.clusters[n]
            
            # Loop through atoms in cluster[n] and assign formula
            for i in self.clusters[n]:
                self.formula[i] = formula
                
            # Tally count of molecules
            self.nmolecules += 1
            
        #---------------------------------------------------------------------------#
        # Find Mw, Mn, Mz, Mz1, and RMW from self.info[molID].mass. Equations from: #
        #    Polymers: Chemistry and Physics of Modern Materials pg 8-9 &           #
        #    Multiscale Modeling for Virtual Manufacturing of Thermoset Composites  #
        #---------------------------------------------------------------------------#
        ni = 0 # tally of number of clusters
        ni_mi1 = 0 # tally of ni_clusters*mi_cluster
        ni_mi2 = 0 # tally of ni_clusters*mi_cluster^2
        ni_mi3 = 0 # tally of ni_clusters*mi_cluster^3
        ni_mi4 = 0 # tally of ni_clusters*mi_cluster^4
        RMW_ni_mi = 0  # tally of ni_clusters*mi_cluster for RMW calculation
        RMW_ni_mi2 = 0 # tally of ni_clusters*mi_cluster^2 for RMW calculation
        for n in self.data:
            mass = self.data[n].mass
            ni += 1
            ni_mi1 += 1*mass
            ni_mi2 += 1*mass**2
            ni_mi3 += 1*mass**3
            ni_mi4 += 1*mass**4
            if n > 1: # skip over 1st largest cluster from tally
                RMW_ni_mi += 1*mass
                RMW_ni_mi2 += 1*mass**2
        # Only compute if div by zero is not possible else leave as default zero's
        if ni > 0: self.Mn = ni_mi1/ni
        if ni_mi1 > 0: self.Mw = ni_mi2/ni_mi1
        if ni_mi2> 0: self.Mz = ni_mi3/ni_mi2
        if ni_mi3 > 0: self.Mz1 = ni_mi4/ni_mi3
        if RMW_ni_mi > 0: self.RMW = RMW_ni_mi2/RMW_ni_mi 



##############################################################
# Function to call cluster_analysis class and update m class #
##############################################################
def add_molecule_data2m(m):
    # Class cluster_analysis class
    molecules = cluster_analysis(m)
    
    # Add molecules attribute to m
    m.molecules = molecules
    
    # Loop through m.atoms and add molecule instance and update molid
    class Molecule: pass # .mass .size .formula
    for i in m.atoms:
        atom = m.atoms[i]
        
        # Find molid and update
        molid = molecules.molids[i]
        atom.molid = molid
        
        # Create a molecule instance for each atom
        M = Molecule()
        M.formula = molecules.formula[i]
        M.mass = molecules.data[molid].mass
        M.size = molecules.data[molid].size
        atom.molecule = M
    return m


##############################################################################
# Class to recreate m class with desired info, but leave out small molecules #
##############################################################################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment .element .molecule
class Bond: pass  # .type .atomids = [atom1id, atom2id]

class cluster_removal:
    def __init__(self, m, delete_atoms, method):
        # Provide attributes of original m class that are not atoms or bonds related
        # ReaxFF Specific info
        self.reaxff_flag = m.reaxff_flag
        if self.reaxff_flag: self.reaxff = m.reaxff # Only add reaxff attribute if reaxff_flag
        
        # Bond dists Specific info
        self.bonddist_flag = m.bonddist_flag
        if self.bonddist_flag: self.bonds_via_dist = m.bonds_via_dist # Only add bond_stats attribute if bonddist_flag
        
        # General info
        self.molecules = m.molecules
        self.elements = m.elements
        self.filename = m.filename
        self.header = m.header
        self.xbox_line = m.xbox_line
        self.ybox_line = m.ybox_line
        self.zbox_line = m.zbox_line
        self.xy = m.xy
        self.xz = m.xz
        self.yz = m.yz
        self.masses = m.masses
        self.bond_coeffs = m.bond_coeffs
        self.nbondtypes = m.nbondtypes
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.natoms = 0 # Initialize and update later
        self.nbonds = 0 # Initialize and update later
        self.natomtypes = m.natomtypes
        self.kept_molecules = {} # {molecule-formula : count of molecule formula}
        self.delete_atoms = delete_atoms
        self.total_system_mass = 0 # Initialize and update later
        self.total_system_size = 0 # Initialize and update later
        
        #------------------------------#
        # Find atoms the meet criteria #
        #------------------------------#
        kept_atoms = []; idcounter = 0; # kept atoms log and idcounter to set new ids
        atomid_map = {} # {old-atomid:new-atomid}
        m.atoms = dict(OrderedDict(sorted(m.atoms.items()))) # Re-order m.atoms to keep ids as close as possible
        for i in m.atoms:
            atom = m.atoms[i]
            size = atom.molecule.size
            mass = atom.molecule.mass
            
            # If method is 'extract' use >= for selection
            if method == 'extract':
                # If delete method is mass search using mass
                if delete_atoms['method'] == 'mass':
                    if mass >= delete_atoms['criteria']:
                        idcounter += 1; atomid_map[i] = idcounter; kept_atoms.append(i); 
                # If delete method is mass search using mass
                elif delete_atoms['method'] == 'size':
                    if size >= delete_atoms['criteria']:
                        idcounter += 1; atomid_map[i] = idcounter; kept_atoms.append(i);
                else: raise Exception(f"ERROR delete_atoms['method'] method not supported:  {delete_atoms['method']}")
                
            # If method is 'by-product' use < for selection
            elif method == 'by-product':
                # If delete method is mass search using mass
                if delete_atoms['method'] == 'mass':
                    if mass < delete_atoms['criteria']:
                        idcounter += 1; atomid_map[i] = idcounter; kept_atoms.append(i);
                # If delete method is mass search using mass
                elif delete_atoms['method'] == 'size':
                    if size < delete_atoms['criteria']:
                        idcounter += 1; atomid_map[i] = idcounter; kept_atoms.append(i);
                else: raise Exception(f"ERROR delete_atoms['method'] method not supported:  {delete_atoms['method']}")
                
            # else cluster_removal method not supported
            else: raise Exception(f"ERROR cluster_removal method not supported:  {method}")
            
        #-------------------------------------------#
        # Build self.atoms instance with kept atoms #
        #-------------------------------------------#
        # sort kept atoms to keep atomids as close as possible when reseting atomids
        kept_atoms = sorted(kept_atoms) # If no volatiles are removed the ids will stay the same
        for i in kept_atoms:
            self.atoms[atomid_map[i]] = m.atoms[i] # Recreate self.atoms from m.atoms and reset atomid
            
        #-------------------------------------------#
        # Build self.bonds instance with kept atoms #
        #-------------------------------------------#
        idcounter = 0;
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            # Keep bond if id1 and id2 are in kept_atoms
            if id1 in kept_atoms and id2 in kept_atoms:
                idcounter += 1; id1_mapped = atomid_map[id1];  id2_mapped = atomid_map[id2];
                
                # Recreate self.bonds from m.bonds
                B = Bond()
                B.type = m.bonds[i].type
                B.atomids = sorted([id1_mapped, id2_mapped])
                self.bonds[idcounter] = B
        
        #------------------------------------#
        # Update self.natoms and self.nbonds #
        #------------------------------------#
        self.natoms = len(self.atoms); self.nbonds = len(self.bonds);
        
        #---------------------#
        # Find kept_molecules #
        #---------------------#
        kept_molids = sorted({self.atoms[i].molid for i in self.atoms})
        kept_formulas = [self.molecules.data[i].formula for i in kept_molids]
        unqiue_formulas = sorted(set(kept_formulas), key=len) # Find all unqiue formulas
        self.nmolecules = len(kept_molids)
        for i in unqiue_formulas:
            self.kept_molecules[i] = kept_formulas.count(i)
            
        #----------------------------------------------------------#
        # Update self.total_system_mass and self.total_system_size #
        #----------------------------------------------------------#
        for i in self.atoms:
            self.total_system_mass += self.masses[self.atoms[i].type].coeffs[0]
            self.total_system_size += 1