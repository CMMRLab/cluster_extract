# -*- coding: utf-8 -*-
"""
Evolution of molecular weight anlaysis
Original author  : Matt Radue
Revision author 1: Prathamesh Deshpande
Revision author 2: Josh Kemppainen modified this to work with reaxff codes 

Revision 1.4
January 4th, 2022
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
class Data:
    pass # .size .mass .pmass .psize .formula

class cluster_analysis:                                                                                                           
    def __init__(self, m):
        self.data = {} # { cluster-id : Data object }
        self.atoms = {} # { cluster-id : {set of atoms in cluster id }}
        self.formula = {} # { atom-id : formula atom belongs too }
        self.molids = {} # {atom-id : molid}
        self.clusters = set([]) # { (tuple of atoms in cluster1), (nclusters) }
        self.nmolecules = 0 # count of molecules in system
        self.Mw = 0  # Weight-average molar mass
        self.Mn = 0  # Number-average molar mass
        self.Mz = 0  # Higher-average molar mass 
        self.Mz1 = 0 # Higher-average z+1 molar mass 
        self.RMW = 0 # Weight-averge reduced molecular weight
        nmol_initial_zeros = 10 # Set max molIDs to intialize w/zeros
        
        # Function to find cluster formula
        def find_cluster_formula(cluster, m):
            elements = [m.atoms[atomid].element for atomid in cluster]
            base_elements = list(sorted({i for i in elements}))
            formula = ''
            for element in base_elements:
                formula += '{}{}-'.format(element, elements.count(element))
            return formula[:-1]
            
        # Generate graph
        graph = {i:[] for i in m.atoms}
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            graph[id1].append(id2)
            graph[id2].append(id1)
            
        # Find clusters
        print('Finding molecules ....')
        checked = {ID:False for ID in m.atoms}
        for ID in graph:
            if checked[ID]: continue
            visited=set([ID]); queue=[ID];
            while queue:
                s = queue.pop(0) 
                for neighbor in graph[s]:
                    if checked[neighbor]: continue
                    visited.add(neighbor)
                    queue.append(neighbor)
                    checked[neighbor]=True
            self.clusters.add( tuple(sorted(visited)) )
        self.clusters = sorted(self.clusters, key=lambda x: x[0])    # Sort all clusters based on 1st atomID in cluster
        self.clusters = sorted(self.clusters, key=len, reverse=True) # Sort all clusters by number of atoms
        
        # Find molids after sorting clusters. molid will be count + 1 where the largest
        # cluster will be molid = 1 and smallet cluster will be molid = nclusters in system
        self.natoms = [];  self.cmass = []; 
        for molID, cluster in enumerate(self.clusters, 1):
            self.cmass.append(sum([m.masses[m.atoms[i].type].coeffs[0] for i in cluster])); self.natoms.append(len(cluster)); 
            for ID in cluster:
                self.molids[ID] = molID

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
            # Compute pmass and psize if applicable and find cluster formula
            psize = 0; pmass = 0; # Intialize as zeros and update
            if self.size_total > 0: psize = round(100*size/self.size_total, 2)
            if self.mass_total > 0: pmass = round(100*mass/self.mass_total, 2)
            formula = find_cluster_formula(self.clusters[n], m)
            
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
            self.atoms[ n + 1 ] = self.clusters[n]
            
            # Loop through atoms in cluster[n] and assign formula
            for i in self.clusters[n]:
                self.formula[i] = formula
                
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
        atom.formula = molecules.formula[i]

    return m


##############################################################################
# Class to recreate m class with desired info, but leave out small molecules #
##############################################################################
class Atom:
    pass # .type .molid .charge .x .y .z .ix. iy .iz .comment .element .molecule

class Bond:
    pass  # .type .atomids = [atom1id, atom2id]
    
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
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.bond_coeffs = {} # { bondtype-id : coeffs class}
        self.natoms = 0 # Initialize and update later
        self.nbonds = 0 # Initialize and update later
        self.nbondtypes = 1 # No unique bond-typing will be performed
        self.natomtypes = m.natomtypes
        self.kept_molecules = {} # {molecule-formula : count of molecule formula}
        self.delete_atoms = delete_atoms
        self.total_system_mass = 0 # Initialize and update later
        self.total_system_size = 0 # Initialize and update later
        
        
        ################################
        # Find atoms the meet criteria #
        ################################
        kept_atoms = []; idcounter = 0; # kept atoms log and idcounter to set new ids
        kept = {i:False for i in m.atoms}
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
                        idcounter += 1
                        atomid_map[i] = idcounter
                        kept_atoms.append(i)
                        kept[i]=True
                        
                # If delete method is mass search using mass
                elif delete_atoms['method'] == 'size':
                    if size >= delete_atoms['criteria']:
                        idcounter += 1
                        atomid_map[i] = idcounter
                        kept_atoms.append(i)
                        kept[i]=True
                        
                # else raise exception
                else: raise Exception(f"ERROR delete_atoms['method'] method not supported:  {delete_atoms['method']}")
                
            # If method is 'by-product' use < for selection
            elif method == 'by-product':
                # If delete method is mass search using mass
                if delete_atoms['method'] == 'mass':
                    if mass < delete_atoms['criteria']:
                        idcounter += 1
                        atomid_map[i] = idcounter
                        kept_atoms.append(i)
                        kept[i]=True
                        
                # If delete method is mass search using mass
                elif delete_atoms['method'] == 'size':
                    if size < delete_atoms['criteria']:
                        idcounter += 1
                        atomid_map[i] = idcounter
                        kept_atoms.append(i)
                        kept[i]=True
                        
                # else raise exception
                else: raise Exception(f"ERROR delete_atoms['method'] method not supported:  {delete_atoms['method']}")
                
            # else cluster_removal method not supported
            else: raise Exception(f"ERROR cluster_removal method not supported:  {method}")
            
            
        #############################################
        # Build self.atoms instance with kept atoms #
        #############################################
        # sort kept atoms to keep atomids as close as possible when reseting atomids
        kept_atoms = sorted(kept_atoms) # If no volatiles are removed the ids will stay the same
        kept_molids = set()
        for i in kept_atoms:
            # Recreate self.atoms from m.atoms and reset atomid
            self.atoms[atomid_map[i]] = m.atoms[i]
            kept_molids.add(self.atoms[atomid_map[i]].molid)
            
            
            
        #############################################
        # Build self.bonds instance with kept atoms #
        #############################################
        idcounter = 0;
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            
            # Keep bond if id1 and id2 are in kept_atoms
            #if id1 in kept_atoms and id2 in kept_atoms:
            if kept[id1] and kept[id2]:
                idcounter += 1
                id1_mapped = atomid_map[id1]
                id2_mapped = atomid_map[id2]
                
                # Recreate self.bonds from m.bonds (reseting all bondtypes to 1)
                B = Bond()
                B.type = 1 # set all bondtypes as 1
                B.atomids = sorted([id1_mapped, id2_mapped])
                self.bonds[idcounter] = B
       
        
        ######################################
        # Update self.natoms and self.nbonds #
        ######################################
        self.natoms = len(self.atoms)
        self.nbonds = len(self.bonds)
        
        
        #######################
        # Find kept_molecules #
        #######################
        kept_molids = sorted(kept_molids); self.nmolecules = len(kept_molids);
        kept_formulas = [self.molecules.data[i].formula for i in kept_molids]
        unqiue_formulas = sorted(set(kept_formulas), key=len) # Find all unqiue formulas
        for i in unqiue_formulas:
            self.kept_molecules[i] = kept_formulas.count(i)
            
            
        ############################################################
        # Update self.total_system_mass and self.total_system_size #
        ############################################################
        for i in self.atoms:
            atom = self.atoms[i]
            mass = self.masses[atom.type].coeffs[0]
            self.total_system_mass += mass
            self.total_system_size += 1