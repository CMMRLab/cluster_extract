##############################
# Import Necessary Libraries #
##############################
import math



# Find hybridization by partial topological considerations and partial geometric considerations (Works for periodic systems as well)
def partially_topological_partially_geometric(mm):
    
    # Find set simulataion cell bounds if periodicity is needed to be removed via the minimum image convention
    xline = mm.xbox_line.split(); yline = mm.ybox_line.split(); zline = mm.zbox_line.split();
    lx = float(xline[1])-float(xline[0])
    ly = float(yline[1])-float(yline[0])
    lz = float(zline[1])-float(zline[0])
    
    # Set VSEPR_approx_angles (degrees) to find the minimized difference of avg_angle to approx_angle to find hybridization state
    VSEPR_approx_angles = {'Sp1':     180,  # Will be used to set potenial Sp1 hybridized geometries
                           'Sp2':     120,  # Will be used to set potenial Sp2 hybridized geometries
                           'Sp3':     109.5 # Will be used to set potenial Sp3 hybridized geometries
                           }
    
    # Set elements that are known to be terminating (IE only 1-bonded neighbor)
    terminal = ['H', 'F', 'Cl', 'Ca', 'Br']
    
    # Set max in plane distance to be considered in plane (angstroms)
    #max_planar_oop_distance = 1e-2
    
    
    ######################################################
    # Loop through atoms and start finding hybridization #
    ######################################################
    for i in mm.atoms:
        
        # Find general atom info
        atom = mm.atoms[i]; nb = int(atom.nb);
        ring = int(atom.ring); element = atom.element;
        neighs1 = atom.neighbor_ids[1] # extract out 1st neigh atomIDs
        neighs2 = atom.neighbor_ids[2] # extract out 2bd neigh atomIDs
        
        
        # Get local angles by searching atomID:i to be the center atom in the angle if nb >= 2: (use only 1st neighs)
        if nb >= 2:
            local_angles = find_local_angles(i, neighs1, neighs2, method='center')
        # else get local angles by searching for atomID:i as being and outside atom of the angle (use 1st and 2nd neighs)
        else:
            local_angles = find_local_angles(i, neighs1, neighs2, method='outside')
        
        
        # Find avg_angle for all angles atom may belong too, based on how local angles are found above
        avg_angle = compute_avg_angle_degrees(mm, lx, ly, lz, local_angles)
        
        
        # Initialize atom.hybridization and update if found. Also add atom.avg_angle and atom.oop_distance
        atom.hybridization = 'unknown'; atom.avg_angle = avg_angle;
        
        
        # Find hybridization start with the most sure topological methods and move to geometry when needed
        # If element is known to be terminating and nb == 1 set 'Terminal' as hybridized state
        if element in terminal:
            atom.hybridization = 'Terminal'
            
        # If either 'C' or 'N' has 4 bonded neighbors it is guaranteed to be Sp3
        elif element in ['C', 'N'] and nb >= 4: 
            atom.hybridization = 'Sp3'
            
        # if ring is greater then 0 and less then or equal to 4 it is guaranteed to be Sp3
        elif ring > 0 and ring <= 4:
            atom.hybridization = 'Sp3'
        
        # if ring is greater then or equal to 5 and nb <= 3 it is guaranteed to be Sp2
        elif ring >= 5 and nb <= 3: 
            atom.hybridization = 'Sp2'
            
        # if atom is oxygen and nb == 2 and ring == 0 it is Sp3
        elif element == 'O' and nb == 2 and ring == 0:
            atom.hybridization = 'Sp3'
            
        # if atom is oxygen and nb == 1 it is Sp2
        # elif element == 'O' and nb == 1 and avg_angle >= 110 and avg_angle <= 130:
        #     atom.hybridization = 'Sp2'
            
        # if atom is oxygen and nb == 1 it is Sp1
        elif element == 'O' and nb == 1 and avg_angle > 130:
            atom.hybridization = 'Sp1'
            
        # Moving to geometric methods which will only be used if the above logic did not find the hybridized state
        else:
            atom.hybridization = hybridization_via_minimized_diff_of_angles(avg_angle, VSEPR_approx_angles)
            
        # Debug print
        #print(i, element, ring, nb, atom.hybridization, avg_angle, local_angles, neighs1, neighs2)
    return mm



# Function to minimized difference between avg_angles and VSEPR_approx_angles dict
def hybridization_via_minimized_diff_of_angles(avg_angle, VSEPR_approx_angles):
    # Initialize hybridization as 'Unknown' and update when minimized angle is found
    hybridization = 'unknown'
    
    # log to keep track of possible hybridization to search for minimum difference later on
    diff_log = {} # { hybridization : abs(ideal_angle - avg_angle) }
    
    # Iteration through VSEPR_approx_angles and tally difference from avg to ideal
    for hybrid in VSEPR_approx_angles:
        ideal_angle = VSEPR_approx_angles[hybrid]
        diff_log[hybrid] = abs(ideal_angle - avg_angle)
    
    # Update hybridization with the minimized difference of angles
    hybridization = min(diff_log, key=diff_log.get)
    return hybridization

# Function to find local angles
def find_local_angles(atomID, neighs1, neighs2, method):
    bonds = set([]); angles = set([]);
    
    # Find bonds from atomID and neighs1
    for i in neighs1:
        bond = tuple( sorted([atomID, i]) )
        bonds.add(bond)
    
    # Find bonds from neighs1 and neighs2 if method is to search for outside atoms
    if method == 'outside':
        for i in neighs1:
            for j in neighs2:
                bond = tuple( sorted([i, j]) )
                bonds.add(bond)
            
        
    # Find angles from bonds
    for bond1 in bonds:
        for bond2 in bonds:

            # Find angle atoms set
            bond12 = set(bond1 + bond2)
            
            # If set is 3 in length it is an angle atom set
            if len(bond12) == 3:
                
                # Find outside atoms
                outside_atom1 = set(bond1) - set(bond2)
                outside_atom2 = set(bond2) - set(bond1)
                
                # Find center atom
                center_atom = bond12 - (outside_atom1 | outside_atom2)
                
                # Set atom names
                atom1 = list(outside_atom1)[0]
                atom2 = list(center_atom)[0]
                atom3 = list(outside_atom2)[0]
                angle = (atom1, atom2, atom3)

                # Add angle to angles set
                if angle not in angles and tuple(reversed(angle)) not in angles:
                    angles.add(angle)
    
    # Find desired angles via method
    desired_angles = []
    for angle in angles:
        outside1, center, outside2 = angle
        
        # Search for central or both method
        if method == 'center':
            if atomID == center:
                desired_angles.append(angle)
                
        # Search for outside or both method
        elif method == 'outside':
            if atomID == outside1 or atomID == outside2:
                desired_angles.append(angle)
                
    return desired_angles


# Function to compute average angle of all bonded atoms
def compute_avg_angle_degrees(mm, lx, ly, lz, angles):
    # Initialize avg_angle as zero and update later on
    avg_angle = 0; angles_degrees = []; # lst to append found angles in degrees too
    
    # Loop through angles and find positon vectors and dot products
    for angle in angles:
        
        # Find center atom and bonded atoms the center to shifr positions
        center = angle[1]; bonded = [angle[0], angle[2]]
        
        # Shift any periodic atoms
        shifted = shift_bonded_atoms(mm, center, bonded, lx, ly, lz)
        
        v1, v2 = get_pos_vectors(shifted, angle)
        dotproduct = sum([i*j for (i, j) in zip(v1, v2)])
        mag_v1 = math.sqrt(sum([i**2 for i in v1]))
        mag_v2 = math.sqrt(sum([i**2 for i in v2]))
        
        # magnitudes could be zero if so error will be triggered so try computing
        try:
            angle_rad = math.acos(dotproduct/(mag_v1*mag_v2))
            angle_deg = math.degrees(angle_rad)
            angles_degrees.append(angle_deg)
        except:
            pass
        
    # if angles_degrees is not empty re-compute avg_angle
    if angles_degrees:
        avg_angle = sum(angles_degrees)/len(angles_degrees)
    return avg_angle


    
# Function to get position vectors from angle
def get_pos_vectors(shifted_positions, angle):
    # Initialize v1 and v2 as zeros and update when found
    v1 = [0, 0, 0]; v2 = [0, 0, 0];
    
    # Find atomIDs from angle with ID2 as the common atom
    ID1 = angle[0]; ID2 = angle[1]; ID3 = angle[2];

    x1 = shifted_positions[ID1][0]; y1 = shifted_positions[ID1][1]; z1 = shifted_positions[ID1][2];
    x2 = shifted_positions[ID2][0]; y2 = shifted_positions[ID2][1]; z2 = shifted_positions[ID2][2];
    x3 = shifted_positions[ID3][0]; y3 = shifted_positions[ID3][1]; z3 = shifted_positions[ID3][2];
    
    # Update position vectors moving from ID2 to ID1 or ID3
    v1 = [x2-x1, y2-y1, z2-z1]; v2 = [x2-x3, y2-y3, z2-z3];
    return v1, v2
    


# Function to shift any atoms that are found to be periodic based on the minimum image convention
def shift_bonded_atoms(mm, atomID, bonded_atoms, lx, ly, lz):
    shifted_positions = {} # {atomID : [x, y, z]}
    
    # Find central atom position and add atomID to shifted_positions
    x2 = mm.atoms[atomID].x; y2 = mm.atoms[atomID].y; z2 = mm.atoms[atomID].z;
    shifted_positions[atomID] = [x2, y2, z2]

    # Loop through bonded atoms
    for i in bonded_atoms:
        
        # Find bonded atom position
        x1 = mm.atoms[i].x; y1 = mm.atoms[i].y; z1 = mm.atoms[i].z;
        
        # if atom is near PBC shift one atom by half of the box dimension
        diffx = x2 - x1
        if diffx > 0.5*lx:
            x1 = x1 + lx
        elif diffx < -0.5*lx:
            x1 = x1 - lx 
        else:
            x1 = x1
            
        diffy = y2 - y1
        if diffy > 0.5*ly:
            y1 = y1 + ly 
        elif diffy < -0.5*ly:
            y1 = y1 - ly 
        else:
            y1 = y1
            
        diffz = z2 - z1
        if diffz > 0.5*lz:
            z1 = z1 + lz 
        elif diffz < -0.5*lz:
            z1 = z1 - lz 
        else:
            z1 = z1
            
        # save newly shifted positions
        shifted_positions[i] = [x1, y1, z1]
    return shifted_positions

  


###############################################
# Function to create hybridization data table #
###############################################
class Data:
    pass #  .size .mass .pmass .psize
    
def hybridization_data(mm):
        
    """
    EXAMPLE DATA STRUCT:
    
    DATA of Sp2 C member rings:
        mm.hybridization['C']['Sp2'].size
        mm.hybridization['C']['Sp2'].mass
        mm.hybridization['C']['Sp2'].pmass
        mm.hybridization['C']['Sp2'].psize
    
    """
    
    # Intialize new attribute of mm as hybridization
    mm.hybridization = {} # { element symbol: {'Sp1': data object, 'Sp2': data object, 'Sp3': data object, 'unknown': data object } }  mass_Sp2_C = data['C']['Sp2'].mass
    
    ###############################
    # Initialization of self.data #
    # structure to hold all info  #
    ###############################
    # All categoreis of sub dictionary
    elems2skip = ['H'] # element types to skip hybridization and stick in 'all' key
    hybrid = ['Sp1', 'Sp2', 'Sp3', 'unknown', 'all']
    
    # Create data dict in dict with instances
    for i in mm.elements:
        tmp_dict = {};
        for j in hybrid:
            d = Data()
            d.size = 0
            d.mass = 0
            d.pmass = 0
            d.psize = 0
            tmp_dict[j] = d
        
        # Create dict in dict with instances as values of sub dict
        mm.hybridization[i] = tmp_dict   
        
    #####################################################
    # Build self.data object with re-hybridized results #
    #####################################################
    for i in mm.atoms:
        atom = mm.atoms[i]
        element = atom.element
        hybridization = atom.hybridization
        mass = mm.masses[mm.atoms[i].type].coeffs[0]
        
        # add element to specific hybridization
        if element not in elems2skip:
            mm.hybridization[element][hybridization].size += 1
            mm.hybridization[element][hybridization].mass += mass 
        
        # add elememt data to all
        mm.hybridization[element]['all'].size += 1
        mm.hybridization[element]['all'].mass += mass  
            
    # Find percents and update
    for i in mm.hybridization:
        for j in mm.hybridization[i]:
            # Find pmass and psize
            pmass = 0; psize = 0;
            if mm.total_system_mass > 0: pmass = round(100*mm.hybridization[i][j].mass/mm.total_system_mass, 2)
            if mm.total_system_size > 0: psize = round(100*mm.hybridization[i][j].size/mm.total_system_size, 2)
            
            # Update mm.hybridization['element']['hybridization']
            mm.hybridization[i][j].pmass = pmass
            mm.hybridization[i][j].psize = psize
    return mm