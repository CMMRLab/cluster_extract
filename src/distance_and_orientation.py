# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
Octboer 13th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
from tqdm import tqdm
import numpy as np
import math

###################################
# Function to fit plane to points #
###################################
def findcenter(p):
    c = [0, 0, 0]
    x = []; y = []; z = [];
    for i in p:
        x.append(i[0])
        y.append(i[1])
        z.append(i[2])
    c = [np.mean(x), np.mean(y), np.mean(z)]   
    return c
def fitplane(p):
    c = findcenter(p)
    a = np.zeros((3, len(p)))    
    for i in range(0, 3):
        for j in range(0, len(p)):
            a[i][j] = p[j][i]-c[i]
    u,s,v = np.linalg.svd(a)
    normal = u[:,2]
    return list(c), list(normal)

#########################################################################################################################        
# Finding eqn of plane from normal vector and center point                                                              #
# https://stackoverflow.com/questions/3461869/plot-a-plane-based-on-a-normal-vector-and-a-point-in-matlab-or-matplotlib #
#########################################################################################################################
def eqn_of_plane(center_point, normal_vector):
    point  = np.array(center_point)
    normal = np.array(normal_vector)
    
    # a plane is a*x+b*y+c*z+d=0
    # [a,b,c] is the normal.
    a = normal[0]; b = normal[1]; c = normal[2]
    
    # Thus, we have to calculate d and we're set
    d = -point.dot(normal) 
    return [a, b, c, d]

#############################################################
# Function to find angle in x, y, and z from normal vectors #
#############################################################
def findangle(ref_vec, test_vec):    
    dp = np.dot(ref_vec, test_vec)   
    angle = math.acos(1*1*dp) # 1*1 assumes both ref_vec and test_vec are normalized  
    angle = math.degrees(angle)
    return angle

###################################################################################
# Function to compute distance when math.dist is not available. math.dist will be #
# quicker but for some reason not all python has math.dist .... compute_distance  #
# will be time optimized but still will be slower then math.dist.                 #
###################################################################################        
def compute_distance(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return math.sqrt(dx*dx + dy*dy + dz*dz)

########################################################
# Function to unwrap atoms via minimum image conention #
########################################################
def unwrap_atoms_via_minimum_image_convention(m, xyz):
    # Find box dimensions
    xline = m.xbox_line.split();
    yline = m.ybox_line.split();
    zline = m.zbox_line.split();
    
    # Lx, Ly, and Ly box dimensions
    lx = float(xline[1])-float(xline[0])
    ly = float(yline[1])-float(yline[0])
    lz = float(zline[1])-float(zline[0])
    
    # Find Cx, Cy, Cz box centers
    Cx = (float(xline[1])+float(xline[0]))/2
    Cy = (float(yline[1])+float(yline[0]))/2
    Cz = (float(zline[1])+float(zline[0]))/2
    
    # set max_x, max_y, max_z w/ minimum image convention
    max_x = lx/2; max_y = ly/2; max_z = lz/2;
    
    # Find atom closest to the box center to use as the anchor
    distances = [compute_distance(x1, y1, z1, Cx, Cy, Cz) for x1, y1, z1 in xyz]
    clostest2center = distances.index(min(distances))
    basex, basey, basez = xyz[clostest2center]
    
    # Find newly unwrapped atoms
    xyz1 = []
    for x, y, z in xyz:
        # Find vectors in each direction
        diffx = basex - x
        diffy = basey - y
        diffz = basez - z
        
        # Apply minimum image convention
        if abs(diffx) > max_x:
            x += np.sign(diffx)*lx
        if abs(diffy) > max_y:
            y += np.sign(diffy)*ly
        if abs(diffz) > max_z:
            z += np.sign(diffz)*lz
            
        # Save new data
        xyz1.append([x, y, z])
    return xyz1

#######################################################################
# Function to add dummy planes/centers to ref_molecules and molecules #
#######################################################################
def generate_dummy_ref_planes(ref_molecules, datatable, c, normal, angles, formula, ID):
    # Save info into m.ref_molecules
    r = Mol()
    r.IDs = []
    r.xyz = []
    r.atoms = []
    r.refIDs = []
    r.center = c
    r.normal = normal
    r.angles = angles
    r.center = c
    r.plane = eqn_of_plane(c, normal)
    r.formula = '{}-{}'.format(formula, 'dummy-plane')
    ref_molecules[ID] = r
    
    # Update molecules.data
    d = Mol()
    d.atoms = []
    d.nrings = 0
    d.prings = 0
    d.size = 0
    d.mass = 0
    d.cluster = []
    d.pmass = 0
    d.psize = 0
    d.formula = '{}-{}'.format(formula, 'dummy-plane')
    d.center = c
    datatable[ID] = d
    return ref_molecules, datatable


######################################################################################
# Function to "unwrap" computed perpendicular distance by subtracting the box length #
# in the most dominate direction, which will account for cross periodic interactions #
######################################################################################
def rm_periodic_distance(perpendicular_distance, normal_ref, center_ref, m):
    # Find box dimensions
    xline = m.xbox_line.split();
    yline = m.ybox_line.split();
    zline = m.zbox_line.split();
    
    # Lx, Ly, and Ly box dimensions
    lx = float(xline[1])-float(xline[0])
    ly = float(yline[1])-float(yline[0])
    lz = float(zline[1])-float(zline[0])
    
    # 0=X-dir; 1=Y-dir; 2=Z-dir
    dominant_dir = normal_ref.index(max(normal_ref))
    
    # Set Li based on dominant_dir
    if dominant_dir == 0: Li = lx
    elif dominant_dir == 1: Li = ly
    elif dominant_dir == 2: Li = lz
    else: Li = 0
    
    # Find shortest distance
    new_distance = perpendicular_distance
    if abs(abs(perpendicular_distance) - Li) < abs(perpendicular_distance):
        new_distance = abs(abs(perpendicular_distance) - Li)
    else: new_distance = perpendicular_distance
    return new_distance

##############################################################################
# Python program to find the Perpendicular(shortest)                         # 
# distance between a point and a Plane in 3 D.                               #  
# https://www.geeksforgeeks.org/distance-between-a-point-and-a-plane-in-3-d/ #
############################################################################## 
def shortest_distance(x1, y1, z1, a, b, c, d):
    d = abs((a * x1 + b * y1 + c * z1 + d))
    e = (math.sqrt(a * a + b * b + c * c))
    normal_distance = abs(d/e)
    return normal_distance
def set_shortest_distance_sign(distance, n_ref, c_ref, n_test, c_test):
    dominant_dir = n_ref.index(max(n_ref))
    if c_ref[dominant_dir] > c_test[dominant_dir]: sign = -1
    elif c_ref[dominant_dir] < c_test[dominant_dir]: sign = 1
    else: sign = 1
    new_distance = sign*distance
    return new_distance

#################################################################
# Python program to find the Angle between                      #
# two Planes in 3 D.                                            #
# https://www.geeksforgeeks.org/angle-between-two-planes-in-3d/ #
#################################################################
def planar_angle(plane1, plane2):
    a1 = plane1[0]; b1 = plane1[1];  c1 = plane1[2];
    a2 = plane2[0]; b2 = plane2[1];  c2 = plane2[2];
     
    d = ( a1 * a2 + b1 * b2 + c1 * c2 )
    e1 = math.sqrt( a1 * a1 + b1 * b1 + c1 * c1)
    e2 = math.sqrt( a2 * a2 + b2 * b2 + c2 * c2)
    d = d / (e1 * e2)    
    try: angle = math.degrees(math.acos(d)) # angle  
    except: angle = 0
    return angle


#################################################################################
# Function to find molecule distances and orientations from reference molecules #
#################################################################################
class Mol: pass # .center .normal .angles .IDs .xyz .plane .formula
def analysis(m):
    print('Distance and orientaion analysis')
    print(m.dist_orient)
    
    #-------------------------------#
    # Set datatable based on method #
    #-------------------------------#
    if m.dist_orient['method'] == 'molIDs':
        datatable = m.molecules.data
    elif m.dist_orient['method'] == 'ringIDs':
        datatable = m.rings.clusters
    elif m.dist_orient['method'] == 'fusedringIDs':
        if not m.find_rings['fused-rings']:
            raise Exception(f"ERROR distance and orientation analysis method: {m.dist_orient['method']}, but find_rings['fused-rings'] = False")
        else: datatable = m.rings.fused.data
    else: raise Exception(f"ERROR distance and orientation analysis method: {m.dist_orient['method']} not supported")
    m.datatable = datatable
    
    #-------------------------------------------------------#
    # Find, sort, and assign IDs to all unique formulas in  #
    # the system and add attribute to m.atoms[ID].formulaID #
    #-------------------------------------------------------#
    unique_formulas = {datatable[i].formula for i in datatable}
    unique_formulas.add('N/A')
    unique_formulas = sorted(unique_formulas)
    m.formulaIDs_forward = {n+1:i for n, i in enumerate(unique_formulas)}
    m.formulaIDs_reverse = {i:n+1 for n, i in enumerate(unique_formulas)}
    for ID in m.atoms:
        m.atoms[ID].formulaID = 'N/A'
    for i in datatable:
        data = datatable[i]
        for ID in data.atoms:
            m.atoms[ID].formulaID = data.formula
    
    #-------------------------------------------------------------#
    # Find Reference molIDs center, normal vector and orientation #
    #-------------------------------------------------------------#
    m.ref_molecules = {} # { molID : mol object }
    for ID in datatable:
        mol = datatable[ID]
        if ID in m.dist_orient['reference']:
            ref_IDs = sorted([i for i in mol.atoms])
            xyz = [[m.atoms[i].x, m.atoms[i].y, m.atoms[i].z] for i in ref_IDs]
            xyz = unwrap_atoms_via_minimum_image_convention(m, xyz)
            c, normal = fitplane(xyz)
            
            # Find x, y, and z angles
            xangle = findangle([1, 0, 0], normal)
            yangle = findangle([0, 1, 0], normal)
            zangle = findangle([0, 0, 1], normal)
            angles = [xangle, yangle, zangle]
            
            # Save info into m.ref_molecules
            r = Mol()
            r.xyz = xyz
            r.IDs = ref_IDs
            r.center = c
            r.normal = normal
            r.angles = angles
            r.center = c
            r.plane = eqn_of_plane(c, normal)
            r.formula = mol.formula
            m.ref_molecules[ID] = r
            datatable[ID].center = c
            
    #----------------------------------------------------------------------------------------------------#
    # Find box dimensions and Lx, Ly, and Ly and Add dummy planes/centers to molecules and ref_molecules #
    #----------------------------------------------------------------------------------------------------#
    xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
    lxx = float(xline[1])+float(xline[0]); lyy = float(yline[1])+float(yline[0]); lzz = float(zline[1])+float(zline[0]);
    dummyID = 0
    def isfloat(num):
        try: int(num); return False
        except: return True
        
    # Add "Dummy" X-plane options
    if 'x' in m.dist_orient['reference']:
        a = [0, 90, 90]; c = [lxx/2, lyy/2, lzz/2]; n = [1, 0, 0]; formula = 'x'; dummyID += 0.1;
        m.ref_molecules, datatable = generate_dummy_ref_planes(m.ref_molecules, datatable, c, n, a, formula, round(dummyID, 3))
    if 'xlo' in m.dist_orient['reference']:
        a = [0, 90, 90]; c = [float(xline[0]), lyy/2, lzz/2]; n = [1, 0, 0]; formula = 'xlo'; dummyID += 0.1;
        m.ref_molecules, datatable = generate_dummy_ref_planes(m.ref_molecules, datatable, c, n, a, formula, round(dummyID, 3))
    if 'xhi' in m.dist_orient['reference']:
        a = [0, 90, 90]; c = [float(xline[1]), lyy/2, lzz/2]; n = [1, 0, 0]; formula = 'xhi'; dummyID += 0.1;
        m.ref_molecules, datatable = generate_dummy_ref_planes(m.ref_molecules, datatable, c, n, a, formula, round(dummyID, 3))
    
    # Add "Dummy" Y-plane options
    if 'y' in m.dist_orient['reference']:
        a = [90, 0, 90]; c = [lxx/2, lyy/2, lzz/2]; n = [0, 1, 0]; formula = 'y'; dummyID += 0.1;
        m.ref_molecules, datatable = generate_dummy_ref_planes(m.ref_molecules, datatable, c, n, a, formula, round(dummyID, 3))
    if 'ylo' in m.dist_orient['reference']:
        a = [90, 0, 90]; c = [lxx/2, float(yline[0]), lzz/2]; n = [0, 1, 0]; formula = 'ylo'; dummyID += 0.1;
        m.ref_molecules, datatable = generate_dummy_ref_planes(m.ref_molecules, datatable, c, n, a, formula, round(dummyID, 3))
    if 'yhi' in m.dist_orient['reference']:
        a = [90, 0, 90]; c = [lxx/2, float(yline[1]), lzz/2]; n = [0, 1, 0]; formula = 'yhi'; dummyID += 0.1;
        m.ref_molecules, datatable = generate_dummy_ref_planes(m.ref_molecules, datatable, c, n, a, formula, round(dummyID, 3))
    
    # Add "Dummy" Z-plane options
    if 'z' in m.dist_orient['reference']:
        a = [90, 90, 0]; c = [lxx/2, lyy/2, lzz/2]; n = [0, 0, 1]; formula = 'z'; dummyID += 0.1;
        m.ref_molecules, datatable = generate_dummy_ref_planes(m.ref_molecules, datatable, c, n, a, formula, round(dummyID, 3))
    if 'zlo' in m.dist_orient['reference']:
        a = [90, 90, 0]; c = [lxx/2, lyy/2, float(zline[0])]; n = [0, 0, 1]; formula = 'zlo'; dummyID += 0.1;
        m.ref_molecules, datatable = generate_dummy_ref_planes(m.ref_molecules, datatable, c, n, a, formula, round(dummyID, 3))
    if 'zhi' in m.dist_orient['reference']:
        a = [90, 90, 0]; c = [lxx/2, lyy/2, float(zline[1])]; n = [0, 0, 1]; formula = 'zhi'; dummyID += 0.1;
        m.ref_molecules, datatable = generate_dummy_ref_planes(m.ref_molecules, datatable, c, n, a, formula, round(dummyID, 3))
        
    #---------------------------------------------------------------------#
    # Find other molecule distances and orientations from m.ref_molecules #
    #---------------------------------------------------------------------#
    m.dist_orientation_molecules = {} # { tuple(ref_molID, test_molID) : [ normal dist, [Rx, Ry, Rz] ]}
    print('\nAnalyzing molecule orientations and distances ....')
    for ID in tqdm(datatable):
        data = datatable[ID]
        #print(ID, data.formula)
        if ID not in m.ref_molecules and len(data.atoms) >= m.dist_orient['mincluster']:
            xyz = [[m.atoms[i].x, m.atoms[i].y, m.atoms[i].z] for i in data.atoms]
            xyz = unwrap_atoms_via_minimum_image_convention(m, xyz)
            center_test, normal_test = fitplane(xyz)
            a_test, b_test, c_test, d_test = eqn_of_plane(center_test, normal_test)
            datatable[ID].center = center_test
            datatable[ID].plane = [a_test, b_test, c_test, d_test]
            
            # Loop through m.ref_molecules to get plane eqn to compute distance and orientaion of molID
            tmp = {} # { ID_ref : [[perpendicular_distance, plane_angle, relative_angles]] }
            nowned = 1 # let each molecule be owned by this many reference planes (1=each molecule is owned by 1 ref plane; 2=each molecule is owned by 2 ref planes; ....)
            for ID_ref in m.ref_molecules:
                a_ref, b_ref, c_ref, d_ref = m.ref_molecules[ID_ref].plane
                ref_xangle, ref_yangle, ref_zangle = m.ref_molecules[ID_ref].angles
                normal_ref = m.ref_molecules[ID_ref].normal
                center_ref = m.ref_molecules[ID_ref].center
                ref_IDs = m.ref_molecules[ID_ref].IDs
                ref_xyz = m.ref_molecules[ID_ref].xyz
                
                # Compute perpendicular distance, then set the sign/ direction and then rm PBC distance (If it makes distance smaller)
                perpendicular_distance = shortest_distance(center_test[0], center_test[1], center_test[2], a_ref, b_ref, c_ref, d_ref)
                perpendicular_distance = set_shortest_distance_sign(perpendicular_distance, normal_ref, center_ref, normal_test, center_test)
                perpendicular_distance = rm_periodic_distance(perpendicular_distance, normal_ref, center_ref, m)
                
                # Compute the angle between the two planes true x, y, and z angles
                plane_angle = planar_angle([a_test, b_test, c_test, d_test], [a_ref, b_ref, c_ref, d_ref])
                xangle = findangle([1, 0, 0], normal_ref)
                yangle = findangle([0, 1, 0], normal_ref)
                zangle = findangle([0, 0, 1], normal_ref)
                
                # Compute relative x, y, and z angles and add to m.dist_orientation_molecules
                relative_xangle = abs(ref_xangle - xangle)
                relative_yangle = abs(ref_yangle - yangle)
                relative_zangle = abs(ref_zangle - zangle)
                relative_angles = [relative_xangle, relative_yangle, relative_zangle]
                tmp[(ID_ref, ID)] = [perpendicular_distance, plane_angle, relative_angles]
            
            # Ordering dictionary by value (distance) in ascending order
            tmp = dict(sorted(tmp.items(), key=lambda x:abs(x[1][1]) )) # [0=keys;1=values][1=index position in value lst]
            tmp = dict(sorted(tmp.items(), key=lambda x:abs(x[1][0]) )) # [0=keys;1=values][0=index position in value lst]
            
            # log as many as nlog to m.dist_orientation_molecules
            for n, i in enumerate(tmp):
                if n+1 <= nowned:
                    m.dist_orientation_molecules[i] = tmp[i]
                    
    # Ordering dictionary by value (distance) in ascending order
    m.dist_orientation_molecules = dict(sorted(m.dist_orientation_molecules.items(), key=lambda x:abs(x[1][1]) )) # [0=keys;1=values][1=index position in value lst]
    m.dist_orientation_molecules = dict(sorted(m.dist_orientation_molecules.items(), key=lambda x:abs(x[1][0]) )) # [0=keys;1=values][0=index position in value lst] 
    return m