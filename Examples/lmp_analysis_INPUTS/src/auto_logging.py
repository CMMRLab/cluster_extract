# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 3.0
May 18th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import os


#########################################
# Function for writing convergence file #
#########################################
def write_convergence(logger, convergencefile):
    # Find write mode as 'w' or 'a' based on if file exists
    if os.path.isfile(convergencefile): mode = 'a'
    else: mode = 'w'
    
    # open file based on mode
    with open(convergencefile, mode) as f:
        
        # Find titles and data to write
        titles = []; data = [];
        for i in logger:
            titles.append('{}'.format(i))
            data.append('{}'.format(logger[i]))
            
        # Join with comma's
        titles = ', '.join(titles); data = ', '.join(data)

        # Write title and data if file does not exists 
        if mode == 'w':
            f.write('{}\n'.format(titles));
            f.write('{}\n'.format(data));
        else: f.write('{}\n'.format(data));
    return
            

            
########################################
# Function for writing byproducts file #
########################################           
def write_by_products(time, byproducts, byproductsfile):
    molecules = byproducts.kept_molecules
    
    # Find byproducts string to write like Jakes species/del
    byproducts_str = '{} {}'.format('Timestep', time)
    for i in molecules:
        byproducts_str += ' {} {}'.format(byproducts[i], i)
    
    # append string to file each time write_by_products is called
    with open(byproductsfile, 'a') as f:
        f.write('{}\n'.format(byproducts_str))
    return