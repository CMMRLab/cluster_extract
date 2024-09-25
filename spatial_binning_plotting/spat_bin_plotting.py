# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
July 8th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import matplotlib.pyplot as plt


# File to analyze
txtfile = 'graphite_composite_spatial_binning_z.txt'



##################
# Start analysis #
##################
# Function to convert strings to strings or ints or floats
def convert(string):
    try:
        output = float(string)
        if '.' not in string: output = int(output)
    except: output = string
    return output

# Function to read file
def readfile(filename):
    with open(filename,'r') as f:
        data = {} # {'column-name': column-values}
        header = {} # {'column-name': column-values}
        index2column = {} # {column:index in list}
        data_flag = False
        header_flag = False
        for line in f:
            line = line.strip()
            line = line.split()
            if 'HEADER' in line:
                header_flag = True
                data_flag = False
                continue
            elif 'DATA' in line:
                header_flag = False
                data_flag = True
                continue
            if header_flag:
                if 'ITEM:' in line:
                    index2column = {n-1:i for n, i in enumerate(line) if i != 'ITEM:'}
                    header = {index2column[i]:'' for i in index2column}
                    continue
                if len(line) == len(index2column): 
                    for n, i in enumerate(line):
                        header[index2column[n]] = convert(i)
            if data_flag:
                if 'ITEM:' in line:
                    index2column = {n-1:i for n, i in enumerate(line) if i != 'ITEM:'}
                    data = {index2column[i]:[] for i in index2column}
                    continue
                if len(line) == len(index2column): 
                    for n, i in enumerate(line):
                        data[index2column[n]].append(convert(i))
    return data, header

# Read file and show columns to plot
data, header = readfile(txtfile)
print('Available columns to plot:')
for i in data:
    print(i)
    
    
# Plot some data (could be done in a for loop as well)
# Create new plot
fig, ax = plt.subplots()
plt.plot(data['center'], data['C|6'], '.', mfc='white', ms=5, markeredgecolor='tab:blue', lw=0.01, label='Carbon in 6 member rings')
plt.plot(data['center'], data['C|3'], '.', mfc='white', ms=5, markeredgecolor='tab:red', lw=0.01, label='Carbon in 3 member rings')
plt.plot(data['center'], data['H'], '.', mfc='white', ms=5, markeredgecolor='tab:green', lw=0.01, label='Hydrogen')
plt.plot(data['center'], data['C11-H18-N2|[6]|C6'], '.', mfc='white', ms=5, markeredgecolor='tab:orange', lw=0.01, label='Aromatic carbons in\nDETDA molecule')

ax.set_xlabel('Bin center (A)', fontsize=12)
ax.set_ylabel('Values (%)', fontsize=12)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), fancybox=True, ncol=2, fontsize=12)

# otherwise the legend gets clipped off of saved image
#               width, height
fig.set_size_inches(6, 5, forward=True)
fig.tight_layout()

# Save figure
basename = txtfile[:txtfile.rfind('.')]
fig.savefig(basename+'.jpeg', dpi=800)
