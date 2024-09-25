# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
May 4th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import matplotlib.pyplot as plt
import matplotlib as mpl  
mpl.rc('font',family='Times New Roman')

# Function to round to nearest base and add 1-unit of base in correct direction. Based on:
# https://stackoverflow.com/questions/2272149/round-to-5-or-other-number-in-python
def round_nearest(num, base=5):
    nearest = base*round(num/base)
    if num < 0: nearest -= base
    elif num > 0: nearest += base
    else: nearest = 0
    return nearest

# Function to find spatial count of atoms
def find(m, basename, spatial):
    # Find box info to compute ilo/ihi; i=x, y, and z
    xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
    xlo = float(xline[0]); xhi = float(xline[1]);
    ylo = float(yline[0]); yhi = float(yline[1]);
    zlo = float(zline[0]); zhi = float(zline[1]);
    
    # Generate bin criteria
    if spatial['direction'] == 'x': minbin = xlo; maxbin = xhi;
    elif spatial['direction'] == 'y': minbin = ylo; maxbin = yhi;
    elif spatial['direction'] == 'z': minbin = zlo; maxbin = zhi;
    else: raise Exception('ERROR direction not supported')
    
    # Set/find binning criteria
    binspan = spatial['bin-span'] # Each bin will span N angstroms (IE [0, N], [N, N+1], [N+1, N+2], ...)
    minbin = round_nearest(minbin, base=binspan); maxbin = round_nearest(maxbin, base=binspan); # round to nearest binspan

    # Loop through atoms and add count to data
    data = { i:{(n, n+binspan):0 for n in range(minbin, maxbin, binspan)} for i in m.formulaIDs_reverse } # { formula : bins }
    for i in m.atoms:
        formula = m.atoms[i].formulaID
        
        # Find postion based in desired direction
        if spatial['direction'] == 'x': pos = m.atoms[i].x
        elif spatial['direction'] == 'y': pos = m.atoms[i].y
        elif spatial['direction'] == 'z': pos = m.atoms[i].z
            
        # tally postion in data[formula] bins
        for lobin, hibin in data[formula]:
            if pos >= lobin and pos < hibin:
                data[formula][(lobin, hibin)] += 1
        
    ###############################
    # Plot Spatial count of atoms #
    ###############################
    def avg_bin(bindict):
        return [(bindim[1]+bindim[0])/2 for bindim in bindict]
    def count(bindict, div_factor):
        tmp = list(bindict.values())
        if div_factor is not None:
            tmp = [int(round(i/div_factor)) for i in tmp]
        return tmp
    
    # Set master colors
    if not spatial['colors']: colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    else: colors = spatial['colors']
        
    # Create new plot
    fig, ax = plt.subplots()

    # Plot all data found
    for n, i in enumerate(data):
        # find name to plot
        if i in spatial['formula2name']: name = spatial['formula2name'][i]
        else: name = i
        
        # find div_factor if applicable
        if i in spatial['normalize-atom-count']: div_factor = spatial['normalize-atom-count'][i]
        else: div_factor = None
        
        # Get x and y data to plot
        x = avg_bin(data[i]); y = count(data[i], div_factor);
        plt.plot(x, y, linewidth=2, label=name, color=colors[n%len(colors)])


    # Adjust axis and set label axis
    ax.set_xlabel('{}-{}'.format(spatial['direction'], 'Position ($\AA$)'), fontsize=12)
    ax.set_ylabel('{}'.format('Number of atoms per bin'), fontsize=12)
    plt.xticks(fontsize=12) 
    plt.yticks(fontsize=12) 

    # Set legend
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=5, fontsize=12)

    # Set figure size (width, height, forward to GUI)
    fig.set_size_inches(6, 4, forward=True)

    # otherwise the legend gets clipped off of saved image
    fig.tight_layout()

    # Save plot as
    figname = '{}{}'.format(basename, '_spatial_count_direction_' + spatial['direction'] + '.jpeg')
    fig.savefig(figname, dpi=300)
    
    
    #####################
    # Write data to csv #
    #####################
    # Generate headers_index and headers list
    n = 0
    headers_index = { '{}-Postion'.format(spatial['direction']): n }
    for i in data:
        n += 1
        headers_index[i] = n
    headers = len(headers_index)*['intialize'] # [ column1 header, column2 header, Ncolumns, ... ]
    for i in headers_index:
        if headers_index[i] == 0: headers[headers_index[i]] = '{}'.format(i)
        else: headers[headers_index[i]] = 'natoms {}'.format(i)
        
    # Generate rows to write to .csv file
    rows = [] # [ [row1], [row2],  nrows, ... ]; where rowN = [column1, column2, Ncolumns, ... ]
    for i in data:        
        # find div_factor if applicable
        if i in spatial['normalize-atom-count']: div_factor = spatial['normalize-atom-count'][i]
        else: div_factor = None
        
        # Get bin and count data to write to csv
        bin_lst = avg_bin(data[i]); count_lst = count(data[i], div_factor);
        
        # generate empty list in rows and append data in correct locations
        for n in range(len(bin_lst)):
            if len(rows) <= n: rows.append(len(headers)*[str(0)])
            
            # add postion into rows[n][0]
            rows[n][0] = str(bin_lst[n])
            
            # add formula data to rows[n][formula_index]
            rows[n][headers_index[i]] = str(count_lst[n])
    
    # Write info to .csv
    csvname = '{}{}'.format(basename, '_spatial_count_direction_' + spatial['direction'] + '.csv')
    with open(csvname, 'w') as f:
        # Write headers and rows
        f.write('{}\n'.format(', '.join(headers)))
        for row in rows:
            f.write('{}\n'.format(', '.join(row)))
    return