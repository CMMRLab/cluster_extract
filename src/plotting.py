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
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np
import math

import matplotlib as mpl  
mpl.rc('font',family='Times New Roman')


#################################################################################
# Function for plotting multiple bars per one location via stack over flow:     #
# https://stackoverflow.com/questions/14270391/python-matplotlib-multiple-bars  #
# Some slight modifications were preform by Josh Kemppainen to make plot what   #
# was desired for this project.                                                 #  
#################################################################################
def bar_plot(ax, data, colors=None, total_width=0.8, single_width=0.7, legend=True):
    # Check if colors where provided, otherwhise use the default color cycle
    if colors is None: colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Number of bars per group
    n_bars = len(data)

    # The width of a single bar
    bar_width = total_width/n_bars; bars = []
    for i, (name, values) in enumerate(data.items()):
        # The offset in x direction of that bar
        x_offset = (i - n_bars / 2) * bar_width + bar_width / 2

        # Draw a bar for every value of that type
        for x, y in enumerate(values):
            bar = ax.bar(x+1+x_offset, y, width=bar_width*single_width, color=colors[i % len(colors)])
    
        # Add a handle to the last drawn bar, which we'll need for the legend
        bars.append(bar[0])

    # Draw legend if we need
    if legend:
        ax.legend(bars, data.keys(), loc='upper center', bbox_to_anchor=(0.5, 1.2), fancybox=True, ncol=5, fontsize=12)
    return

# Moving average function
def moving_average(x, w=5):
    return list(np.convolve(x, np.ones(w), 'valid') / w)


# Function to round to nearest base and add 1-unit of base in correct direction. Based on:
# https://stackoverflow.com/questions/2272149/round-to-5-or-other-number-in-python
def round_nearest(num, base=5):
    nearest = base*round(num/base)
    if num < 0: nearest -= base
    elif num > 0: nearest += base
    else: nearest = 0
    return nearest



##########################################
# Mast function to plot all desired data #
##########################################
def master(m, basename, plot_options):
    ########################################################
    # Find Distance and orientation for different formulas #
    ########################################################
    unique_formulas = [i for i in m.formulaIDs_reverse]
    dist_formula = {i:[] for i in unique_formulas} # { formula : [lst of distances ] }
    angle_formula = {i:[] for i in unique_formulas} # { formula : [lst of angles ] }
    ordered_formulas = [] # Will keep formulas in order based on distance 
    for i in m.dist_orientation_molecules:
        dist, plane_angle, angles = m.dist_orientation_molecules[i]
        #formula1 = m.molecules.data[i[0]].formula # ref molecule
        formula2 = m.datatable[i[1]].formula # test molecule
        
        dist_formula[formula2].append(dist)
        angle_formula[formula2].append(plane_angle)
        if formula2 not in ordered_formulas:
            ordered_formulas.append(formula2)
        
        
        
    #####################
    # Set master colors #
    #####################
    if not plot_options['colors']:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    else:
        colors = plot_options['colors']
        
        
        
    ################################
    # Plot Distance vs Orientation #
    ################################
    # Create new plot
    fig, ax = plt.subplots()

    # Plot all data found
    n = 0; dists = []
    for i in ordered_formulas:
        if len(dist_formula[i]) > plot_options['min-molecule-count'] and len(dist_formula[i]) > plot_options['min-molecule-count']:
            if i in plot_options['formula2name']: name = plot_options['formula2name'][i]
            else: name = i
            plt.plot(dist_formula[i], angle_formula[i], '.', mfc='white', markersize=20, linewidth=10, label=name, color=colors[n%len(colors)])
            dists.append(max(dist_formula[i])) # For xlim
            dists.append(min(dist_formula[i])) # For xlim
            n += 1

    # Adjust axis and set label axis
    ax.set_xlabel('{}'.format('Perpendicular Distance ($\AA$)'), fontsize=12)
    ax.set_ylabel('{}'.format('Orientation of Molecule ($^\circ$)'), fontsize=12)
    plt.xticks(fontsize=12) 
    plt.yticks(fontsize=12) 

    # Set axis limits and ticks
    xlo = round_nearest( math.floor(min(dists)), base=plot_options['dist-bin-span'] )
    xhi = round_nearest( math.ceil(max(dists)), base=plot_options['dist-bin-span'] )
    plt.xlim([xlo, xhi])

    # Set legend
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=5, fontsize=12)

    # Set figure size (width, height, forward to GUI)
    fig.set_size_inches(6, 4, forward=True)

    # otherwise the legend gets clipped off of saved image
    fig.tight_layout()

    # Save plot as
    figname = '{}{}'.format(basename, '_dist_vs_orientation.jpeg')
    fig.savefig(figname, dpi=300)
    
    
    ##################################
    # Plot Distance binned Histogram #
    ##################################
    # Create new plot
    fig, ax = plt.subplots()
    
    # Set/find binning criteria
    distances = [j for i in ordered_formulas for j in dist_formula[i]]
    binspan = plot_options['dist-bin-span'] # Each bin will span N angstroms (IE [0, N], [N, N+1], [N+1, N+2], ...)
    minbin = math.floor(min(distances)); maxbin = math.ceil(max(distances)); # Set minbin and find maxbin
    minbin = round_nearest(minbin, base=binspan); maxbin = round_nearest(maxbin, base=binspan); # round to nearest binspan
    
    # For angle plotting domain
    dist_abs = [abs(i) for i in distances] 
    mindist = math.floor(min(dist_abs))
    maxdist = math.ceil(max(dist_abs))
    
    # Build dictionary to hold bar plot data
    multi_bar_data = {} # { formula : [lst of count in each bin] }
    for i in ordered_formulas:
        if len(dist_formula[i]) > plot_options['min-molecule-count']:
            bins = {(n, n+binspan):0 for n in range(minbin, maxbin, binspan)} # intialize bins count as zeros
            for dist in dist_formula[i]:
                for lobin, hibin in bins:
                    if dist >= lobin and dist < hibin:
                        bins[(lobin, hibin)] += 1
            #multi_bar_data[i] = sorted( list(bins.values()), key=lambda x: abs(x-0) )
            multi_bar_data[i] = list(bins.values())
            
    # If zeros switch to 0.1 to have bar show up in plot (comment/uncomment is not-desired/desired)
    multi_bar_data = {i:[j if j != 0 else 0.1 for j in multi_bar_data[i]] for i in multi_bar_data}
            
    # Adjust key names if formula in 'formula2name'
    if plot_options['formula2name']:
        for k, v in list(multi_bar_data.items()):
            multi_bar_data[plot_options['formula2name'].get(k, k)] = multi_bar_data.pop(k)
        
    bar_plot(ax, multi_bar_data, colors=colors, legend=True)
    plt.rc('xtick', labelsize=12)    # fontsize of the xticks
    plt.rc('ytick', labelsize=12)    # fontsize of the yticks
    plt.rc('axes', labelsize=12)     # fontsize of the axes title
    plt.rc('figure', titlesize=12)   # fontsize of the figure title
    plt.xlabel('{}'.format('Binned Distances ($\AA$)'), fontsize=12)
    plt.ylabel('{}'.format('Count of Atoms in Molecule'), fontsize=12)
    
    bins = [(n, n+binspan) for n in range(minbin, maxbin, binspan)]
    string_bins = ['{} to {}'.format(i[0], i[1]) for i in bins]
    string_bins.insert(0, '')
    plt.xticks(range(len(string_bins)), string_bins, fontsize=12, rotation=45, ha='right') 
    plt.yticks(fontsize=12) 
    ax.yaxis.set_major_locator(MaxNLocator(integer=True)) # Force all y-axis values to be ints
    
    # Set figure size (width, height, forward to GUI)
    fig.set_size_inches(6, 4, forward=True)

    # otherwise the legend gets clipped off of saved image
    fig.tight_layout()

    figname = '{}{}'.format(basename, '_binned_dist_bar_chart.jpeg')
    fig.savefig(figname, dpi=300)
    
    
    ########################################################################
    # Generate and Plot distance psuedo-RDF termed as a DDF via bin center #
    ########################################################################
    distances = [abs(j) for i in ordered_formulas for j in dist_formula[i]]
    binspan = plot_options['ddf-bin-span'] # Each bin will span N angstroms (IE [0, N], [N, N+1], [N+1, N+2], ...)
    minbin = 0; maxbin = math.ceil(max(distances)); # Set minbin and find maxbin
    nbins = round((maxbin-minbin)/binspan)
    
    # Generate ddf data
    ddf = {} # d=distance;d=distribution;f=function { formula name : bins }  where bins = { (0,2):count, (2,4):count, ... }
    for i in ordered_formulas:
        if len(dist_formula[i]) > plot_options['min-molecule-count']:
            bins = {( n*binspan, (n+1)*binspan ):0 for n in range(nbins)} # intialize bins count as zeros
            for dist in dist_formula[i]:
                dist = abs(dist)
                for lobin, hibin in bins:
                    if dist >= lobin and dist < hibin:
                        bins[(lobin, hibin)] += 1
            ddf[i] = bins
            
    # Create new plot
    fig, ax = plt.subplots()
    
    bins = [( n*binspan, (n+1)*binspan ) for n in range(nbins)] # generate bins for plotting
    avgs = [(i[0]+i[1])/2 for i in bins] # center of bins
    if plot_options['ddf-max'] != 'all':
        maxdist = float(plot_options['ddf-max'])
    else: maxdist = max(distances)
    for n, formula in enumerate(ddf):
        count = list(ddf[formula].values())
        if formula in plot_options['formula2name']: name = plot_options['formula2name'][formula]
        else: name = formula
        
        # Apply moving average
        x = avgs
        y = moving_average(count, w=plot_options['ddf-mv-avg'])
        
        # reduce data down by max length
        x1 = []; y1 = [];
        for xx, yy in zip(x, y):
            if xx <= maxdist:
                x1.append(xx); y1.append(yy)
                
        # Plot anaylzed data
        plt.plot(x1, y1, ls='-', linewidth=3, label=name, color=colors[n%len(colors)])
        
    # Adjust axis and set label axis
    ax.set_xlabel('{}'.format('Perpendicular Distance ($\AA$)'), fontsize=12)
    ax.set_ylabel('{}'.format('Molecule Count'), fontsize=12)
    plt.xticks(fontsize=12) 
    plt.yticks(fontsize=12) 
    #plt.xlim([0, maxdist])
    ax.yaxis.set_major_locator(MaxNLocator(integer=True)) # Force all y-axis values to be ints
    
    # Set legend
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=5, fontsize=12)

    # Set figure size (width, height, forward to GUI)
    fig.set_size_inches(6, 4, forward=True)

    # otherwise the legend gets clipped off of saved image
    fig.tight_layout()

    # Save plot as
    figname = '{}{}'.format(basename, '_ddf_bin_center.jpeg')
    fig.savefig(figname, dpi=300)
    
    
    #################################################################
    # Generate and Plot distance psuedo-RDF termed as a DDF via RMS #
    #################################################################
    distances = [abs(j) for i in ordered_formulas for j in dist_formula[i]]
    binspan = plot_options['ddf-bin-span'] # Each bin will span N angstroms (IE [0, N], [N, N+1], [N+1, N+2], ...)
    minbin = math.floor(min(distances)); maxbin = math.ceil(max(distances)); # Set minbin and find maxbin
    nbins = round((maxbin-minbin)/binspan)
    
    # Generate ddf data
    ddf = {} # d=distance;d=distribution;f=function { formula name : bins }  where bins = { (0,2):count, (2,4):count, ... }
    for i in ordered_formulas:
        if len(dist_formula[i]) > plot_options['min-molecule-count']:
            bins = {( n*binspan, (n+1)*binspan ):[] for n in range(nbins)} # intialize bins count as zeros
            for dist in dist_formula[i]:
                dist = abs(dist)
                for lobin, hibin in bins:
                    if dist >= lobin and dist < hibin:
                        bins[(lobin, hibin)].append(dist)
            ddf[i] = bins
    
    # Functions need for plotting
    def compute_RMS(lst, binspan):
        if lst: rms = math.sqrt(sum([i*i for i in lst])/len(lst))
        else: rms = (binspan[0]+binspan[1])/2
        return rms

    # Create new plot
    fig, ax = plt.subplots()
    
    bins = [( n*binspan, (n+1)*binspan ) for n in range(nbins)] # generate bins for plotting
    if plot_options['ddf-max'] != 'all':
        maxdist = float(plot_options['ddf-max'])
    else: maxdist = max(distances)
    for n, formula in enumerate(ddf):
        tmp = list(ddf[formula].values())
        count = [len(i) for i in tmp]
        rms = [compute_RMS(i, bins[n]) for n, i in enumerate(tmp)]
        if formula in plot_options['formula2name']: name = plot_options['formula2name'][formula]
        else: name = formula
        
        # Apply moving average to y-data and make the same length as z-data
        count = moving_average(count, w=plot_options['ddf-mv-avg'])
        x = []; y = [];
        for xx, yy in zip(rms, count):
            if xx <= maxdist and yy > 0: # make sure xx small enough and yy is large enough
                x.append(xx); y.append(yy);
        
        # Plot anaylzed data
        plt.plot(x, y, ls='-', linewidth=3, label=name, color=colors[n%len(colors)])

        
    # Adjust axis and set label axis
    ax.set_xlabel('{}'.format('RMS Perpendicular Distance ($\AA$)'), fontsize=12)
    ax.set_ylabel('{}'.format('Molecule Count'), fontsize=12)
    plt.xticks(fontsize=12) 
    plt.yticks(fontsize=12) 
    #plt.xlim([0, maxdist])
    ax.yaxis.set_major_locator(MaxNLocator(integer=True)) # Force all y-axis values to be ints
    
    # Set legend
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=5, fontsize=12)

    # Set figure size (width, height, forward to GUI)
    fig.set_size_inches(6, 4, forward=True)

    # otherwise the legend gets clipped off of saved image
    fig.tight_layout()

    # Save plot as
    figname = '{}{}'.format(basename, '_ddf_RMS.jpeg')
    fig.savefig(figname, dpi=300)
    
    
    #####################
    # Write data to csv #
    #####################
    # Generate headers_index and headers list
    n = 0
    headers_index = { '{}'.format('RMS Perpendicular Distance'): n }
    for i in ddf:
        n += 1
        headers_index[i] = n
    headers = len(headers_index)*['intialize'] # [ column1 header, column2 header, Ncolumns, ... ]
    for i in headers_index:
        if headers_index[i] == 0: headers[headers_index[i]] = '{}'.format(i)
        else: headers[headers_index[i]] = 'count {}'.format(i)
    
    # Generate rows to write to .csv file
    rows = [] # [ [row1], [row2],  nrows, ... ]; where rowN = [column1, column2, Ncolumns, ... ]
    for n, formula in enumerate(ddf):
        tmp = list(ddf[formula].values())
        count = [len(i) for i in tmp]
        rms = [compute_RMS(i, bins[n]) for n, i in enumerate(tmp)]
        
        # Apply moving average
        x = rms
        y = moving_average(count, w=plot_options['ddf-mv-avg'])
        
        # reduce data down by max length
        x1 = []; y1 = [];
        for xx, yy in zip(x, y):
            if xx <= maxdist:
                x1.append(xx); y1.append(yy)
                
        # generate empty list in rows and append data in correct locations
        for n in range(len(x1)):
            if len(rows) <= n: rows.append(len(headers)*[str(0)])
            
            # add postion into rows[n][0]
            rows[n][0] = str(x1[n])
            
            # add formula data to rows[n][formula_index]
            rows[n][headers_index[formula]] = str(y1[n])
    
    # Write info to .csv
    csvname = '{}{}'.format(basename, '_ddf_RMS.csv')
    with open(csvname, 'w') as f:
        # Write headers and rows
        f.write('{}\n'.format(', '.join(headers)))
        for row in rows:
            f.write('{}\n'.format(', '.join(row)))
    
    
    
    ###############################
    # Plot Angle binned Histogram #
    ###############################
    # Create new plot
    fig, ax = plt.subplots()
    
    # Set/find binning criteria
    angles = [j for i in ordered_formulas for j in angle_formula[i]]
    binspan = plot_options['angle-bin-span'] # Each bin will span N angstroms (IE [0, N], [N, N+1], [N+1, N+2], ...)
    minbin = math.floor(min(angles)); maxbin = math.ceil(max(angles)); # Set minbin and find maxbin
    minbin = round_nearest(minbin, base=binspan); maxbin = round_nearest(maxbin, base=binspan); # round to nearest binspan
    nbins = round((maxbin-minbin)/binspan)
    
    # Setting distance range based on user inputs
    if plot_options['dist-range4angles']:
        dist_range = plot_options['dist-range4angles']
    else: dist_range = [mindist, maxdist] # else user entire distance spectrum
    
    # Build dictionary to hold bar plot data
    multi_bar_data = {} # { formula : [lst of count in each bin] }
    for i in ordered_formulas:
        if len(angle_formula[i]) > plot_options['min-molecule-count']:
            bins = {( n*binspan, (n+1)*binspan ):0 for n in range(nbins) if n+1 <= plot_options['max-angle-bins']} # intialize bins count as zeros
            for n, angle in enumerate(angle_formula[i]):
                dist = dist_formula[i][n]
                if abs(dist) >= min(dist_range) and abs(dist) < max(dist_range):
                    for lobin, hibin in bins:
                        if angle >= lobin and angle < hibin:
                            bins[(lobin, hibin)] += 1
            multi_bar_data[i] = list(bins.values())
    
    # If zeros switch to 0.1 to have bar show up in plot (comment/uncomment is not-desired/desired)
    multi_bar_data = {i:[j if j != 0 else 0.1 for j in multi_bar_data[i]] for i in multi_bar_data}
            
    # Adjust key names if formula in 'formula2name'
    if plot_options['formula2name']:
        for k, v in list(multi_bar_data.items()):
            multi_bar_data[plot_options['formula2name'].get(k, k)] = multi_bar_data.pop(k)
            
    # Write data to log file
    with open(basename+'.txt', 'a') as f:
        bins = [( n*binspan, (n+1)*binspan ) for n in range(nbins) if n+1 <= plot_options['max-angle-bins']]
        string_bins = ['{} to {}'.format(i[0], i[1]) for i in bins]
        stringspan = len(string_bins[0])
        string_bins = '        '.join(string_bins) # 8 spaces
        name = 'Orientation binned data'
        column = '{:<20}'.format('formula')
        dashes = ''.join(int((len(string_bins)+len(column)-len(name))/2)*['-'])
        
        # Write header
        f.write('\n\n\n-{}{}{}-\n'.format(dashes, name, dashes))
        f.write('{} {}\n'.format(column, string_bins))
        for i in multi_bar_data:
            tmp = '{:<20}'.format(i)
            for n, j in enumerate(multi_bar_data[i]):
                num = str(j)
                size = int(stringspan+8-round(len(num)/2))
                space = '{:>{s}}'.format(' ', s=size/2)
                tmp += '{:>}{:>}{:>}'.format(space, num, space)
            f.write('{}\n'.format(tmp))
                
            
            
    # Plot based on colors or not
    if not plot_options['colors']:  bar_plot(ax, multi_bar_data, legend=True)
    else: bar_plot(ax, multi_bar_data, colors=colors, legend=True)
    plt.rc('xtick', labelsize=12)    # fontsize of the xticks
    plt.rc('ytick', labelsize=12)    # fontsize of the yticks
    plt.rc('axes', labelsize=12)      # fontsize of the axes title
    plt.rc('figure', titlesize=12)   # fontsize of the figure title
    plt.xlabel('{} {} +- {} to {} $\AA$'.format('Binned Orientations ($^\circ$)', 'within distance range: ', min(dist_range), max(dist_range)), fontsize=12)
    plt.ylabel('{}'.format('Count of Molecules'), fontsize=12)
    
    bins = [(n, n+binspan) for c, n in enumerate(range(minbin, maxbin, binspan)) if c+1 <= plot_options['max-angle-bins']]
    string_bins = ['{} to {}'.format(i[0], i[1]) for i in bins]
    string_bins.insert(0, '')
    plt.xticks(range(len(string_bins)), string_bins, fontsize=12, rotation=45, ha='right') 
    plt.yticks(fontsize=12) 
    ax.yaxis.set_major_locator(MaxNLocator(integer=True)) # Force all y-axis values to be ints
    
    # Set figure size (width, height, forward to GUI)
    fig.set_size_inches(6, 4, forward=True)

    # otherwise the legend gets clipped off of saved image
    fig.tight_layout()

    figname = '{}{}'.format(basename, '_binned_orientation_bar_chart.jpeg')
    fig.savefig(figname, dpi=300)
    return