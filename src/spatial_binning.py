# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
July 16th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import math


# Function to generate bins
def generate_bins(lo, hi, increment):
    bins = {} # {binID: (lo, hi, center), ...}
    length = hi - lo
    number = math.ceil(length/increment)
    delta = length/number
    for n in range(number):
        lower = (n)*delta + lo
        higher = (n+1)*delta + lo
        center = (higher + lower)/2
        bins[n+1] = (lower, higher, center)
    return bins, delta, number

# Function to generate data
def generate_data(bins, logger):
    data = {} # {binID: {log1:QTY, log2:QTY, ...}, ...}
    for i in bins:
        data[i] = {log:0 for log in logger}
    return data

# Function to find binID
def assign_to_bin(position, bins):
    binID = 0
    for binID in bins:
        lo, hi, center = bins[binID]
        if lo <= position <= hi: break
    return binID 

# Function to generate attribute coupling
def attribute_coupling(binning, atom, return_type='dict'):
    attributes = {i:set() for i in binning['logging']}
    for log in binning['logging']:
        if binning['delimiter'] in log:
            tmp = log.split(binning['delimiter'])
            tmp = [str(j).strip() for j in tmp]
            string = ''
            for n, j in enumerate(tmp):
                string += str(getattr(atom, j))
                if n+1 < len(tmp): string += str(binning['delimiter'])
            attributes[log].add(string)
        else: 
            string = str(getattr(atom, log))
            attributes[log].add(string)
    if return_type == 'list':
        lst = []
        for log in attributes:
            lst.extend(list(attributes[log]))
        return lst
    else:
        return attributes
    
# Function to update measure is it is a percent
def update_binned_measure(data, div, attributes, local=False):
    for i in data:
        if local: 
            divs = {i:0 for i in attributes} # { attribute or attribute pair : div factor }  
            attr = {} # { logged-attribute : logging-string}
            for k in data[i]:
                for attribute in attributes:
                    if k in attributes[attribute]: break
                divs[attribute] += data[i][k]
                attr[k] = attribute
        for j in data[i]:
            value = data[i][j]
            if local: div = divs[attr[j]]
            if div != 0:
                data[i][j] = 100*value/div
    return data

# Function to write files of logged information
def write_file(direction, rootname, binning, lo, hi, span, nbins, data, bins):
    filename = '{}_{}.txt'.format(rootname, direction)
    with open(filename,'w') as f:
        f.write('HEADER\n')
        f.write('ITEM:    {}lo    {}hi    bin-span    bin-count    measure    delimiter\n'.format(direction, direction))
        f.write('{:^10.6f} {:^10.6f} {:^10.6f} {:^10} {:^10} {:^10}\n'.format(lo, hi, span, nbins, binning['measure'], binning['delimiter']))
        
        f.write('\nDATA\n')
        sections = []; datas = [];
        for i in data:
            lower, upper, center = bins[i]
            values1 = [i, lower, upper, center]
            headers1 = ['binID', 'lo', 'hi', 'center']
            values2 = list(data[i].values())
            headers2 = list(data[i].keys())
            values = []
            headers = '    '.join(headers1 + headers2)
            float_ints = values1 + values2
            for n, i in enumerate(float_ints):
                tmp = str(i)
                if '.' in tmp:
                    tmp = '{:.4f}'.format(i)
                values.append(tmp)
            values = '    '.join(values)
            sections.append(headers)
            datas.append(values)
        f.write('ITEM: {}\n'.format(sections[0]))
        for i in datas:
            f.write('{}\n'.format(i))
    return


# Function to spatially bin quantities
def analysis(mm, binning, basename):    
    print('Finding spatially binned measures ....')

    # Find attributes to log
    if 'delimiter' not in binning:
        binning['delimiter'] = '|'
    attributes = {i:set() for i in binning['logging']}
    for i in mm.atoms:
        atom = mm.atoms[i]
        per_atom_attributes = attribute_coupling(binning, atom, return_type='dict')
        for log in per_atom_attributes:
            attributes[log].update(per_atom_attributes[log])
    
    # Create sorted logger
    logger = []
    for i in attributes:
        lst = list(attributes[i])
        lst = sorted(lst)
        logger.extend(lst)
    
    # Find set simulataion cell bounds
    xline = mm.xbox_line.split(); yline = mm.ybox_line.split(); zline = mm.zbox_line.split();
    xlo = float(xline[0]); xhi = float(xline[1])
    ylo = float(yline[0]); yhi = float(yline[1])
    zlo = float(zline[0]); zhi = float(zline[1])
    
    # Generate bins [ibins = {binID: (lo, hi, center), ...}]
    xbins, xspan, nx = generate_bins(xlo, xhi, binning['increment'])
    ybins, yspan, ny = generate_bins(ylo, yhi, binning['increment'])
    zbins, zspan, nz = generate_bins(zlo, zhi, binning['increment'])
    
    # Generate data [idata = {binID: {log1:QTY, log2:QTY, ...}, ...}]
    xdata = generate_data(xbins, logger)
    ydata = generate_data(ybins, logger)
    zdata = generate_data(zbins, logger)
        
    # tally atom attributes based on bins
    for i in mm.atoms:
        atom = mm.atoms[i]

        # Find measure to log
        measure = 0
        if binning['measure'] in ['size', 'psize', 'psize-local']:
            measure = 1
        elif binning['measure'] in ['mass', 'pmass', 'pmass-local']:
            measure = mm.masses[atom.type].coeffs[0]
        else:
            raise Exception(f"ERROR binning['measure'] {binning['measure']} is not supported. Supported measures 'size' or 'mass' or 'psize' or 'pmass' ")
        
        # Get per-atom attributes and update idata
        xbinID = assign_to_bin(atom.x, xbins)
        ybinID = assign_to_bin(atom.y, ybins)
        zbinID = assign_to_bin(atom.z, zbins)
        per_atom_attributes = attribute_coupling(binning, atom, return_type='list')
        for attribute in per_atom_attributes:
            xdata[xbinID][attribute] += measure
            ydata[ybinID][attribute] += measure
            zdata[zbinID][attribute] += measure
            
    # if measure is 'psize' or 'pmass' update quantities
    if binning['measure'] in ['psize', 'pmass', 'psize-local', 'pmass-local']:
        if binning['measure'] in ['psize-local', 'pmass-local']:
            local = True
            div = 1
        else: local = False
        if binning['measure'] == 'psize': div = mm.total_system_size
        if binning['measure'] == 'pmass': div = mm.total_system_mass
        xdata = update_binned_measure(xdata, div, attributes, local)
        ydata = update_binned_measure(ydata, div, attributes, local)
        zdata = update_binned_measure(zdata, div, attributes, local)
        
    # Write files
    rootname = '{}_spatial_binning'.format(basename)
    if 'x' in binning['direction']: write_file('x', rootname, binning, xlo, xhi, xspan, nx, xdata, xbins)
    if 'y' in binning['direction']: write_file('y', rootname, binning, ylo, yhi, yspan, ny, ydata, ybins)
    if 'z' in binning['direction']: write_file('z', rootname, binning, zlo, zhi, zspan, nz, zdata, zbins)
        
        
    return mm