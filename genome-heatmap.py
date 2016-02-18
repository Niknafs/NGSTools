from __future__ import absolute_import, division
import argparse
import collections
import math
import os
import sys
from itertools import takewhile
import numpy as np

from Bio._py3k import map, range, zip
iteritems = (dict.iteritems if sys.version_info[0] < 3 else dict.items)

MB = 1e-6

# If running headless, use a suitable GUI-less plotting backend
if not os.environ.get('DISPLAY'):
    import matplotlib
    try:
        matplotlib.use("Agg", force=True)
    except TypeError:
        # Older matplotlib doesn't have 'force' argument
        matplotlib.use("Agg")

from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
pyplot.ioff()

#---------------------------------------------------------------#
class InputParser(object):
    def __init__(self):
        parser = argparse.ArgumentParser(usage="genome-wide heatmap module largely borrowing from and inspired by a similar module in CNVkit package (https://github.com/etal/cnvkit). Original function modified to run as a standalone script, , allow visualization of both continuous and categorical data, flexibility of input file format, and include a colorbar.")
        parser.add_argument('-p', '--path',
                    help="file with paths of the profiles to be visualized. In each file, the first three columns (coordinate), and the column indicated by index (signal) will be used [0-based].")

        parser.add_argument('-c', '--chromosome',
                            help="""Chromosome (e.g. 'chr1') or chromosomal range (e.g.
                            'chr1:2333000-2444000') to display. If a range is given,
                            all targeted genes in this range will be shown, unless
                            '--gene'/'-g' is already given.""")
        parser.add_argument('-l','--colormap',
                            help = 'label to color mapping for categorical signal')
        parser.add_argument('-r','--crange',
                            help = 'dynamic range of values in heatmap, in case of numerical signal')
        parser.add_argument('-i','--index',
                            help = 'index of column with signal to be plotted.')
        parser.add_argument('-o', '--output',
                            help="Output PDF file name.")
        args = parser.parse_args()
        self.__dict__.update(args.__dict__)
#---------------------------------------------------------------#        
def ifloat(string):
    try:
        return float(string)
    except:
        return -1000
#---------------------------------------------------------------#
def _cmd_heatmap(args):
    """Plot a signal for multiple samples as a heatmap."""
    filenames = [line.strip() for line in file(args.path)]
    create_heatmap(filenames, args.index, args.colormap, args.crange, args.chromosome)
    if args.output:
        pyplot.savefig(args.output, format='pdf', bbox_inches="tight")
        print "Output: " +  args.output
    else:
        pyplot.show()
#---------------------------------------------------------------#
def create_heatmap(filenames, index, colormap, crange, show_chromosome=None):
    """
    Plot signal for multiple samples as a heatmap.
    Signal can be continuous or categorical value.
    In case of numeric values, a continuous color
    scale will be used to map values to color. 
    In case of categorical value, the user should 
    provide mapping of values to color as an input"""
    
    _fig = pyplot.figure(figsize = (12,8))
    gs = matplotlib.gridspec.GridSpec(2,1,height_ratios = [30,1])
    axis = pyplot.subplot(gs[0])
    axis_aux = pyplot.subplot(gs[1])

    # List sample names on the y-axis
    axis.set_yticks([i + 0.5 for i in range(len(filenames))])
    axis.set_yticklabels(list(map(fbase, filenames)))
    axis.invert_yaxis()
    axis.set_ylabel("Samples")
    axis.set_axis_bgcolor('#DDDDDD')

    # Group each file's probes/segments by chromosome
    sample_data = [collections.defaultdict(list) for _f in filenames]
    
    #-------------------------------------------------------------#
    # read in the signal value in each sample from the input files
    index = int(index)
    for i, fname in enumerate(filenames):
        f_h = file(fname)
        for line in f_h:
            toks = line.strip().split('\t')
            if toks[0].lower() in ['chrom', 'chromosome']:
                # header line
                continue
            if colormap == None:
                # numerical data, convert signal to floating point numbers
                sample_data[i][toks[0]].append((int(toks[1]), int(toks[2]), ifloat(toks[index])))
            else:
                sample_data[i][toks[0]].append((int(toks[1]), int(toks[2]), toks[index]))

        f_h.close()
     #-------------------------------------------------------------#

    # Calculate the size (max endpoint value) of each chromosome
    chrom_sizes = {}
    for row in sample_data:
        for chrom, data in iteritems(row):
            max_posn = max(coord[1] for coord in data)
            chrom_sizes[chrom] = max(max_posn, chrom_sizes.get(chrom, 0))
    chrom_sizes = collections.OrderedDict(sorted(iteritems(chrom_sizes),
                                                 key=sorter_chrom_at(0)))

    if colormap != None:
        cvg2rgb = CVG2RGB(colormap = colormap)
    else:
        if crange != None:
            vmax = float(crange)
            vmin = (-1) * vmax
        else:
            vmax = get_max_abs_value(sample_data)
            vmin = (-1) * vmax

        # this matplotlib color map is appropriate for diverging data
        # ; e.g. copy number values. can be subsituted by any maplotlib
        # colormap
        my_cmap = matplotlib.cm.seismic
        color_norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
    
        cvg2rgb = CVG2RGB(rcmap = my_cmap, color_norm = color_norm)

        
    def plot_rect(y_idx, x_start, x_end, cvg):
        """Draw a rectangle in the given coordinates and color."""
        x_coords = (x_start, x_start, x_end + 1, x_end + 1)
        y_coords = (y_idx, y_idx + 1, y_idx + 1, y_idx)
        if cvg in [-1000.0, 'NA']:
            # missing data, a shade of dark gray
            rgbcolor =  (0.3, 0.3, 0.3)
        else:
            rgbcolor = cvg2rgb.get_color(cvg)
        axis.fill(x_coords, y_coords, color=rgbcolor)

    if show_chromosome:
        # Lay out only the selected chromosome
        chrom_offsets = {show_chromosome: 0.0}
        # Set x-axis the chromosomal positions (in Mb), title as the chromosome
        axis.set_xlim(0, chrom_sizes[show_chromosome] * MB)
        axis.set_title(show_chromosome)
        axis.set_xlabel("Position (Mb)")
        axis.tick_params(which='both', direction='out')
        axis.get_xaxis().tick_bottom()
        axis.get_yaxis().tick_left()
        # Plot the individual probe/segment coverages
        for i, row in enumerate(sample_data):
            for start, end, cvg in row[show_chromosome]:
                plot_rect(i, start * MB, end * MB, cvg)

    else:
        # Lay out chromosome dividers and x-axis labels
        # (Just enough padding to avoid overlap with the divider line)
        chrom_offsets = plot_x_dividers(axis, chrom_sizes, 1)
        # Plot the individual probe/segment coverages
        for i, row in enumerate(sample_data):
            for chrom, curr_offset in iteritems(chrom_offsets):
                for start, end, cvg in row[chrom]:
                    plot_rect(i, start + curr_offset, end + curr_offset, cvg)

    plot_y_dividers(axis, sample_data, color = (0.4, 0.4, 0.4))
    
    if colormap != None:
        cols = cvg2rgb.label2color.items()
        
        cols = [(index, cols[index][0], cols[index][1]) \
                for index in range(len(cols))]

        colors = [list(x[2]) for x in cols]
        labels = [x[1] for x in cols]
        cmap = matplotlib.colors.ListedColormap(colors)
        bounds = [col[0] + 0.5 for col in cols]
        bounds = [-0.5] + bounds
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        cb2 = matplotlib.colorbar.ColorbarBase(axis_aux, cmap = cmap,\
                                               boundaries = bounds, \
                                               norm = norm, \
                                               ticks = bounds, \
                                               orientation = 'horizontal')
        cb2.set_ticks([x + 0.5 for x in bounds])
        cb2.set_ticklabels(labels)

    else:
        cb1 = matplotlib.colorbar.ColorbarBase(axis_aux, cmap = my_cmap, norm = color_norm,\
                                               orientation = 'horizontal')
        pass
#---------------------------------------------------------------#
def get_max_abs_value(sample_data):
    val = 0
    for i in range(len(sample_data)):
        for chrom in sample_data[i].keys():
            vals = [data[2] for data in sample_data[i][chrom]]
            vals = filter(lambda x: x!= -1000, vals) 
            val = max(val, max(map(abs, vals)))
    return val
#---------------------------------------------------------------#
def read_class_color(path):
    # in case of categorical data
    # the use should provide a table with mappings of category labels
    # to rgb value
    # e.g.
    # positive 0.9 0.1 0.1
    # negative 0.1 0.1 0.9
    cols = collections.OrderedDict()
    for line in file(path):
        toks = line.strip().split('\t')
        cols[toks[0]] = tuple(map(float, toks[1:]))
    return cols
#---------------------------------------------------------------#
class CVG2RGB(object):
    def __init__(self, colormap = None, rcmap = None, color_norm = None ):
        self.colormap = colormap
        self.rcmap = rcmap
        self.color_norm = color_norm
        if self.colormap != None:
            self.label2color = read_class_color(self.colormap)        
    def get_color(self, cvg ):
        if self.colormap != None:
            return self.label2color.get(cvg)
        # if not, map the numerical value to a continuous color scale
        else:
            return self.rcmap(self.color_norm(cvg))
#-----------------------------------------------------------------#
def fbase(fname):
    """Strip directory and all extensions from a filename."""
    return os.path.basename(fname).split('.', 1)[0]
#-----------------------------------------------------------------#
def sorter_chrom_at(index):
    """Create a sort key function that gets chromosome label at a list index."""
    return lambda row: sorter_chrom(row[index])
#-----------------------------------------------------------------#
def sorter_chrom(label):
    """Create a sorting key from chromosome label.

    Sort by integers first, then letters or strings. The prefix "chr"
    (case-insensitive), if present, is stripped automatically for sorting.

    E.g. chr1 < chr2 < chr10 < chrX < chrY < chrM
    """
    # Strip "chr" prefix
    chrom = (label[3:] if label.lower().startswith('chr')
             else label)
    if chrom in ('X', 'Y'):
        key = (1000, chrom)
    else:
        # Separate numeric and special chromosomes
        nums = ''.join(takewhile(str.isdigit, chrom))
        chars = chrom[len(nums):]
        nums = int(nums) if nums else 0
        if not chars:
            key = (nums, '')
        elif len(chars) == 1:
            key = (2000 + nums, chars)
        else:
            key = (3000 + nums, chars)
    return key
#-----------------------------------------------------------------#
def plot_y_dividers(axis, sample_data, color = (0.45, 0.45, 0.45)):
    '''
    plot horizontal lines to mark the split between samples
    '''
    for i in range(len(sample_data)):
        axis.axhline(i, color = color)
#-----------------------------------------------------------------#
def plot_x_dividers(axis, chromosome_sizes, pad):
    """Plot vertical dividers and x-axis labels given the chromosome sizes.

    Returns a table of the x-position offsets of each chromosome.

    Draws vertical black lines between each chromosome, with padding.
    Labels each chromosome range with the chromosome name, centered in the
    region, under a tick.
    Sets the x-axis limits to the covered range.
    """
    assert isinstance(chromosome_sizes, collections.OrderedDict)

    x_dividers = []
    x_centers = []
    x_starts = collections.OrderedDict()
    curr_offset = pad
    for label, size in chromosome_sizes.items():
        x_starts[label] = curr_offset
        x_centers.append(curr_offset + 0.5 * size)
        x_dividers.append(curr_offset + size + pad)
        curr_offset += size + 2 * pad

    axis.set_xlim(0, curr_offset)
    for xposn in x_dividers[:-1]:
        axis.axvline(x=xposn, color='k')
    # Use chromosome names as x-axis labels (instead of base positions)
    axis.set_xticks(x_centers)
    axis.set_xticklabels(chromosome_sizes.keys(), rotation=60)
    axis.tick_params(labelsize='small')
    axis.tick_params(axis='x', length=0)
    axis.get_yaxis().tick_left()

    return x_starts
#-----------------------------------------------------------------#
def main():
    args = InputParser()
    _cmd_heatmap(args)

if __name__ == '__main__':
    main()

