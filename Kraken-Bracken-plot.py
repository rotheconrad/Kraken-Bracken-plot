#!/usr/bin/env python

'''
Builds bar plots from Bracken results with optional adjustment for
unclassified reads from Kraken results.

Bracken modifies Kraken results for a more acurate abundance estimate
which is great but Bracken does not include the unclassified reads. I
think it is important for a thorough analysis to report the fraction of
the metagenome that is not classified compared to what is classified.
Environmental samples are frequently >50% unclassified.

Therefore, this script has the option of accepting the Kraken reports
that the Bracken reports were generated from. It will read the total
number of unclassified reads from the Kraken report and adjust the
"Fraction of Total Reads" reported by Bracken because this "Fraction of
Total Reads" is only the fractional total of reads assigned to the
"Level ID" specified (the unclassified reads are not counted). With this
optional adjustment, our plots will show the fraction of total reads
(including unclassified).

input file: List of Bracken report files, SampleName

-i, input_file_specificiations:
    This is a .txt file with one entry per line in the order you want.
    The output plot will follow the order in this file.

    Example:
                Path/to/BrackenFilename\tSampleA
                Path/to/BrackenFilename\tSampleB
                Path/to/BrackenFilenmae\tSampleC

optional input:

-k, kraken_file_specifications:
    This is a .txt file with one entry per line. The Sample names must
    match the input file sample neames and the Kraken files should 
    correspond to those use for the Bracken files. The order of this
    file does not matter but I find it is easiest if the order is the
    same as the input file.

    Example:
                Path/to/KrakenFilename\tSampleA
                Path/to/KrakenFilename\tSampleB
                Path/to/KrakenFilename\tSampleC

-l, custom_legend_file:
    This is a .txt file with one entry per line of classification, color.
    The classification must match those used in the input file.
    The colors can be whatever you want. I prefer to use hex colors.
    Example: 
                Alteromonadales\t#2ca25f
                Sphingobacteriales\t#8856a7
                Rickettsiales\t#43a2ca

output: stacked bar plots in vector pdf format. Publication ready.

Ref Links:
https://ccb.jhu.edu/software/kraken2
https://github.com/DerrickWood/kraken2
https://ccb.jhu.edu/software/bracken
https://github.com/jenniferlu717/Bracken

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: December 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from os import listdir
from os.path import isfile, join
from collections import defaultdict
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import pandas as pd


# Set the default sans-serif font to Helvetica
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# Set default font to sans-serif
matplotlib.rcParams['font.family'] = "sans-serif"

class TickRedrawer(matplotlib.artist.Artist):
    # this is just to get the stupid ticks to draw right on the plot
    #https://stackoverflow.com/questions/19677963/
    #matplotlib-keep-grid-lines-behind-the-graph-but-the-y-and-x-axis-above
    """Artist to redraw ticks."""

    __name__ = "ticks"

    zorder = 10

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer: matplotlib.backend_bases.RendererBase) -> None:
        """Draw the ticks."""
        if not self.get_visible():
            self.stale = False
            return

        renderer.open_group(self.__name__, gid=self.get_gid())

        for axis in (self.axes.xaxis, self.axes.yaxis):
            loc_min, loc_max = axis.get_view_interval()

            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if tick.get_visible() and loc_min <= tick.get_loc() <= loc_max:
                    for artist in (tick.tick1line, tick.tick2line):
                        artist.draw(renderer)

        renderer.close_group(self.__name__)
        self.stale = False


def parse_infile(infile):
    # Parses the infiles (Bracken report files)
    # returns dict of {filename: samplename}
    # returns list of samplename order

    file_lookup = {}
    smpl_order = []

    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            filename = X[0]
            samplename = X[1]
            file_lookup[filename] = samplename
            smpl_order.append(samplename)

    return file_lookup, smpl_order


def parse_krakenfile(krakenfile):
    # Parses the infiles (Kraken report files) if provided
    # returns dict of {samplename: unclassified fraction as float}
    kraken_lookup = {}

    with open(krakenfile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            filename = X[0]
            samplename = X[1]

            with open(filename, 'r') as f:
                for l in f:
                    X = l.rstrip().split('\t')
                    fraction = float(X[0]) / 100 # convert percent to fraction
                    taxID = X[-1]
                    if taxID == "unclassified":
                        kraken_lookup[samplename] = fraction
                        break
                        print(l)

    return kraken_lookup


def parse_bracken(infile, krakenfile):
    """ Read the inputs and build a dataframe with index of Sample Names
    and columns for each TaxID """

    # parse infile
    file_lookup, smpl_order = parse_infile(infile)
    # if krakenfile, parse krakenfile
    if krakenfile:
        kraken_lookup = parse_krakenfile(krakenfile)
    # set as none. change later if krakenfile true.
    # will become column with unclassified fraction value
    unclassified = None

    # Initialize variables
    # data stores by {Sample: {class: value}}
    data = {}
    # List of global classifications from all samples
    Classifications = {}

    # Grab list of blast files from blast directory
    file_list = [f for f in file_lookup.keys()]

    # read through all files in input directory
    for file in file_list:
        print(f'\t\tReading File: {file}')
        # Set sample name
        sample = file_lookup[file]
        if krakenfile:
            krakenmod = 1 - kraken_lookup[sample]
        # Store sample data
        smpldata = {}
        # parse file
        with open(file, 'r') as f:
            # retrieve header to skip it
            header = f.readline()
            # read each line of file
            for l in f:
                # Grab the columns of interest
                X = l.rstrip().split('\t')
                classification = X[0]
                fraction = float(X[-1])
                if krakenfile:
                    fraction = fraction * krakenmod
                    #print(X[-1], krakenmod, fraction)
                # Store the data
                smpldata[classification] = fraction
                Classifications[classification] = ''

        data[sample] = smpldata

    # reorganize raw data to prep for data frame
    data2 = defaultdict(list)

    for smpl, classes in data.items():
        data2['Sample'].append(smpl)
        for clss in Classifications.keys():
            data2[clss].append(classes.setdefault(clss, 0))

    # convert to dataframe
    df = pd.DataFrame(data2)
    df = df.set_index('Sample')
    # make sure dataframe is in the specified order
    df = df.reindex(smpl_order)
    # grab unclassified fraction if krakenfile is true
    if krakenfile:
        unclassified = [kraken_lookup[sample] for sample in smpl_order]

    return df, unclassified


def parse_legend(lgnd_file):
    # parses the optional legend file if provided
    # returns lists of legend labels and matching colors.
    lgnd_lbs = []
    lgnd_clr = []

    with open(lgnd_file, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            label = X[0]
            color = X[1]
            lgnd_lbs.append(label)
            lgnd_clr.append(color)

    return lgnd_lbs, lgnd_clr


def plot_prep(df, unclassified, lgnd_file, top_classes):
    # prep the data for easy plotting.
    # Selects classes specified in legend file or top_classes number
    # creates other column with sums of classes not selected
    # creats unclassified column is krakenfile is provided
    # returns data frame to plot plus list of colors.

    # Turn off SettingWithCopyWarning from df2 & df3 copies
    pd.options.mode.chained_assignment = None 

    if lgnd_file:
        # optional
        lgnd_lbs, lgnd_clr = parse_legend(lgndfile)
    else:
        # default behavior
        totals = df.sum().sort_values(ascending=False)
        # Select top_classes number of global most abundant classifications.
        lgnd_lbs = totals.index.tolist()[:top_classes]
        # Generate a matching number of colors to use in the plot.

        ''' Using cmaps. First option I tried but I went a different direction.
        # First we select these 10 distinct colors as a cmap
        lgnd_cmap = plt.cm.get_cmap('tab10', 10)
        # But I want the cmap as a list. convert cmap to list
        # initialize list
        lgnd_clrs = []
        # loop through length of the cmap and add each color to the list
        for c in range(lgnd_cmap.N):
            color = lgnd_cmap(c)
            # oh and conver the colors to hex colors. I like those.
            lgnd_clrs.append(matplotlib.colors.rgb2hex(color))
        '''

        # List of 10 distinct colors generated from:
        # https://mokole.com/palette.html
        lgnd_clrs = [
                        '#2f4f4f','#2e8b57','#800000','#000080','#9acd32',
                        '#ff0000','#ff8c00','#ffd700','#00ff00','#ba55d3',
                        '#00fa9a','#e9967a','#00ffff','#0000ff','#ff00ff',
                        '#1e90ff','#eee8aa','#dda0dd','#ff1493','#7b68ee',
                        ]
        '''
        # List of 10 distinct colors generated from:
        # https://sashamaps.net/docs/resources/20-colors/
        lgnd_clrs = [
                        '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                        '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
                        '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000',
                        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9',
                        ]
        '''

        # Now we repeat these colors to the desired length
        # set counter equal to number of top_classes specified
        counter = top_classes
        # initialize list
        lgnd_clr = []
        # iterate through top_classes appending colors until length is reached.
        while counter >= 0:
            if counter > 20:
                lgnd_clr.extend(lgnd_clrs)
            else:
                lgnd_clr.extend(lgnd_clrs[:counter])
            counter -= 20

    # select lgnd_labels rows from df
    df2 = df[lgnd_lbs]
    # create "other" category by summing remaining classifactions
    df3 = df[df.columns.difference(lgnd_lbs)]
    df2['Other'] = df3.sum(axis=1)
    # add gray color for other
    lgnd_clr.append('#d3d3d3')
    # add unclassified if krakenfile is True
    if unclassified:
        df2['Unclassified'] = unclassified
        # add black color for unclassified
        lgnd_clr.append('#000000')

    print('\n\n', df2)

    return df2, lgnd_clr


def bracken_plot(df, outfile, lgnd_clmns, lgnd_clr):
    """ Plot Stacked Bar Charts by Sample """

    print('\n\nBuilding the plot...')

    # plot the data.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5,10), dpi=300)

    ax1 = df.plot.barh(
                    stacked=True,
                    ax=ax1,
                    width=.98,
                    color=lgnd_clr,
                    legend=False,
                    #alpha=0.6
                    )

    # set the axis parameters / style
    ax1.set_ylabel('Sample', fontsize=12)
    ax1.set_xlabel('Fractional Read Abundance', fontsize=12)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    if 'Unclassified' in df.columns:
        x = round(1-df['Unclassified'].min(), 2) + 0.01
        ax1.set_xlim([0,x])
    else:
        ax1.set_xlim([0,1])
    ax1.minorticks_on()
    ax1.tick_params(axis='y', which='minor', left=False)
    ax1.tick_params(labelsize=12, direction='inout', width=2, length=6)
    # set grid style
    ax1.xaxis.grid(
        which="minor", color='#737373', linestyle='--', linewidth=1
        )
    ax1.xaxis.grid(
        which="major", color='#737373', linestyle='--', linewidth=1
        )
    #ax1.set_axisbelow(True)
    ax1.add_artist(TickRedrawer())
    ax1.invert_yaxis()
    # Plot legend
    ax2.set_axis_off()
    handles, labels = ax1.get_legend_handles_labels()
    ax2.legend(
        handles,
        labels,
        ncol=lgnd_clmns,
        loc='center',
        fontsize=12,
        markerscale=2,
        frameon=False
        )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.subplots_adjust(wspace=0.05)
    plt.savefig(outfile)
    plt.close()

    return True


def main():
    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file_specifications',
        help='Text file with path/to/filename.bracken\tSampleName.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-k', '--kraken_file_specifications',
        help='(Optional) Text file with path/to/filename.kraken\tSampleName.',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-l', '--custom_legend_file',
        help='(Optional) List of classification, color for the legend.',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-n', '--number_legend_columns',
        help='(Optional) Number of columns for the legend (Default=2).',
        metavar='',
        type=int,
        default=2,
        required=False
        )
    parser.add_argument(
        '-t', '--number_top_classes',
        help='(Optional) Number most abundant classes to plot (Default=20).',
        metavar='',
        type=int,
        default=20,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n\n')

    # define input params
    infile = args['input_file_specifications']
    outfile = args['output_file']
    krakenfile = args['kraken_file_specifications']
    lgnd_file = args['custom_legend_file']
    lgnd_clmns = args['number_legend_columns']
    top_classes = args['number_top_classes']

    # Read in the Kraken files and build a dataframe with columns:
    # Sample, Classification, Fraction
    df, unclassified = parse_bracken(infile, krakenfile)

    # setup the df for plotting
    df, lgnd_clr = plot_prep(df, unclassified, lgnd_file, top_classes)

    # build the plot
    _ = bracken_plot(df, outfile, lgnd_clmns, lgnd_clr)

    print('\n\nComplete success space cadet! The script without errors.\n\n')

if __name__ == "__main__":
    main()
