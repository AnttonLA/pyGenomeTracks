# -*- coding: utf-8 -*-
from .GenomeTrack import GenomeTrack
import numpy as np
import pandas as pd


# Expects .gwas file


class GwasTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['gwas']  # this is used by make_tracks_file to guess the type of track based on file name
    TRACK_TYPE = 'gwas'
    OPTIONS_TXT = """
# Title of the track. Usually displayed to the right as a label
title =
# Height of the track    
height =
# Type of file
file_type = gwas
# File containing the data. We expect an IGV .gwas format file with the columns: CHR, BP, SNP and P. Optionally, extra
# annotation columns CS and INT can be added.
file =
# Y label text
ylabel =
# Fontsize of the labels
fontsize =
"""
    DEFAULTS_PROPERTIES = {'fontsize': 6,
                           'orientation': None,
                           'color': 'grey',
                           'border_color': 'black',
                           'labels': True,
                           'line_width': 0.5,
                           'max_labels': 60,
                           'max_value': 1,
                           'min_value': 0,
                           'fontstyle': 'normal',
                           'ylabel': None}

    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['title', 'file_type', 'file', 'ylabel']
    FLOAT_PROPERTIES = {'height': [0, np.inf], 'fontsize': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwarg):
        super(GwasTrack, self).__init__(*args, **kwarg)

    def plot(self, ax, chrom, region_start, region_end):
        """
        Plot the given title at a fixed location in the axis.
        The chrom, region_start and region_end variables are not used at the moment.

        :param ax: matplotlib axis to plot
        :param chrom_region: chromosome name
        :param start_region: start coordinate of genomic position
        :param end_region: end coordinate
        """

        f = self.properties['file']
        df = pd.read_csv(f, sep='\t')
        # Check the columns of the df include the required columns (CHR, BP, SNP, P)
        required_columns = ['CHR', 'BP', 'SNP', 'P']
        for c in required_columns:
            if c not in df.columns:
                raise ValueError(f"File {f} does not contain required column {c}")

        x = df['BP']

        def floor(val):
            """
            Alter values to be between 0.02 and 0.95 to avoid cropping by matplotlib.
            TODO: This is a hacky and temporary solution, and does not apply when P is a P-value and not a PP.
            """
            if val > 0.95:
                return 0.95
            elif val < 0.02:
                return 0.02
            else:
                return val

        y = df['P'].apply(floor)

        ax.scatter(x, y, s=10, color='grey')

        ax.set_ylim(bottom=0, top=1)

        # PLOT LEAD VARIANTS, CS COL MARKED AS 1
        if 'CS' in df.columns:
            sub = df[df.CS == 'YES']

            x = sub['BP'].tolist()
            y = sub['P'].apply(floor).tolist()

        if 'INT' in df.columns:  # TODO: this is temporary. Needs to be far more robust.
            # Names will be a list with '' in the positions where INT is not YES and the value of the SNP column otherwise.
            names = sub.apply(lambda row: row['SNP'] if row['INT'] == 'YES' else '', axis=1).tolist()

        ax.scatter(x, y, s=40, color='red', marker='h')

        for i, n in enumerate(names):
            xy = (x[i], y[i])
            ax.text(xy[0], xy[1], n, fontsize=12, ha='center', va='bottom')

    def plot_y_axis(self, ax, plot_axis, transform='no', log_pseudocount=0, y_axis='transformed', only_at_ticks=False):
        """
        Plot the y-axis of the track. This overwrites the GenomeTrack.plot_y_axis method.
        """
        # Add y axis label
        # TODO: this does not work. The idea is to add a label to the y axis so that we can specify the unit.
        # Currently, I am trying to plot the title of the track, but it does not work:(
        GenomeTrack.plot_y_axis(self, ax, plot_axis, transform, log_pseudocount, y_axis, only_at_ticks, add_ylabel=True,
                                ylabel_text=self.properties['ylabel'])

    #def plot_label(self, label_ax, width_dpi, h_align='left'):
    #    pass
