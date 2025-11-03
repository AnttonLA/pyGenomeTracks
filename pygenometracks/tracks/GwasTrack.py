# -*- coding: utf-8 -*-
from .GenomeTrack import GenomeTrack
import numpy as np
import \
    pandas as pd  # TODO: package is not otherwise pandas dependent (uses custom scripts to read data in). Should probably do the same.
import matplotlib.pyplot as plt  # Matplotlib is already a dependency of pyGenomeTracks. We need it to get the cmap.


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
# Color of accentuated points (CS = 1)
color =
"""
    DEFAULTS_PROPERTIES = {'fontsize': 6,
                           'color': 'red',
                           'border_color': 'black',
                           'labels': True,
                           'line_width': 0.5,
                           'max_labels': 60,
                           'max_value': 1,
                           'min_value': 0,
                           'fontstyle': 'normal',
                           'ylabel': None,
                           'y_values_format': 'pval',
                           'y_axis_max_val': None,
                           'id_fontsize': 12,
                           'id_position': 'center',
                           'cs_dotsize': 45,
                           'ytick_fontsize': 8}

    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {'y_values_format': ['PP', '-log10', 'pval'], 'id_position': ['center', 'left', 'right']}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['title', 'file_type', 'file', 'ylabel', 'y_values_format', 'color']
    FLOAT_PROPERTIES = {'height': [0, np.inf], 'fontsize': [0, np.inf], 'id_fontsize': [0, np.inf],
                        'cs_dotsize': [0, np.inf], 'y_axis_max_val': [0, np.inf]}
    INTEGER_PROPERTIES = {'ytick_fontsize': [0, 100]}

    def __init__(self, *args, **kwarg):
        super(GwasTrack, self).__init__(*args, **kwarg)

    def process_gwas(self):
        """
        Process input .gwas file to a pandas dataframe.
        Read the file, check it has all the required columns. Add CS and INT columns if not present.
        If y_values_format is -log10, transform the P-values to -log10(P).

        :return: pandas dataframe with the data, and the maximum value of the y-axis
        """
        f = self.properties['file']
        df = pd.read_csv(f, sep='\t')
        # Check the columns of the df include the required columns (CHR, BP, SNP, P)
        required_columns = ['CHR', 'BP', 'SNP', 'P']
        for column in required_columns:
            if column not in df.columns:
                raise ValueError(f"File {f} does not contain required column {column}")

        # Add columns for CS and INT if they don't exist (standard IGV .gwas format)
        if 'CS' not in df.columns:
            df['CS'] = '0'  # TODO: these should be set to 0 when we make the change to 0/1 instead of NO/YES. /done now
        if 'INT' not in df.columns:
            df['INT'] = 0


        # For the -log10 scale, calculate the -log10 of the P-value
        if self.properties['y_values_format'] == '-log10':
            df['P'] = df['P'].apply(lambda v: -np.log10(v))

        # Maximum value of the y-axis, rounded up
        maxval = np.ceil(df['P'].max())

        return df, maxval

    def plot(self, ax, chrom, region_start, region_end):
        """
        Plot the given title at a fixed location in the axis.
        The chrom, region_start and region_end variables are not used at the moment.

        :param ax: matplotlib axis to plot
        :param chrom_region: chromosome name
        :param start_region: start coordinate of genomic position
        :param end_region: end coordinate
        """

        df, max_y = self.process_gwas()
        if self.properties['y_axis_max_val']:
            max_y = self.properties['y_axis_max_val']

        x = df['BP']

        def clamp_pp_vals(val):
            """
            Alter values to be between 0.02 and 0.95 to avoid cropping by matplotlib.
            TODO: This is a hacky and temporary solution, and does not apply when P is a P-value and not a PP.
            """
            if val > 0.95:
                return 0.95
            elif val < max_y / 20:
                return max_y / 20
            else:
                return val

        if self.properties['y_values_format'] == 'PP':
            y = df['P'].apply(clamp_pp_vals)
            if max_y > 1:
                print("Y axis cannot be bigger than 1 for Posterior Probabilities!")
                max_y = 1
            ax.set_ylim(bottom=0, top=max_y)
        elif self.properties['y_values_format'] == '-log10':
            y = df['P']  # Values are already -log10 transformed by process_gwas()
            ax.set_ylim(bottom=0, top=max_y)
        elif self.properties['y_values_format'] == 'pval':
            print(self.properties['y_values_format'])
            y = df['P']
        else:
            raise ValueError(f"y_values_format {self.properties['y_values_format']} not recognized.")

        ax.scatter(x, y, s=10, color='grey')  # The color of the non-CS points is grey

        # Plot the Credible Set variants on top of the grey points
        if 'CS' in df.columns:
            cs_values = df[df['CS'] != 0]['CS'].unique()
            cs_values = sorted(cs_values)

            cmap = plt.get_cmap('tab10')  # change this to 'Set1', 'tab20' or whatever as needed

            for i, cs_val in enumerate(cs_values):
                cs_subset = df[df['CS'] == cs_val]
                x_subset = cs_subset['BP']
                y_subset = cs_subset['P'].apply(clamp_pp_vals) if self.properties['y_values_format'] == 'PP' else cs_subset['P']
                color = cmap(i % cmap.N)  # Wrap around if there are more cs_values than colors in the colormap

                ax.scatter(x_subset,
                           y_subset,
                           s=self.properties['cs_dotsize'],
                           facecolor=color,
                           marker='o',
                           edgecolors='black',
                           linewidths=.66)

                # Label the INT=1 variants in that credible set
                if 'INT' in cs_subset.columns:
                    int_subset = cs_subset[cs_subset['INT'] != 0]
                    for _, row in int_subset.iterrows():
                        ax.text(row['BP'],
                                row['P'].apply(clamp_pp_vals) if self.properties['y_values_format'] == 'PP' else row['P'] + row['INT'],
                                row['SNP'],
                                fontsize=self.properties['id_fontsize'],
                                ha=self.properties['id_position'],
                                va='bottom',
                                snap=True)

    def plot_y_axis(self, ax, plot_axis, transform='no', log_pseudocount=0, y_axis='transformed', only_at_ticks=False,
                    add_ylabel=True, ylabel_text=None, ytick_fontsize=8):
        """
        Override the GenomeTrack.plot_y_axis method to add a y-axis label to the normally label-less y-axis.
        """
        GenomeTrack.plot_y_axis(self, ax, plot_axis, transform, log_pseudocount, y_axis, only_at_ticks, add_ylabel=True,
                                ylabel_text=self.properties['ylabel'], ytick_fontsize=self.properties['ytick_fontsize'])

    # def plot_label(self, label_ax, width_dpi, h_align='left'):
    #    pass
