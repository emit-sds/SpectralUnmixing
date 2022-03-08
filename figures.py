import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sys import platform
import os
from p_tqdm import p_map
from functools import partial

if not "win" in platform:
    plt.switch_backend('Agg')

# bad wavelength regions
bad_wv_regions = [[0, 440], [1310, 1490], [1770, 2050], [2440, 2880]]


def get_good_bands_mask(wavelengths, wavelength_pairs):
    wavelengths = np.array(wavelengths)
    if wavelength_pairs is None:
        wavelength_pairs = bad_wv_regions
    good_bands = np.ones(len(wavelengths)).astype(bool)

    for wvp in wavelength_pairs:
        wvl_diff = wavelengths - wvp[0]
        wvl_diff[wvl_diff < 0] = np.nanmax(wvl_diff)
        lower_index = np.nanargmin(wvl_diff)

        wvl_diff = wvp[1] - wavelengths
        wvl_diff[wvl_diff < 0] = np.nanmax(wvl_diff)
        upper_index = np.nanargmin(wvl_diff)
        good_bands[lower_index:upper_index + 1] = False
    return good_bands


def plot_individual_em(wvl, em_name, bn, row):
    index = str(row[0])
    spectra = np.asarray(row[1])
    plt.plot(wvl, spectra[6:], label=em_name + ": " + str(index) + ' lib :' + spectra[0])
    plt.legend()
    outfname = os.path.join(bn, em_name + '_' + str(index) + '.png')
    plt.savefig(outfname, bbox_inches='tight')
    plt.clf()
    plt.close()


class figures:
    def __init__(self, base_directory: str, instrument: str):
        self.base_directory = base_directory
        self.output_directory = os.path.join(base_directory, 'output')
        self.fig_directory = os.path.join(base_directory, "figures")

        # check for figure directory
        if os.path.isdir(self.fig_directory):
            pass
        else:
            os.mkdir(self.fig_directory)

        # data paths - libraries
        self.em_lib = os.path.join(self.output_directory, 'endmember_library.csv')
        self.sim_lib = os.path.join(self.output_directory, 'simulation_library.csv')

        # set wavelegnths
        if instrument == 'aviris':
            self.wvls = np.loadtxt(os.path.join('wavelengths', instrument + "_wavelengths.txt"), usecols=1) * 10 ** 3
        elif instrument == 'neon':
            self.wvls = np.loadtxt(os.path.join('wavelengths', instrument + "_wavelengths.txt"), usecols=0)
        elif instrument == 'emit':
            self.wvls = np.loadtxt(os.path.join('wavelengths', instrument + "_wavelengths.txt"), usecols=1) * 10 ** 3

    def endmember_library(self):
        print('building endmember libraries...')
        df_em = pd.read_csv(self.em_lib)
        df_sim = pd.read_csv(self.sim_lib)

        dfs = [df_em, df_sim]
        wvls = self.wvls
        good_bands = get_good_bands_mask(wvls, bad_wv_regions)
        wvls[~good_bands] = np.nan

        for _df, df in enumerate(dfs):
            ems = df.level_1.unique()

            # create figure
            fig = plt.figure(constrained_layout=True, figsize=(6, 6))
            ncols = 1
            nrows = len(ems)
            gs = gridspec.GridSpec(ncols=ncols, nrows=nrows, wspace=0, hspace=0, figure=fig)

            for _em, em in enumerate(ems):
                df_select = df.loc[df['level_1'] == em].copy()

                # add plots
                ax = fig.add_subplot(gs[_em, 0])
                ax.set_ylim(0, 1)
                ax.set_xticks(np.arange(0, 2500, 100), minor=True)
                ax.set_yticks(np.arange(0, 1, 0.1), minor=True)
                ax.set_xlim(np.nanmin(self.wvls), np.nanmax(self.wvls))

                ax.set_ylabel('Reflectance (%)')
                if _em == 2:
                    ax.set_xlabel('Wavelength (nm)')

                for _row, row in df_select.iterrows():
                    spectra = row[6:]
                    ax.plot(wvls, spectra)

            if _df == 0:
                out_name = 'analysis'
            else:
                out_name = 'simulation'

            # save figure
            plt.savefig(os.path.join(self.fig_directory, out_name + '.png'), dpi=300, bbox_inches='tight')
            plt.clf()
            plt.close()
        print("\t success!")

    def individual_em(self):
        print("plotting individual endmembers...")
        df_em = pd.read_csv(self.em_lib)
        df_sim = pd.read_csv(self.sim_lib)

        dfs = [df_em, df_sim]
        wvls = self.wvls
        good_bands = get_good_bands_mask(wvls, bad_wv_regions)
        wvls[~good_bands] = np.nan

        for _df, df in enumerate(dfs):
            ems = df.level_1.unique()

            if _df == 0:
                out_name = 'analysis'
            else:
                out_name = 'simulation'

            # check for directory
            if os.path.isdir(os.path.join(self.fig_directory, out_name)):
                pass
            else:
                os.mkdir(os.path.join(self.fig_directory, out_name))

            # loop through spectra
            for _em, em in enumerate(ems):
                df_select = df.loc[df['level_1'] == em].copy()
                # check for directory
                if os.path.isdir(os.path.join(self.fig_directory, out_name, em)):
                    pass
                else:
                    os.mkdir(os.path.join(self.fig_directory, out_name, em))

                base_name = os.path.join(self.fig_directory, out_name, em)

                p_map(partial(plot_individual_em, wvls, em, base_name),
                      [row for row in df_select.iterrows()],
                      **{"desc": "\t plotting individual em: " + em, "ncols": 150})

        print('\t\t success!')
