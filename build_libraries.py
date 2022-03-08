import os
import json
import pandas as pd
from sklearn.model_selection import train_test_split
import itertools
import numpy as np
from p_tqdm import p_map
from functools import partial
from spectral.io import envi
import time
from tqdm import tqdm


def reflectance_row(data):

    row_spectra = []
    row_fractions = []
    for seed in data[1]:
        np.random.seed(seed)
        fractions = np.random.dirichlet(np.ones(3))

        spectra = np.array(data[0][0]).astype(dtype=float) * fractions[0] + \
                  np.array(data[0][1]).astype(dtype=float) * fractions[1] + \
                  np.array(data[0][2]).astype(dtype=float) * fractions[2]
        row_spectra.append(spectra)
        row_fractions.append(fractions)

    return row_spectra, row_fractions


def get_meta(lines: int, samples: int, bands, wvls: bool):

    if wvls:
        meta = {
            'lines': lines,
            'samples': samples,
            'bands': len(bands),
            'wavelength': bands,
            'interleave': 'bil',
            'data type': 4,
            'file_type': 'ENVI Standard',
            'byte order': 0,
            'header offset': 0,
            'wavelength units': 'nm'
        }
    else:
        meta = {
            'lines': lines,
            'samples': samples,
            'bands': len(bands),
            'interleave': 'bil',
            'data type': 4,
            'file_type': 'ENVI Standard',
            'byte order': 0,
            'header offset': 0
        }

    return meta


def save_envi(output_file, meta, grid):
    outDataset = envi.create_image(output_file, meta, ext='', force=True)
    mm = outDataset.open_memmap(interleave='bip', writable=True)
    mm[...] = grid
    del mm


class endmembers:
    def __init__(self, base_directory: str, configs: str, instrument: str):
        self.base_directory = base_directory
        self.output_directory = os.path.join(self.base_directory, "output")

        # Load configs
        f = open(configs)
        self.data = json.load(f)
        self.instrument = instrument

        # define seed
        np.random.seed(13)

    def percentage_base(self):
        """ Create a analysis (endmember) and simulation library..."""
        print('Splitting library...')
        table = os.path.join(self.output_directory, "convolved", 'all_data_' + self.instrument + '.csv')
        df = pd.read_csv(table)
        ems = df['level_1'].unique()

        analysis_list = []
        simulation_list = []

        for _em, em in enumerate(ems):
            df_filter = df.loc[df['level_1'] == em].copy()
            simulation, analysis = train_test_split(df_filter, test_size=0.4, random_state=13)
            analysis_list.append(analysis)
            simulation_list.append(simulation)

        df_analysis = pd.concat(analysis_list)
        df_simulation = pd.concat(simulation_list)

        # clean the spectral tables with bad indices
        for _df_spec, df_spec in enumerate([df_analysis, df_simulation]):
            if _df_spec == 0:
                out_name = 'endmember_library.csv'
            else:
                out_name = 'simulation_library.csv'

            df_spec.to_csv(os.path.join(self.output_directory, out_name), index=False)

        print('\t success!')

    def reflectance_bootstrap(self, cols: int, level: str, combinations:int):
        ts = time.time()
        simulation_table = os.path.join(self.base_directory, 'output', 'simulation_library.csv')
        df = pd.read_csv(simulation_table)

        class_names = df.level_1.unique()
        class_lists = []

        for em in class_names:
            df_select = df.loc[df[level] == em].copy()
            df_select = df_select.values.tolist()
            class_lists.append(df_select)

        all_combinations = list(itertools.product(*class_lists))
        bootsrap_index = np.random.choice(len(all_combinations), replace=False, size=combinations)
        bootstrap_spectra = [all_combinations[i] for i in bootsrap_index]

        # Import wavelengths and genereate grids
        if self.instrument == 'aviris':
            wvls = np.loadtxt(os.path.join('wavelengths', self.instrument + "_wavelengths.txt"), usecols=1) * 10 ** 3

        # need FWHM for neon; ask Phil if not calculate
        elif self.instrument == 'neon':
            wvls = np.loadtxt(os.path.join('wavelengths', self.instrument + "_wavelengths.txt"), usecols=0)

        elif self.instrument == 'emit':
            wvls = np.loadtxt(os.path.join('wavelengths', self.instrument + "_wavelengths.txt"), usecols=1) * 10 ** 3

        fraction_grid = np.zeros((len(bootsrap_index), cols, len(class_names)))
        spectra_grid = np.zeros((len(bootsrap_index), cols, len(wvls)))
        index_grid = np.zeros((len(bootsrap_index), len(class_names), len(wvls)))

        # spectra array
        start = 0
        spec_array = np.array(bootstrap_spectra)

        # Process row in parallel
        seeds = list(range(0, len(bootsrap_index) * cols))
        seeds_array = np.asarray(seeds)
        seeds_array = seeds_array.reshape(len(bootsrap_index), cols)

        # parallel spectra processes
        process_spectra = p_map(reflectance_row, [(spec_array[_row, :, 5:], seeds_array[_row, :]) for _row, row in
                                                  enumerate(spectra_grid)], **{"desc": "\t\t generating fractions...",
                                                                               "ncols": 150})

        # Populate results row by row
        for _row, row in enumerate(tqdm(process_spectra, ncols=150, desc='Building reflectance file...')):
            for _col, (refl, frac) in enumerate(zip(row[0], row[1])):
                spectra_grid[_row, _col, :] = refl
                fraction_grid[_row, _col, :] = frac

            index_grid[_row, :, :] = spec_array[_row, :, 5:].astype(dtype=float)

        # save the datasets
        refl_meta = get_meta(lines=len(bootsrap_index), samples=cols, bands=wvls, wvls=True)
        index_meta = get_meta(lines=len(bootsrap_index), samples=len(class_names), bands=wvls, wvls=False)
        fraction_meta = get_meta(lines=len(bootsrap_index), samples=cols, bands=class_names, wvls=False)

        # save index, spectra, fraction grid
        output_files = [os.path.join(self.output_directory, 'index_grid_bootstrap.hdr'),
                        os.path.join(self.output_directory, 'spectra_grid_bootstrap.hdr'),
                        os.path.join(self.output_directory, 'fraction_grid_bootstrap.hdr')]

        meta_docs = [index_meta, refl_meta, fraction_meta]
        grids = [index_grid, spectra_grid, fraction_grid]

        p_map(save_envi, output_files, meta_docs, grids, **{"desc": "\t\t saving envi files...",
                                                            "ncols": 150})
        print('Time in parallel:', time.time() - ts)
        print('\t -Success!')
