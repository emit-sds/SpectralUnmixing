import pandas as pd
import numpy as np
import os
from glob import glob
import geopandas as gp
import struct
import matplotlib.pyplot as plt
from tqdm import tqdm
import json
import isofit.core.common as isc
import multiprocessing as mp
import time
from p_tqdm import p_map
from functools import partial

results = []


def plot_individual_em(wvl, bn, row):
    index = str(row[0])
    spectra = np.asarray(row[1])
    plt.plot(wvl, spectra[6:], label="idx: " + str(index) + ' em :' + spectra[1])
    plt.legend()
    outfname = os.path.join(bn,  spectra[1] + '_' + str(index) + '.png')
    plt.savefig(outfname, bbox_inches='tight')
    plt.clf()
    plt.close()


def parallel_convole(em_spectra, em_wvls, wvls, fwhm, classifications):
    refl_convolve = isc.resample_spectrum(x=em_spectra, wl=em_wvls, wl2=wvls, fwhm2=fwhm)
    result = list(classifications) + list(refl_convolve)
    return (result)


def get_convolved_spectrum(result):
    global results
    results.append(result)


class clean:
    def __init__(self, base_directory: str, configs: str):
        self.base_directory = base_directory
        self.output_directory = os.path.join(self.base_directory, "output")

        # Load configs
        f = open(configs)
        self.data = json.load(f)
        self.metadata = self.data['geo_metadata']  # geodata from spec library
        self.emit_mask = os.path.join(self.base_directory, 'gis', 'emit_mask.shp')  # EMIT dust mask
        self.spectra_files = self.data['spectra_file']  # spectral files
        self.datasets = self.data['datasets']  # spectral datasets
        self.fig_directory = os.path.join(base_directory, "figures")

        # create directory for outputs
        if os.path.isdir(self.output_directory):
            pass
        else:
            os.mkdir(self.output_directory)

        # check for figure directory
        if os.path.isdir(self.fig_directory):
            pass
        else:
            os.mkdir(self.fig_directory)

    def all_data(self):
        # check if directory for all data exists:
        if os.path.isdir(os.path.join(self.output_directory, "all_data")):
            pass
        else:
            os.mkdir(os.path.join(self.output_directory, "all_data"))

        # loop for each dataset
        for i in self.datasets:
            ds_name = os.path.basename(i)
            lat_lon_keys = self.data['lat_lon_keys']
            print("loading... " + ds_name)
            data_master = []

            if ds_name == "ASTER":
                spectral_subset = self.data[i + "_subset"]
                wavelengths = np.linspace(350, 2500, 2151).tolist()
                col_names = ["dataset", "level_1", "level_2", "level_3", 'longitude', 'latitude'] + wavelengths
                for subset in spectral_subset:
                    spectral_files = glob(os.path.join(self.base_directory, 'raw_data', i, subset + '*spectrum.txt'))

                    for file in tqdm(spectral_files, ncols=80, ascii=True, desc=' \t loading...' + subset):
                        name = os.path.basename(file)
                        spec_class = name.split(".")[0]
                        ems_portion = name.split(".")[4]

                        # skip tir data; keep only vswir
                        if ems_portion == 'tir':
                            pass
                        else:
                            # define classifications of spectra
                            if spec_class == 'vegetation':
                                spec_class_1 = 'pv'
                                spec_class_2 = 'pv'
                                spec_class_3 = name.split(".")[2] + "_" + name.split(".")[3]

                            if spec_class == 'nonphotosyntheticvegetation':
                                spec_class_1 = 'npv'
                                spec_class_2 = 'npv'
                                spec_class_3 = name.split(".")[2] + "_" + name.split(".")[3]

                            if spec_class == 'soil':
                                spec_class_1 = 'soil'
                                spec_class_2 = 'soil'
                                spec_class_3 = name.split(".")[1]

                            # load spectra and sort
                            spectra = np.loadtxt(file, skiprows=20)
                            spectra = spectra[spectra[:, 0].argsort()]
                            lat_long = open(file, 'r')
                            lines = lat_long.readlines()
                            long = lines[8].split()[2].replace(";", "")
                            lat = lines[8].split()[1].replace(";", "")

                            if np.min(spectra[:, 0] != 0.35) or spec_class == 'soil':
                                # print(np.round(spectra[1135:1144,0], 4))
                                idx = np.where(np.round(spectra[:, 0], 2) == 2.5)
                                # print(idx[0])
                            else:
                                spectra = (spectra[:, 1] / 100).flatten().tolist()

                                try:
                                    flt_long = float(long)
                                    flt_lat = float(lat)
                                    data_master.append([ds_name, spec_class_1, spec_class_2, spec_class_3, flt_long, flt_lat] + spectra)
                                except:
                                    pass
                # Create master dataframe
                df_master = pd.DataFrame(data_master, columns=col_names)
                df_master.to_csv(os.path.join(self.output_directory, "all_data", "all_data_" + i + ".csv"), index=False)

            elif ds_name == "MEYER_OKIN":
                spectral_table = os.path.join(self.base_directory, 'raw_data', i, self.spectra_files[i])
                df = pd.read_csv(spectral_table, usecols=['file_name', 'description', 'PV', 'NPV', 'LITTER', 'SOIL'])
                coords = self.data['MEYER_GEODATA']
                wavelengths = np.linspace(350, 2500, 2151).tolist()
                meyer_cols = ["dataset", "level_1", "level_2", "level_3", 'longitude', 'latitude'] + wavelengths

                df['em_sum'] = df.PV + df.NPV + df.LITTER + df.SOIL
                df_select = df.loc[df['em_sum'] == 1].copy()

                for index, row in tqdm(df_select.iterrows(), ncols=80, ascii=True, desc='\t loading... ASD files'):
                    file_name = os.path.join(self.base_directory, 'raw_data', i, row[0])
                    spec_class_3 = row[1]

                    if row[2] == 1:
                        spec_class_1 = 'pv'
                        spec_class_2 = 'pv'
                    elif row[3] == 1:
                        spec_class_1 = 'npv'
                        spec_class_2 = 'npv'
                    elif row[4] == 1:
                        spec_class_1 = 'npv'
                        spec_class_2 = 'litter'
                    elif row[5] == 1:
                        spec_class_1 = 'soil'
                        spec_class_2 = 'soil'

                    site_number = row[0][0:3]

                    try:
                        lon_lat = coords[site_number.upper()]

                    except:
                        site_number = site_number[:2]
                        lon_lat = coords[site_number.upper()]

                    data = open(file_name, "rb").read()
                    spectrum = data[484:]
                    spectra = np.array(list(struct.iter_unpack('<f', spectrum)), dtype=float).flatten()

                    # fix first splice of spectrum
                    spectra[:651] *= spectra[651] / spectra[650]

                    spectra = spectra.tolist()
                    data_master.append(
                        [ds_name, spec_class_1, spec_class_2, spec_class_3, lon_lat[0], lon_lat[1]] + spectra)

                # Create master dataframe and filter out bad data
                df_master = pd.DataFrame(data_master, columns=meyer_cols)
                df_master = df_master[df_master.level_3 != 'error']

                # check for figure directory
                if os.path.isdir(os.path.join(self.fig_directory, i)):
                    pass
                else:
                    os.mkdir(os.path.join(self.fig_directory, i))

                base_name = os.path.join(self.fig_directory, i)
                p_map(partial(plot_individual_em, wavelengths, base_name),
                      [row for row in df_master.iterrows()],
                      **{"desc": "\t plotting meyer & okin 2015", "ncols": 150})

                bad_indices = self.data['meyer_bad']
                bad_df = df_master.index.isin(bad_indices)
                df_master = df_master[~bad_df]

                df_master.to_csv((os.path.join(self.output_directory, "all_data", "all_data_" + i + ".csv")),
                                 index=False)

            elif ds_name == 'PILOT-2':
                spectral_table = os.path.join(self.base_directory, 'raw_data', i, self.spectra_files[i])
                col_wls = ["X" + str(x) for x in range(350, 2501)]
                wavelengths = np.linspace(350, 2500, 2151).tolist()
                cols = ['Soil_type_USDA', 'ID', lat_lon_keys[ds_name][1], lat_lon_keys[ds_name][0]] + col_wls
                df = pd.read_csv(spectral_table, usecols=cols)

                # sort order of columns in df
                df = df[cols]

                # rename latitude/longitude
                df = df.rename(columns={lat_lon_keys[ds_name][0]: 'latitude', lat_lon_keys[ds_name][1]: 'longitude'})

                # add level 1 classification
                df.insert(0, 'level_1', 'soil')
                df.insert(0, 'dataset', ds_name)
                col_pilot2 = ["dataset", "level_1", "level_2", "level_3", 'longitude', 'latitude'] + wavelengths
                df.columns = col_pilot2

            elif ds_name == 'OSSL':
                col_wls = ["scan_visnir." + str(x) + '_pcnt' for x in range(350, 2501, 2)]
                meta_cols = ['id.layer_uuid_c', lat_lon_keys[ds_name][1], lat_lon_keys[ds_name][0]]
                spectra_cols = ['id.layer_uuid_c'] + col_wls

                metadata_table = pd.read_csv(os.path.join(self.base_directory, 'raw_data', i, self.metadata[ds_name]),
                                             usecols=meta_cols)
                spectral_table = pd.read_csv(os.path.join(self.base_directory, 'raw_data', i, self.spectra_files[i]),
                                             usecols=spectra_cols)

                # merge dataframes
                df = pd.merge(metadata_table, spectral_table, on='id.layer_uuid_c')
                df = df.rename(columns={lat_lon_keys[ds_name][0]: 'latitude', lat_lon_keys[ds_name][1]: 'longitude'})
                df.insert(0, 'level_1', 'soil')
                df.insert(1, 'level_2', 'soil')
                df.insert(0, 'dataset', ds_name)

                ossl_wvls = [x for x in range(350, 2501, 2)]
                df.columns = ["dataset", "level_1", "level_2", "level_3", 'longitude', 'latitude'] + ossl_wvls

                # convert % to decimal
                for wvl in ossl_wvls:
                    df[wvl] = df[wvl] / 100

            elif ds_name == 'NGSA':
                spectral_table = os.path.join(self.base_directory, 'raw_data', i, self.spectra_files[i])
                df = pd.read_csv(spectral_table, index_col='Wavelength_(nm)')
                df = df.T
                df = df.reset_index()
                df.insert(0, 'level_1', 'soil')
                df.insert(1, 'level_2', 'soil')
                df.insert(0, 'dataset', ds_name)
                df.insert(4, 'longitude', 'unk')
                df.insert(5, 'latitude', 'unk')
                df = df.rename(columns={'index': 'level_3'})

                # save data
            if ds_name == 'MEYER_OKIN' or ds_name == 'ASTER':
                pass
            else:
                df.to_csv((os.path.join(self.output_directory, "all_data", "all_data_" + i + ".csv")), index=False)

    def geo_data(self):
        print("applying geo-filter...", self.emit_mask)
        tables = glob(os.path.join(self.output_directory, "all_data", '*.csv'))
        shp = gp.read_file(self.emit_mask)
        for i in tables:
            ds_name = os.path.basename(i).split(".")[0].split("_")[2]
            df = pd.read_csv(i)
            df = df.drop_duplicates()
            df = df.dropna()

            if ds_name == 'NGSA':
                pass
            else:
                points = gp.GeoDataFrame(df, geometry=gp.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
                within_points = gp.sjoin(points, shp, op='within')

                if within_points.empty:
                    continue
                else:
                    df = pd.DataFrame(within_points)
                    df.drop(columns=df.columns[-14:], axis=1, inplace=True)

            df.to_csv(os.path.join(self.output_directory, 'geofilter', 'geofilter_' + ds_name + '.csv'),
                      index=False)
        print("\t success!")

    def convolve(self, wavelengths: str, geo_filter: bool, level: str):
        # define radiometric resolution for spectra - 4 decimal places ??

        print("loading band convolution...", wavelengths)

        # ask phil about FWHM procedures
        if wavelengths == 'aviris':
            wvls = np.loadtxt(os.path.join('wavelengths', wavelengths + "_wavelengths.txt"), usecols=1) * 10 ** 3
            fwhm = np.loadtxt(os.path.join('wavelengths', wavelengths + "_wavelengths.txt"), usecols=2) * 10 ** 3

        # need FWHM for neon; ask Phil if not calculate
        elif wavelengths == 'neon':
            wvls = np.loadtxt(os.path.join('wavelengths', wavelengths + "_wavelengths.txt"), usecols=0)
            fwhm = np.loadtxt(os.path.join('wavelengths', wavelengths + "_wavelengths.txt"), usecols=1)

        elif wavelengths == 'emit':
            wvls = np.loadtxt(os.path.join('wavelengths', wavelengths + "_wavelengths.txt"), usecols=1) * 10 ** 3
            fwhm = np.loadtxt(os.path.join('wavelengths', wavelengths + "_wavelengths.txt"), usecols=2) * 10 ** 3

        if geo_filter:
            tables = glob(os.path.join(self.output_directory, "geofilter", '*.csv'))

        else:
            tables = glob(os.path.join(self.output_directory, "all_data", '*.csv'))

        spectra_start = 6
        output_cols = ["dataset", "level_1", "level_2", "level_3", 'longitude', 'latitude'] + wvls.tolist()

        # Begin parallel processing
        print("Number of CPUs available:", mp.cpu_count())
        ts = time.time()
        pool = mp.Pool(mp.cpu_count())

        for i in tables:
            df = pd.read_csv(i)
            em_wvls = np.array(df.columns[spectra_start:], dtype=float)
            msg_progess_bar = "\t loading convolution... " + os.path.basename(i).split(".")[0]

            for index, row in tqdm(df.iterrows(), ncols=100, total=df.shape[0], ascii=True, desc=msg_progess_bar):
                classifications = row.values[:spectra_start]
                em_spectra = np.array(row.values[spectra_start:])
                pool.apply_async(parallel_convole, args=(em_spectra, em_wvls, wvls, fwhm, classifications),
                                 callback=get_convolved_spectrum)

        # close the processing tool
        pool.close()
        pool.join()

        # create dataframe, remove duplicates, no data, and save
        df_merge = pd.DataFrame(results, columns=output_cols)
        df_merge = df_merge.drop_duplicates()
        df_merge = df_merge.dropna()

        if geo_filter:
            df_merge.to_csv(os.path.join(self.output_directory, "convolved", "all_data_" + wavelengths + ".csv"),
                            index=False)
        else:
            df_merge.to_csv(os.path.join(self.output_directory, "convolved", "geofilter_" + wavelengths + ".csv"),
                            index=False)

        unique_counts = df_merge[level].value_counts()
        print(unique_counts)
        print('Time in parallel:', time.time() - ts)
        print('\t -Success!')
