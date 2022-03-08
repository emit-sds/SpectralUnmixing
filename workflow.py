import pandas as pd
import numpy as np
import os
import argparse
import clean_spectra
import build_libraries
import download_data
from figures import figures
from run_unmix import unmmix


def dwnld(base_directory: str, configs: str, wavelengths: str):
    dd = download_data.create_environment(base_directory=base_directory, configs=configs)
    dd.download()


def run_cs(base_directory:str, configs:str, wavelengths:str, geofilter: bool, level: str):
    cs = clean_spectra.clean(base_directory=base_directory, configs=configs)
    cs.all_data()
    cs.geo_data()
    cs.convolve(wavelengths=wavelengths, geo_filter=geofilter, level=level)


def build_em(base_directory:str, configs:str, wavelengths:str, cols:int, level:str, combinations:int):
    bl = build_libraries.endmembers(base_directory=base_directory, configs=configs, instrument=wavelengths)
    bl.percentage_base()
    bl.reflectance_bootstrap(cols=cols, level=level, combinations=combinations)


# setting the dry-run to False will execute the subprocess call
def spectral_unmix(base_directory:str):
    su = unmmix(base_directory=base_directory)
    su.debug(dry_run=True)


def create_figures(base_directory: str, wavelengths:str):
    figs = figures(base_directory=base_directory, instrument=wavelengths)
    figs.endmember_library()
    figs.individual_em()


def main():
    parser = argparse.ArgumentParser(description='Run spectra clean workflow')
    parser.add_argument('-bd', type=str, help='Specify base directory')
    parser.add_argument('-configs', type=str, help='Specify configuration file', default='configs.json')
    parser.add_argument('-wvls', type=str, help='Specify instrument wavelengths', default='emit')
    parser.add_argument('-dwnld', type=bool, help='Download data?', default=False)
    parser.add_argument('-geo', type=bool, help='Apply geo-filter with EMIT shapefile', default=True)
    parser.add_argument('-mode', type=str, help='Specify mode to run')
    parser.add_argument('-cols', type=int, help='# of columns in reflectance files', default=5)
    parser.add_argument('-level', type=str, help='level of classification to use', default='level_1')
    parser.add_argument('-comb', type=int, help='number of spectral combinations to use eg bootsrap', default=50000)

    args = parser.parse_args()

    if args.dwnld:
        dwnld(base_directory=args.bd, configs=args.configs, wavelengths=args.wvls)

    if args.mode in ['process', 'all']:
        run_cs(base_directory=args.bd, configs=args.configs, wavelengths=args.wvls, geofilter=args.geo,
               level=args.level)

    if args.mode in ['build', 'all']:
        build_em(base_directory=args.bd, configs=args.configs, wavelengths=args.wvls,
                 cols=args.cols, level=args.level, combinations=args.comb)

    if args.mode in ['unmix', 'all']:
        spectral_unmix(base_directory=args.bd)

    if args.mode in ['figs', 'all']:
        create_figures(base_directory=args.bd, wavelengths=args.wvls)


if __name__ == '__main__':
    main()
