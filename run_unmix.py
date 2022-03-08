import os
from sys import platform
import subprocess


def execute_call(cmd_list, dry_run=False):
    if dry_run:
        print(cmd_list)
    else:
        subprocess.call(cmd_list)


class unmmix:
    def __init__(self, base_directory: str):
        self.base_directory = base_directory
        self.em_file = os.path.join(self.base_directory, 'output', 'endmember_library.csv')
        self.reflectance_file = os.path.join(self.base_directory, 'output', 'spectra_grid_bootsrap.hdr')

        # create results directory
        for i in ["sma", "mesma", "debug"]:
            if os.path.isdir(os.path.join(self.base_directory, "output", i)):
                pass
            else:
                os.mkdir(os.path.join(self.base_directory, "output", i))

        # get platform for processing
        if "win" in platform:
            self.level_arg = ['Level_1']
            self.n_cores = '15'
        else:
            self.level_arg = 'Level_1'
            self.n_cores = '40'

    #def sma(self, normalization=True):
    #def mesma(self, normalization=True):

    def debug(self, dry_run=True):
        execute_call(['julia', '-p', self.n_cores, 'unmix.jl', self.reflectance_file, self.em_file, self.level_arg,
                      os.path.join(self.base_directory, 'output', 'debug', "debug_test"), "--mode", "sma",
                      "--num_endmembers", "-1", "--normalization", "brightness"], dry_run)




        # run options for sma/mesma
        # normalization_options = ["brightness", "none", "1070", "1756", "2030", "1500"]
        # num_mc = ["5", "10", "50", "100", "200"]
        # num_em = ["5", "10", "30"]
        # max_comb = ["10", "100", "500", "1000"]
        #
        # dry_run = False




