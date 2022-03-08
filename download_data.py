import os
from tqdm import tqdm
import requests
import json

class create_environment:
    def __init__(self, base_directory: str, configs: str):
        self.base_directory = base_directory
        self.output_directory = os.path.join(self.base_directory, "output")

        # Load configs
        f = open(configs)
        self.data = json.load(f)
        self.metadata = self.data['geo_metadata']  # geodata from spec library
        self.spectra_files = self.data['spectra_file']  # spectral files
        self.datasets = self.data['datasets']  # spectral datasets

        # create directories to download data
        if os.path.isdir(os.path.join(self.base_directory, 'raw_data')):
            pass
        else:
            os.mkdir(os.path.join(self.output_directory, 'raw_data'))


    def download(self):
        print("downloading data...")
        urls = self.data['urls']
        for data_set in self.datasets:
            url = urls[data_set]

            # create directories to download data
            if os.path.isdir(os.path.join(self.base_directory, 'raw_data', data_set)):
                pass
            else:
                os.mkdir(os.path.join(self.base_directory, 'raw_data', data_set))

            if url == 'UNK':
                pass

            else:
                r = requests.get(url, stream=True)

                file_size = int(r.headers.get('Content-Length', 0))

                local_filename = url.split('/')[-1]
                download_location = os.path.join(self.base_directory, 'raw_data', data_set, local_filename)
                desc = '\t downloading... ' + data_set

                with open(download_location, "wb") as f, tqdm(desc=desc,
                                                              total=file_size, unit='iB',
                                                              unit_scale=True,
                                                              unit_divisor=1024, leave=True, ncols=100
                                                              ) as bar:
                    for data in r.iter_content(chunk_size=1024):
                        size = f.write(data)
                        bar.update(size)

                # will add function to unzip file and keep only the .csv tables

        print("data ready for processing!")