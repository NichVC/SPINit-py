# Official packages
import copy
import errno
import glob
import os
import xml.etree.ElementTree as ET

from typing import List, Dict

# Third-party packages
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

# In-house packages


FILE_PATH_TEMPLATE = {
    "data"  : None,
    "Serie" : None,
    "header": None
}

class SPINitExp(object):
    """
        Basic Class that that read, parses, and stores parameters and data acquired from the SPINit-SpinAligner platform.

        Parameters:
        -----------

        exp_data_path: str
            experimental data path

    """





    
    def __init__(self, exp_dataset_path:str, **kwargs) -> None:
        
        """
        0. verify exp_dataset_path
            confirm the existence of:
            i.  data.dat    : binary data file
            ii. Serie.xml   : parameter file
            iii.header.xml  : parameter file
        1. load parameters from files
        2. load data from file
        
        """
        self._file_paths = self._verify_dataset(exp_dataset_path)

        self.params = self._load_params()
        self.data = self._load_data()

    def _verify_dataset(self, exp_dataset_path):
        file_paths = FILE_PATH_TEMPLATE

        if (os.path.exists(exp_dataset_path)):
            for key, val in file_paths.items():
                for file in glob.glob(exp_dataset_path + '*' + key + '.*'):
                    if file:
                        file_paths[key] = file
                    else:
                        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"No <{key}> file in {exp_dataset_path}")
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"<{exp_dataset_path}> does not exist")
            
        return file_paths

    def _load_params(self)->Dict:
        
        Serie_path = os.path.join(self._file_paths['Serie'])
        Serie_dict = self._generate_param_dict(Serie_path)

        header_path = os.path.join(self._file_paths['header'])
        header_dict = self._generate_param_dict(header_path)


        return (Serie_dict | header_dict)
    
    def _generate_param_dict(self, xml_path:str)->Dict:
        xml = ET.parse(xml_path)
        root = xml.getroot()

        elems = root.findall(".//params//entry")
        params = [ list(elem)[0].text for elem in list(elems)]
        matching_exprs = [".//params//entry[key='" + param + "']//value//value" for param in params]
        vals  = [root.findtext(expr) for expr in matching_exprs]        

        return dict(zip(params, vals))


    def _load_data(self):
        fid = np.fromfile(self._file_paths['data'], dtype='>f' )
        fid = np.asarray(fid[0::2, ...] + 1j * fid[1::2, ...])
        fid.astype(np.complex128)

        dim_1d = int(self.params['ACQUISITION_MATRIX_DIMENSION_1D'])
        dim_2d = int(self.params['ACQUISITION_MATRIX_DIMENSION_2D'])
        dim_3d = int(self.params['ACQUISITION_MATRIX_DIMENSION_3D'])
        dim_4d = int(self.params['ACQUISITION_MATRIX_DIMENSION_4D'])

        fid = np.reshape(fid, (dim_4d, dim_3d, dim_2d, dim_1d))
        fid = np.squeeze(fid)
        return fid

    def plot_data_magnitude(self):
        data_extremum = np.max(np.abs(self.data))
        plt.figure()
        plt.imshow(np.abs(self.data), vmax=data_extremum, vmin= -1*data_extremum, cmap='seismic')
        plt.show()