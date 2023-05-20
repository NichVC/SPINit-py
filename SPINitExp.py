# Official packages
import os
import copy

from typing import List, Dict

# Third-party packages
import numpy as np

# In-house packages

class SPINitExp(object):
    """
        Basic Class that that read, parses, and stores parameters and data acquired from the SPINit-SpinAligner platform.

        Parameters:
        -----------

        exp_data_path: str
            experimental data path

    """

    FILE_PATH_TEMPLATE = {
        "data"  = "data.dat",
        "Serie" = "Serie.xml",
        "header"= "header.xml"
    }



    
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
        return NotImplemented

    def _load_params(self):
        return NotImplemented

    def _load_data(self):
        return NotImplemented

