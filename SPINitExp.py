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
    
    def __init__(self, exp_dataset_path:str, **kwargs) -> None:
        
        """
        
        """

