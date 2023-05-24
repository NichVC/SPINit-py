"""
Created on Tue May 23 11:07:10 2023

@author: Xiao JI
@email: xiao.ji@ucsf.edu
@phone: +1 (919)-564-6986
@organization: Department of Radiology, UCSF

@author: Nichlas Vous Christensen
@email: nvc@clin.au.dk
@phone: +45 23464522
@organization: Aarhus University, The MR Research Centre
"""

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

from scipy.optimize import curve_fit

import warnings
  
# In-house packages

DATA_PROCESSING_PARAMETERS = {
    'start_timing_at_zero'          : False,
    'signal_processing_method'      : None,
    'signal_processing_mode'        : "Magnitude",
    'frequency_domain_upper_bound'  : None,
    'frequency_domain_lower_bound'  : None,
}

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
        0. update parameters for post-processing        
        1. verify exp_dataset_path
            confirm the existence of:
            i.  data.dat    : binary data file
            ii. Serie.xml   : parameter file
            iii.header.xml  : parameter file
        2. load parameters from files
        3. load data from file
        
        """
        self.data_processing_params = self._update_data_processing_params(kwargs)
        self._file_paths = self._verify_dataset(exp_dataset_path)

        self.params = self._load_params()
        self.data = self._load_data()

    def _update_data_processing_params(self, kwargs):
        """
        parse possible tags for post-processing
        """
        recon_params = copy.deepcopy(DATA_PROCESSING_PARAMETERS)
        recon_params.update((k, kwargs[k]) for k in (recon_params.keys() & kwargs.keys()) )
        return recon_params  

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
        
        Serie_path = self._file_paths['Serie']
        Serie_dict = self._generate_param_dict(Serie_path)

        header_path = self._file_paths['header']
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

    def plot_data_magnitude_2d(self):
        data_extremum = np.max(np.abs(self.data))
        plt.figure()
        plt.imshow(np.abs(self.data), vmax=data_extremum, vmin= -1*data_extremum, cmap='seismic')
        plt.show()

    def _contract_f1_gbe(self):
        """
        Contract the direct dimension (F1) of the FID using <global background estimation> algorithm
        """
        return NotImplemented
    
    def _contract_f1_fds(self):
        """
        Contract the direct dimension (F1) of the FID using <frequency_domain_summation> algorithm
        """
        return NotImplemented    

    def _contract_f1_hfs(self):
        """
        Contract the direct dimension (F1) of the FID using <half_fid_substraction> algorithm

        Under the assuption that the noise power spectrum is stable along the whole FID, 
            the total noise energy in the first half of the FID is the same w.r.t the second half.
        When the T2* is so short and the FID is so long that the observable signal dies out in the first half of the signal, 
            the difference of the Frobius-norm between the first half and second half of the FID is then a good metric that descibes the signal intensity.
        """
        pos_half_f1 = int(np.ceil(int(self.params['ACQUISITION_MATRIX_DIMENSION_1D'])/2))
        return np.sum(np.abs(self.data[...,:pos_half_f1:]), axis=1) - np.sum(np.abs(self.data[...,pos_half_f1::]), axis=1)

    def _contract_f1_dfs(self):
        """
        Contract the direct dimension (F1) of the FID using <direct_fid_summation> algorithm
        """
        if self.data_processing_params['signal_processing_mode'] == 'Magnitude':
            return np.sum(np.abs(self.data), axis=1)
        elif self.data_processing_params['signal_processing_mode'] == 'Real':
            return np.sum(np.real(self.data), axis=1)
        elif self.data_processing_params['signal_processing_mode'] == 'Imaginary':
            return np.sum(np.imag(self.data), axis=1)
        else:
            raise ValueError("signal_process_mode can only be chose from the three following options: Magnitude, Real, Imaginary")

    def _mono_exp_growth(self,x,A,B,C):
        return A * (1 - np.exp(-x*B)) + C

    def _bi_exp_growth(self,x,A,B,C,D,E):
        return A * (1 - np.exp(-x*B)) + C * (1 - np.exp(-x*D)) + E
    
    def _mono_exp_decay(self,x,A,B,C):
        return A * np.exp(-x*B) + C

    def _bi_exp_decay(self,x,A,B,C,D,E):
        return A * (np.exp(-x*B)) + C * (np.exp(-x*D)) + E


class SPINitSpectrum(SPINitExp):
    """
    Container class for simple 1D dataset
    """

    def __init__(self, exp_dataset_path:str, **kwargs):
        super().__init__(exp_dataset_path, **kwargs)
        self.spectrum = self.FT_data()

    def FT_data(self):
        return np.fft.fftshift(np.fft.fft(self.data))
    
class SPINitEvolution(SPINitExp):
    """
    Container class for time dependent multi-dimensional experiments such as solid state build-ups
    The first dimension (outer-most) is the dimension of time-evolution
    """

    def __init__(self, exp_dataset_path:str, **kwargs):
        super().__init__(exp_dataset_path, **kwargs)
        self.bup_curve = self.process_data()
        self.fit_data = self.fit_data

    def process_data(self):
        delta_time = float(self.params['Polarisation_Growth_Delay'])
        
        if self.data_processing_params['signal_processing_method']:
            if self.data_processing_params['signal_processing_method'] == 'global_background_estimation':
                bup_intensity_array = self._contract_f1_gbe()
            if self.data_processing_params['signal_processing_method'] == 'frequency_domain_summation':
                bup_intensity_array = self._contract_f1_fds()
            if self.data_processing_params['signal_processing_method'] == 'half_fid_substraction':
                bup_intensity_array = self._contract_f1_hfs()                
            if self.data_processing_params['signal_processing_method'] == 'direct_fid_summation':
                bup_intensity_array = self._contract_f1_dfs()                
        else:
            bup_intensity_array = self._contract_f1_dfs()
        
        bup_intensity_array = bup_intensity_array[bup_intensity_array != 0]

        dimension_f2 = len(bup_intensity_array) if len(bup_intensity_array) < int(self.params['ACQUISITION_MATRIX_DIMENSION_2D']) else int(self.params['ACQUISITION_MATRIX_DIMENSION_2D'])
        
        bup_time_array = np.cumsum(np.repeat(float(delta_time), dimension_f2))
        if self.data_processing_params['start_timing_at_zero']:
            bup_time_array = bup_time_array - delta_time

        return (bup_time_array, bup_intensity_array)
    
    def fit_data(self,fit_model = 'Mono_Exp_Growth', plot_bup_fit = True):
        
        self.fit_model = fit_model
        
        # Supressing warning in fitting due to exp() of big numbers.
        warnings.filterwarnings('ignore')
        
        # Prefiltering of data in order to remove wrongly captured data (after dissolution).
        bup_time_array, bup_intensity_array = self.bup_curve
        
        bup_max = np.max(bup_intensity_array)
        numb_points_to_remove = np.sum(bup_intensity_array[np.argmax(bup_intensity_array):]<0.5*bup_max)
        bup_intensity_array = bup_intensity_array[:-numb_points_to_remove]
        bup_time_array = bup_time_array[:-numb_points_to_remove]
        
        # Perform fitting and return fitted curve and parameters.
        if self.fit_model == 'Mono_Exp_Growth':
            popt, pcov = curve_fit(self._mono_exp_growth, bup_time_array, bup_intensity_array,p0=[np.abs(bup_max),1e-3,np.abs(bup_max)/100])
            bup_intensity_array_fitted = self._mono_exp_growth(bup_time_array,popt[0],popt[1],popt[2])
        elif self.fit_model == 'Bi_Exp_Growth':
            popt, pcov = curve_fit(self._bi_exp_growth, bup_time_array, bup_intensity_array,p0=[np.abs(bup_max),1e-3,np.abs(bup_max)/100,1e-3,np.abs(bup_max)/100])
            bup_intensity_array_fitted = self._bi_exp_growth(bup_time_array,popt[0],popt[1],popt[2],popt[3],popt[4])   
        else:
            raise ValueError("fit_model can only be chose from the two following options: Mono_Exp_Growth, Bi_Exp_Growth")
            
        self.fit_params = popt
        self.fit_curves = (bup_time_array, bup_intensity_array, bup_intensity_array_fitted)
            
        if plot_bup_fit:
            self._plot_bup_fit()
        
        return (bup_time_array, bup_intensity_array, bup_intensity_array_fitted, popt)
    
    def _plot_bup_fit(self):
        plt.figure()
        if self.fit_model == 'Mono_Exp_Growth':
            half_recovery_m = int(1/self.fit_params[1]/60) # minutes
            half_recovery_s = round(1/self.fit_params[1]%60) # seconds
            plt.plot(self.fit_curves[0]/60,self.fit_curves[1],'b*')
            plt.plot(self.fit_curves[0]/60,self.fit_curves[2],'k')
            plt.xlabel('Time (min)')
            plt.ylabel('Integral (arb. unit)')    
            plt.grid(True)
            plt.legend(['Datapoints',f'Fit: {self.fit_params[0]:.0f} $\cdot$ (1 - exp(-x $\cdot$ {self.fit_params[1]:.6f})) + {self.fit_params[2]:.0f}'])
            plt.title(f'Half recovery: {half_recovery_m:.0f} min and {half_recovery_s:.0f} seconds\nMax polarization = {max(self.fit_curves[2]):.0f}')
        elif self.fit_model == 'Bi_Exp_Growth':
            half_recovery_m_1 = int(1/self.fit_params[1]/60) # minutes
            half_recovery_s_1 = round(1/self.fit_params[1]%60) # seconds
            half_recovery_m_2 = int(1/self.fit_params[3]/60) # minutes
            half_recovery_s_2 = round(1/self.fit_params[3]%60) # seconds
            plt.plot(self.fit_curves[0]/60,self.fit_curves[1],'b*')
            plt.plot(self.fit_curves[0]/60,self.fit_curves[2],'k')
            plt.xlabel('Time (min)')
            plt.ylabel('Integral (arb. unit)')    
            plt.grid(True)
            plt.legend(['Datapoints',f'Fit: {self.fit_params[0]:.0f} $\cdot$ (1 - exp(-x $\cdot$ {self.fit_params[1]:.6f}))\n + {self.fit_params[2]:.0f} $\cdot$ (1 - exp(-x $\cdot$ {self.fit_params[3]:.6f})) + {self.fit_params[4]:.0f}'])
            plt.title(f'Half recovery [1]: {half_recovery_m_1:.0f} min and {half_recovery_s_1:.0f} seconds\nHalf recovery [2]: {half_recovery_m_2:.0f} min and {half_recovery_s_2:.0f} seconds\nMax polarization = {max(self.fit_curves[2]):.0f}')
       

    

class SPINitSweep(SPINitExp):
    """
    Container class for multi-dimensional experiments with parameter sweeps, such as Microwave Frequency/Power sweeps
    The first dimension (outer-most) is the dimension of parameter sweep
    """
    def __init__(self, exp_dataset_path:str, **kwargs):
        super().__init__(exp_dataset_path, **kwargs)
        self.sweep_params = self._load_sweep_params()
        self.sweep_profile = self.process_data()

    def _load_sweep_params(self):
        xml = ET.parse(self._file_paths['Serie'])
        root = xml.getroot()
        elems = root.findall(".//setUpQueue//commands//parameters//name")
        params = [ elem.text for elem in list(elems)]
        elems = root.findall(".//setUpQueue//commands//parameters//value")
        vals  = [ elem.text for elem in list(elems)]
        return dict(zip(params, vals))

    def process_data(self):
        sweep_initial    = self.sweep_params['Initial frequency']
        sweep_step_size  = self.sweep_params['Frequency step']
        
        if self.data_processing_params['signal_processing_method']:
            if self.data_processing_params['signal_processing_method'] == 'global_background_estimation':
                sweep_intensity_array = self._contract_f1_gbe()
            if self.data_processing_params['signal_processing_method'] == 'frequency_domain_summation':
                sweep_intensity_array = self._contract_f1_fds()
            if self.data_processing_params['signal_processing_method'] == 'half_fid_substraction':
                sweep_intensity_array = self._contract_f1_hfs()                
            if self.data_processing_params['signal_processing_method'] == 'direct_fid_summation':
                sweep_intensity_array = self._contract_f1_dfs()                
        else:
            sweep_intensity_array = self._contract_f1_dfs()
        
        sweep_intensity_array = sweep_intensity_array[sweep_intensity_array != 0]

        dimension_f2 = len(sweep_intensity_array) if len(sweep_intensity_array) < int(self.params['ACQUISITION_MATRIX_DIMENSION_2D']) else int(self.params['ACQUISITION_MATRIX_DIMENSION_2D'])
        
        sweep_param_array = (np.repeat(float(sweep_initial), dimension_f2)+np.cumsum(np.repeat(float(sweep_step_size), dimension_f2)))-float(sweep_step_size)

        return (sweep_param_array, sweep_intensity_array)

