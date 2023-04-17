# -*- coding: utf-8 -*-
"""
This module will load in 2D data from the SPINit software and integrate center peaks. 
This is meant for use in analyzing hyporpolarization 13C build-up and frequency 
sweep curves from a SpinAligner setup.

Created on Fri Jun 10 10:05:24 2022
Version 1.0 - Latest update: 24-11-2022

@author: Nichlas Vous Christensen
@email: nvc@clin.au.dk
@phone: +45 23464522
@organization: Aarhus University, The MR Research Centre
"""

from xml.etree.ElementTree import parse
from struct import unpack
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
import os
from datetime import datetime

def load_SpinAligner_data(dir_path):
    # Function to search and save parameters from the XML metadata documents.
    def get_header_param(xml, param_name):
        expr = "./params/entry[key='" + param_name + "']/value/value"
        return xml.getroot().findtext(expr)
    def get_serie_param(xml, param_name):
        expr = ".//commands/[name='Polarization frequency sweep']/parameters[name='" + param_name + "']/value"
        return xml.getroot().findtext(expr)
    
    header_path = dir_path + '/header.xml' 
    serie_path = dir_path + '/Serie.xml'
    data_path = dir_path + '/data.dat' 
    
    date_time = datetime.fromtimestamp(os.path.getmtime(dir_path + '/data.dat'))
    
    header_xml = parse(header_path)
    serie_xml = parse(serie_path)
    
    dimension1d = int(get_header_param(header_xml, 'MATRIX_DIMENSION_1D')) 
    dimension2d = int(get_header_param(header_xml, 'MATRIX_DIMENSION_2D')) 
    measure_time = float(get_header_param(header_xml, 'Polarisation_Growth_Delay'))
    freq_initial = get_serie_param(serie_xml, 'Initial frequency')
    freq_step_size = get_serie_param(serie_xml, 'Frequency step')
    
    # If it is a frequency sweep, generate frequency step array.
    if freq_initial != None and freq_step_size != None:
        freq_steps = (np.repeat(float(freq_initial),dimension2d)+np.cumsum(np.repeat(float(freq_step_size),dimension2d)))-float(freq_step_size)
    else:
        freq_steps = None
        
    data_raw = np.zeros((dimension2d,dimension1d),dtype=complex)
    
    # Open data file and save real and imaginary 2D data.
    data_file = open(data_path, mode='rb') 
    for i in range(0, dimension2d): 
        for j in range(0, dimension1d): 
            real, imaginary = unpack('>ff', data_file.read(8)) 
            data_raw[i,j] = complex(real,imaginary)
                         
    # Fourier transform and roll peak to center of spectrum.
    midpoint = int(data_raw.shape[1]/2)
    data = fft(data_raw,axis=1)
    data = np.roll(data,midpoint)
    
    # Generate timelist 't' from 'measure_time' parameter.
    t = np.cumsum(np.repeat(measure_time,dimension2d))
    
    # All spectra are integrated around the center point.
    # For 512 spectral points, 103 points around the center are integrated.
    # From testing on pyruvate data, this fits approximately with the peak base.
    integral_range = [int(midpoint-midpoint/5), int(midpoint+midpoint/5)]
    integrals = abs(data)[:,integral_range[0]:integral_range[1]].sum(axis=1)
    
    return data, integrals, t, freq_steps, date_time

if __name__ == '__main__':
    # Data folder which includes data.dat, header.xml, serie.xml files.
    # Example files: "1033" (frequency sweep) and "1106" (build-up monitering)
    dir_path = 'data_example\\1033'
    [data, integrals, t, freq_steps, date_time] = load_SpinAligner_data(dir_path)
    
    # Plots for testing:
    #plt.plot(abs(data[13,:])) 
    if freq_steps is None:
        plt.plot(t,integrals)
        plt.xlabel('Time (s)')
        plt.ylabel('Integral (arb. unit)')      
    else:
        plt.plot(freq_steps/1e9,integrals)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Integral (arb. unit)')    

           












