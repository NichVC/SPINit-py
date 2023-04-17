# -*- coding: utf-8 -*-
"""
This module will load in 2D data from the SPINit software and save it according to 
data type (build-up monitoring / frequency sweep) as a report.

Created on Fri Jun 10 13:03:43 2022
Version 1.0 - Latest update: 24-11-2022

@author: Nichlas Vous Christensen
@email: nvc@clin.au.dk
@phone: +45 23464522
@organization: Aarhus University, The MR Research Centre
"""

from load_SpinAligner_data import load_SpinAligner_data
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import tkinter.filedialog
import tkinter
import time

def make_plot(dir_path,save_dir,sample):
    # Generate integrals
    [data, integrals, t, freq_steps, date_time] = load_SpinAligner_data(dir_path)
    
    figure = plt.figure(figsize=[10,6])
    # Build-up curve or frequency sweep.
    if freq_steps is None:
        # If present, remove wrongly captured data at the end (after sample has been dissolved).
        max_index = np.argmax(integrals)
        cut_off_index = len(data)-1
        for i in range(max_index+1,len(integrals)):
            if integrals[max_index] * 0.9 > integrals[i]:
                cut_off_index = i
                break
        integrals = integrals[0:cut_off_index-1]
        t = t[0:cut_off_index-1]
        
        # Do exponential fit.
        def func(x,A,B,C):
            return A * (1 - np.exp(-x*B)) + C
        
        popt, pcov = curve_fit(func, t, integrals,p0=[100e6,1e-3,30e5])
        half_recovery_m = int(1/popt[1]/60) # minutes
        half_recovery_s = round(1/popt[1]%60) # seconds
        
        plt.plot(t/60,integrals,'b*')
        plt.plot(t/60,func(t,popt[0],popt[1],popt[2]),'k')
        plt.xlabel('Time (min)')
        plt.ylabel('Integral (arb. unit)')    
        plt.grid(True)
        plt.legend(['Datapoints',f'Fit: {popt[0]:.0f} $\cdot$ (1 - exp(-x $\cdot$ {popt[1]:.6f})) + {popt[2]:.0f}'])
        plt.title(f'Sample: {sample}\nHalf recovery: {half_recovery_m:.0f} min and {half_recovery_s:.0f} seconds\nMax polarization = {max(integrals):.0f}')
        figure_type = 'Build-up-curve'
    else:  
        plt.plot(freq_steps/1e9,integrals,'b*')
        plt.plot(freq_steps/1e9,integrals,'b-')
        plt.xlabel('Transmitter Frequency (GHz)')
        plt.ylabel('Integral (arb. unit)') 
        plt.title(f'Sample: {sample}')
        plt.grid(True)
        figure_type = 'Frequency-sweep'

    # Save plot as image file in the correct folder.
    year = date_time.year
    month = date_time.strftime("%b")
    month_num = date_time.month
    day = date_time.day
    time_stamp = date_time.strftime("%Y%m%d-%H%M%S")
    if os.path.isdir(f'{save_dir}\\{year}_reports\\{month_num}-{month}\\{day}'):
        pass
    else:
        os.makedirs(f'{save_dir}\\{year}_reports\\{month_num}-{month}\\{day}')
        
    figure.savefig(f'{save_dir}\\{year}_reports\\{month_num}-{month}\\{day}\\{time_stamp}_{figure_type}.png', format='png')

def main():
    samples = ['18 ul C1-pyr', '150 ul C1-pyr', 'C2-pyr', 'bicarbonate', 'fumerate']
    bioprobe = int(input('\nBioprobe?\n18 ul pyr [0], 150 ul pyr [1], alanine [2], bicarbonate [3], fumerate [4], custom sample name [5]: '))
    if bioprobe != 5:
        sample = samples[0]
    else:
        sample = input('Write custom sample name: ')
    
    current_dir = os.getcwd()
    
    root = tkinter.Tk()
    root.withdraw()
    dir_path = tkinter.filedialog.askdirectory(initialdir=f"{current_dir}\\data_example",title='Select a data folder')
    save_dir = current_dir
    
    make_plot(dir_path,save_dir,sample)
    
    print('\nReport has been generated, program will exit in 5 seconds.')
    time.sleep(5)

if __name__ == "__main__":
    main()


