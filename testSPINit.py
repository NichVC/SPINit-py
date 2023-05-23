# -*- coding: utf-8 -*-
"""
Created on Tue May 23 11:07:10 2023

@author: Nichlas Vous Christensen
@email: nvc@clin.au.dk
@phone: +45 23464522
@organization: Aarhus University, The MR Research Centre
"""

import numpy as np
import matplotlib.pyplot as plt

from SPINitExp import SPINitEvolution, SPINitSweep

dataset_path_1106 = "./data_example/1106/"

dataset_1106 = SPINitEvolution(dataset_path_1106, signal_processing_mode="Magnitude",fitting_method='mono')

plt.figure()
time_pts, intensity_pts = dataset_1106.bup_curve
plt.scatter(x=time_pts,y=intensity_pts)
time_pts_fit, intensity_pts, intensity_pts_fit, popt = dataset_1106.fit_model
plt.plot(time_pts_fit,intensity_pts_fit,linewidth=3,color='k')
plt.show()














