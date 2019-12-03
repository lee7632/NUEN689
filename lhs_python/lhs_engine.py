from scipy.stats import norm
import random
import numpy as np
#
# 
#
param_dict = dict()
#
param_dict['sp_angle'] = dict()
param_dict['inlet_L']  = dict()
param_dict['Vx_scale_1'] = dict()
param_dict['Vy_scale_1'] = dict()
param_dict['Vz_scale_1'] = dict()
param_dict['k_scale_1']  = dict()
param_dict['Vx_scale_2'] = dict()
param_dict['Vy_scale_2'] = dict()
param_dict['Vz_scale_2'] = dict()
param_dict['k_scale_2']  = dict()
param_dict['nu']       = dict()
param_dict['CmuL']     = dict()
param_dict['beta']     = dict()
#
#
param_dict['sp_angle']['cdf_funcs'] = norm(3.0,0.3)
param_dict['inlet_L']['cdf_funcs']  = norm(50,5)
param_dict['Vx_scale_1']['cdf_funcs'] = norm(1.0, 0.2e-2)
param_dict['Vy_scale_1']['cdf_funcs'] = norm(1.0, 0.2e-2)
param_dict['Vz_scale_1']['cdf_funcs'] = norm(1.0, 0.2e-2)
param_dict['k_scale_1']['cdf_funcs']  = norm(1.0, 5.e-2)
param_dict['Vx_scale_2']['cdf_funcs'] = norm(1.0, 0.2e-2)
param_dict['Vy_scale_2']['cdf_funcs'] = norm(1.0, 0.2e-2)
param_dict['Vz_scale_2']['cdf_funcs'] = norm(1.0, 0.2e-2)
param_dict['k_scale_2']['cdf_funcs']  = norm(1.0, 5.e-2)
param_dict['nu']['cdf_funcs']       = norm(1.0e-6,5e-2*1.0e-6)
param_dict['CmuL']['cdf_funcs']     = norm(0.000315, 7.009251029889001e-05)
param_dict['beta']['cdf_funcs']     = norm(0.075, 0.02)
#
param_dict['sp_angle']['values'] = []
param_dict['inlet_L']['values']  = []
param_dict['Vx_scale_1']['values'] = []
param_dict['Vy_scale_1']['values'] = []
param_dict['Vz_scale_1']['values'] = []
param_dict['k_scale_1']['values']  = []
param_dict['Vx_scale_2']['values'] = []
param_dict['Vy_scale_2']['values'] = []
param_dict['Vz_scale_2']['values'] = []
param_dict['k_scale_2']['values']  = []
param_dict['nu']['values']       = []
param_dict['CmuL']['values']     = []
param_dict['beta']['values']     = []
#
bins = 20
for key in param_dict.keys():
    sample = random.sample(range(bins), bins)
    param_dict[key]['sampling_integers'] = sample.copy()
#
#
dx = 1.0/bins
#
for i in range(bins):
    model_values = []
    for key in param_dict.keys():
        p_bin = param_dict[key]['sampling_integers'][i]
        x_min = dx*p_bin
        x_max = x_min + dx
        while(True):
            ran_val = np.random.random()
            if ran_val>=x_min and ran_val<=x_max:
                break
        #
        pdf_sample = param_dict[key]['cdf_funcs'].ppf(ran_val)
        #
        param_dict[key]['values'].append(pdf_sample)
        #
#
#

