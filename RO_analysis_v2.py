#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 10:37:47 2024

@author: alexandriaudenkwo
"""
#Key
#QB is Q1
#QA is Q2

import numpy as np
import matplotlib.pyplot as plt
import fit_functions as fitfun
import analysis
from scipy.stats import norm

importdir = r"/Users/alexandriaudenkwo/Desktop/Research/two_qubit_engine/"
  
#Analyze RO sweep data

#Data here is raw data that needs to be analyzed; example 20 Mhz
raw_data_A = np.loadtxt(importdir+"20_mhz/integratpn_sweep2_20_mhz_detuning_piA_piB_ro_B.txtA", float)
raw_data_B = np.loadtxt(importdir+"20_mhz/integratpn_sweep2_20_mhz_detuning_piA_piB_ro_B.txtB", float)

#We actually have the data for all detunings so below is 13 MHz
#raw_data_A = np.loadtxt(importdir+"13_mhz/integratpn_sweep4_13_mhz_detuning_piA_piB_ro_B.txtA", float)
#raw_data_B = np.loadtxt(importdir+"13_mhz/integratpn_sweep4_13_mhz_detuning_piA_piB_ro_B.txtB", float)


start_RO = 100
stop_RO = 4100
steps_RO = 41
delta = (stop_RO - start_RO)/(steps_RO - 1)

#dictionary to store axes for each plot
d = {}

binwidth =20

#for loop to analyze RO data from 100 to 4100
for i in range(steps_RO):
    #Make axes for each plot
    plt.figure(i, figsize=(8, 8))
    d["ax_hist"+str(i)] = plt.axes()
    d["ax_hist"+str(i)].tick_params(direction='in', labelbottom=False)
    d["ax_hist"+str(i)].set_title('Pi_A Pi_B RO_B RO duration: ' +str(start_RO+i*delta))
    
    #Define bins
    pi_a_min = np.min(raw_data_A[i])
    pi_a_max = np.max(raw_data_A[i])
    bins_pi_a = np.arange(pi_a_min, pi_a_max, binwidth)
    
    pi_b_min = np.min(raw_data_B[i])
    pi_b_max = np.max(raw_data_B[i])
    bins_pi_b = np.arange(pi_b_min, pi_b_max, binwidth)
    
   
    #Histogram for pi on QA
    hist_pi_a= d["ax_hist"+str(i)].hist(raw_data_A[i], bins=bins_pi_a, histtype='step', orientation='vertical', color = "firebrick",density = True)
    data_A = np.sort(raw_data_A[i])
    
    #Histogram for pi on QB
    hist_pi_b= d["ax_hist"+str(i)].hist(raw_data_B[i], bins=bins_pi_b, histtype='step', orientation='vertical', color = "mediumpurple", density = True)
    data_B = np.sort(raw_data_B[i])
    
    #two lines of code below from ChaatGPT; centers the bin so that bin array and count array are of same dimension
    bin_centers_pi_a = [(bins_pi_a[j] + bins_pi_a[j+1]) / 2 for j in range(len(bins_pi_a) - 1)]
    bin_centers_pi_b = [(bins_pi_b[j] + bins_pi_b[j+1]) / 2 for j in range(len(bins_pi_b) - 1)]
    
    #fit
    gaussian_pi_a_vals,hist_fit_a = analysis.fit_gaussian(data_A, bin_centers_pi_a, hist_pi_a[0])
    gaussian_pi_b_vals,hist_fit_b = analysis.fit_gaussian(data_B, bin_centers_pi_b, hist_pi_b[0])

    #plot result
    plt.plot(bin_centers_pi_a, hist_fit_a,'k', linewidth=2)
    plt.plot(bin_centers_pi_b, hist_fit_b,'k', linewidth=2)
    plt.show()
