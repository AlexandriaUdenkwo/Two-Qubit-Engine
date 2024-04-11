# -*- coding: utf-8 -*-
"""
Created on Fri May 27 15:14:48 2022

@author: crow104
"""

import numpy as np
import sys
file_dir= r"C:\Users\crow104\Documents\Python Scripts\sequence_generator"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Users\crow104\Documents\Python Scripts\data_acquisition_scripts"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Data\2020\201021_thermal_readout\Ramsey vs thermal noise"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Users\crow104\Documents\Python Scripts\sequence_generator\py_sequences"
if file_dir not in sys.path: sys.path.append(file_dir)

from generator import *
import os
pi = np.pi
import pyvisa
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
from tkinter import Tk
import tkinter as tki
from tkinter.filedialog import askopenfilename, asksaveasfilename, asksaveasfile
import wx_programs
import scipy.io as io
import time
from IPython import get_ipython
from scipy import optimize as opt
from matplotlib.patches import Rectangle
import warnings
from scipy.optimize import fsolve

import sys
import numpy as np
import matplotlib.pyplot as plt
import FPJPA_sequence as Fs
import analysis
import IQ_cluster
import bnc
import math

instrument_path = r"C:\Users\Crow108\Documents\Python\instr\analyzer"
if instrument_path not in sys.path: sys.path.append(instrument_path )

analyzer_path = r"C:\Users\Crow108\Documents\Python\instr\python_interface\python_without_WX2184C"
if analyzer_path not in sys.path: sys.path.append(analyzer_path )
import daq_programs

save_dir = r"C:\Data\2022\2022-10-17_SQUILL_NH_4Kfilter_BBN_JPA/"
target_bnc_black_address = 'USB0::0x03EB::0xAFFF::621-016100000-0126::INSTR'
target_bnc_address =  'USB0::0x03EB::0xAFFF::141-216340000-0292::INSTR'

steps = 61
start_val = 0#6.888895-0.0006 #8.77
stop_val = 360#6.888895+0.0006 #12.77
sweep_vals = np.linspace(start_val,stop_val,steps)

reps = 1000
ro_dur = 3000
pi_ge = 51#52
pi_ef = 34#33
ssm_ge = -0.126#-0.125#-0.2006
ssm_ef = ssm_ge - 0.15
wx_amps = [1, 1, 1, 1]

snr_ge_vals = []
snr_ef_vals = []
snr_gf_vals = []
dist_ge_vals = []
dist_ef_vals = []
dist_gf_vals = []

state_assignment = np.zeros((1000,3))

#
#for i,ss in enumerate(sweep_vals):
#    print(i)
##    print(i*(stop_val-start_val)/(steps-1)+start_val)
#    print(ss)
##    ro_dur = np.int(ss)
#    IQ_angle = ss
##    bnc.set_bnc_output(ss, bnc_addr=target_bnc_black_address) #in GHz
##    bnc.set_bnc_output(12.7159,ss,bnc_addr=target_bnc_address) #in dBm for power
##    
#    ##run nopi, pi_ge, pi_ge + pi_ef sequences
#    Fs.pi_nopi(off = 0, coef=0,coefpief=0,pi_ge=pi_ge,pi_ef=pi_ef,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
#    wx_programs.wx_set_and_amplitude_and_offset(amp=wx_amps)
#    daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_nopi,counts_nopi,rec_readout_nopi,rec_avg_all_nopi,rec_allnopi = daq_programs.run_daq(3,reps,ro_dur,IQangle=IQ_angle); #np.int(ss)
#    Fs.pi_nopi(off = 0,coef=1,coefpief=0,pi_ge=pi_ge,pi_ef=pi_ef,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
#    wx_programs.wx_set_and_amplitude_and_offset(amp=wx_amps)
#    daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_pi,counts_pi,rec_readout_pi,rec_avg_all_pi,rec_allpi = daq_programs.run_daq(3,reps,ro_dur,IQangle=IQ_angle);
#    Fs.pi_nopi(off = 0,coef=1,coefpief=1,pi_ge=pi_ge,pi_ef=pi_ef,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
#    wx_programs.wx_set_and_amplitude_and_offset(amp=wx_amps)
#    daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_pipi,counts_pipi,rec_readout_pipi,rec_avg_all_pipi,rec_allpipi = daq_programs.run_daq(3,reps,ro_dur,IQangle=IQ_angle);
##    
##    for i in range(3):
##        for j in range(reps):
##            state,pair_angle = daq_programs.IQ_angle_threshold(IQ_center_x = -80,IQ_center_y =-120,angle_1=260,angle_2=0,angle_3=100,IQ_x=rec_readout_vs_pats[0,j,i],IQ_y =rec_readout_vs_pats[1,j,i])
##            state_assignment[j,i] = state
##            
##            ##TO DO - reassign IQ points to respective blobs from state_assignment and scatter plot to see assignment result
##    g_x = []
##    g_y = []
##    e_x = []
##    e_y = []
##    f_x = []
##    f_y = []
##   
##    for i in range(3):
##       for j in range(reps):
##           if state_assignment[j,i] == 0:
##               g_x.append(rec_readout_vs_pats[0,j,i])   #rec_readout_vs_pats[IorQ,# of reps, seq point (g e or f)]
##               g_y.append(rec_readout_vs_pats[1,j,i])
##               
##           if state_assignment[j,i] == 1:
##               e_x.append(rec_readout_vs_pats[0,j,i])   #rec_readout_vs_pats[IorQ,# of reps, seq point (g e or f)]
##               e_y.append(rec_readout_vs_pats[1,j,i])
##               
##           if state_assignment[j,i] == 2:
##               f_x.append(rec_readout_vs_pats[0,j,i])   #rec_readout_vs_pats[IorQ,# of reps, seq point (g e or f)]
##               f_y.append(rec_readout_vs_pats[1,j,i])
##    ##if state_assignment[j][i] == 0
##    ###g_x.append(rec_readout_vs_pats[0,j,i])
##    ###g_y.append(rec_readout_vs_pats[1,j,i])
##    ##if state_assignment[j][i] == 1
##    ###e_x.append(rec_readout_vs_pats[0,j,i])
##    ###e_y.append(rec_readout_vs_pats[1,j,i])
##    ##if state_assignment[j][i] == 2
##    ###f_x.append(rec_readout_vs_pats[0,j,i])
##    ###f_y.append(rec_readout_vs_pats[1,j,i])
###    plt.show()
##            
##            
#    ##calculate center of blobs
#    g_blob_x = sum(rec_readout_nopi[0])/len(rec_readout_nopi[0])
#    g_blob_y = sum(rec_readout_nopi[1])/len(rec_readout_nopi[1])
#    e_blob_x = sum(rec_readout_pi[0])/len(rec_readout_pi[0])
#    e_blob_y = sum(rec_readout_pi[1])/len(rec_readout_pi[1])
#    f_blob_x = sum(rec_readout_pipi[0])/len(rec_readout_pipi[0])
#    f_blob_y = sum(rec_readout_pipi[1])/len(rec_readout_pipi[1])
#
#    ##calculate difference between blobs
#    ge_blob_distance = math.sqrt((g_blob_x-e_blob_x)**2+(g_blob_y-e_blob_y)**2)
#    ef_blob_distance = math.sqrt((f_blob_x-e_blob_x)**2+(f_blob_y-e_blob_y)**2)
#    gf_blob_distance = math.sqrt((g_blob_x-f_blob_x)**2+(g_blob_y-f_blob_y)**2)
#    
#    dist_ge_vals.append(ge_blob_distance)
#    dist_ef_vals.append(ef_blob_distance)
#    dist_gf_vals.append(gf_blob_distance)
#    
#    print("dist_ge:",ge_blob_distance)
#    print("dist_ef:",ef_blob_distance)
#    print("dist_gf:",gf_blob_distance)
##    
##    ##fit to gaussian and calculate SNR
#    fit_pi, counts_fit_pi = analysis.fit_single_gaussian(bins_pi,counts_pi)
#    fit_pipi, counts_fit_pipi = analysis.fit_single_gaussian(bins_pipi,counts_pipi)
#    fit_nopi, counts_fit_nopi = analysis.fit_single_gaussian(bins_nopi,counts_nopi)
#    snr_ge = ((fit_pi[0][1] - fit_nopi[0][1])/(np.abs(fit_pi[0][2]) + np.abs(fit_nopi[0][2]))*2)**2
#    snr_ef = ((fit_pi[0][1] - fit_pipi[0][1])/(np.abs(fit_pi[0][2]) + np.abs(fit_pipi[0][2]))*2)**2
#    snr_gf = ((fit_pipi[0][1] - fit_nopi[0][1])/(np.abs(fit_pipi[0][2]) + np.abs(fit_nopi[0][2]))*2)**2
#    
#    snr_ge_vals.append(snr_ge)
#    snr_ef_vals.append(snr_ef)
#    snr_gf_vals.append(snr_gf)
#    
#    print("snr_ge:",snr_ge)
#    print("snr_ef:",snr_ef)
#    print("snr_gf:",snr_gf)
##    
##    
###    np.savetxt(save_dir+'sweep_RO_freq_7.0137_by100kHz_SNR_q1',snr_vals)
###    np.savetxt(save_dir+'Dec19_sweep_IQ_angle_SNR_with15dBmamp_q1',snr_vals)
###    np.savetxt(save_dir+'Dec19_sweep_RO_dur500to3000_SNR_with15dBmamp_q1',snr_vals)
###    np.savetxt(save_dir+'Dec21_sweep_JPA8.77to12.77dBm_SNR_with15dBmamp_q1',snr_vals)
##    
##    
##    plt.plot(bins_pi,counts_pi,label='pi')
##    plt.plot(bins_nopi,counts_nopi,label='nopi')
##    plt.plot(bins_pipi,counts_pipi,label='pipi')
##    plt.plot(bins_pi,counts_fit_pi)
##    plt.plot(bins_nopi,counts_fit_nopi)
##    plt.legend()
##    plt.show()
##    
##    ##plot blobs on same plot
###    plt.scatter(rec_readout_pipi[0][0:1000],rec_readout_pipi[1][0:1000],c='black', marker='.', alpha=0.1)
###    plt.scatter(rec_readout_pipi[0][0:1000],rec_readout_pipi[1][0:1000],c='red', marker='.', alpha=0.1)
###    plt.scatter(rec_readout_pipi[0][0:1000],rec_readout_pipi[1][0:1000],c='blue', marker='.', alpha=0.1)
#    plt.show()
#    plt.scatter(rec_readout_nopi[0][0:1000],rec_readout_nopi[1][0:1000],c='black', marker='.', alpha=0.4)
#    plt.scatter(rec_readout_pi[0][0:1000],rec_readout_pi[1][0:1000],c='red', marker='.', alpha=0.4)
#    plt.scatter(rec_readout_pipi[0][0:1000],rec_readout_pipi[1][0:1000],c='blue', marker='.', alpha=0.4)
#    plt.show()
##    
##    plt.scatter(rec_readout_vs_pats[0,0:1000,0],rec_readout_vs_pats[1,0:1000,0],c='black', marker='.', alpha=0.1,label='|g>')
##    plt.scatter(rec_readout_vs_pats[0,0:1000,1],rec_readout_vs_pats[1,0:1000,1],c='red', marker='.', alpha=0.1,label='|e>')
##    plt.scatter(rec_readout_vs_pats[0,0:1000,2],rec_readout_vs_pats[1,0:1000,2],c='blue', marker='.', alpha=0.1,label='|f>')
##    plt.legend()
##    
##    plt.show()
##    
##    
##
##    plt.scatter(g_x,g_y,c='black', marker='.', alpha=0.1,label='|g>')
##    plt.scatter(e_x,e_y,c='red', marker='.', alpha=0.1,label='|e>')
##    plt.scatter(f_x,f_y,c='blue', marker='.', alpha=0.1,label='|f>')
##    plt.legend()
#    
#
#### plot sweep results ###
#plt.plot(sweep_vals,snr_ge_vals,label='SNR_ge')
#plt.plot(sweep_vals,snr_ef_vals,label='SNR_ef')
#plt.plot(sweep_vals,snr_gf_vals,label='SNR_gf')
#
#SNRgemax = np.amax(snr_ge_vals)
#sweep_SNRgemax = sweep_vals[np.where(snr_ge_vals == SNRgemax)]
#plt.plot(sweep_SNRgemax,SNRgemax,'*',label='SNR_ge max')
#SNRefmax = np.amax(snr_ef_vals)
#sweep_SNRefmax = sweep_vals[np.where(snr_ef_vals == SNRefmax)]
#plt.plot(sweep_SNRefmax,SNRefmax,'*',label='SNR_ef max')
#SNRgfmax = np.amax(snr_gf_vals)
#sweep_SNRgfmax = sweep_vals[np.where(snr_gf_vals == SNRgfmax)]
#plt.plot(sweep_SNRgfmax,SNRgfmax,'*',label='SNR_gf max')
#
#plt.xlabel('IQ angle (degrees)')
#plt.xlabel('readout duration (ns)')
#plt.xlabel('readout frequency (GHz)')
##plt.xlabel('JPA power (dBm)')
#plt.ylabel('SNR')
#plt.legend()
#plt.show()
##
##
#plt.plot(sweep_vals,dist_ge_vals,label='dist_ge')
#plt.plot(sweep_vals,dist_ef_vals,label='dist_ef')
#plt.plot(sweep_vals,dist_gf_vals,label='dist_gf')
#
#distgemax = np.amax(dist_ge_vals)
#sweep_distgemax = sweep_vals[np.where(dist_ge_vals == distgemax)]
#plt.plot(sweep_distgemax,distgemax,'*',label='dist_ge max')
#distefmax = np.amax(dist_ef_vals)
#sweep_distefmax = sweep_vals[np.where(dist_ef_vals == distefmax)]
#plt.plot(sweep_distefmax,distefmax,'*',label='dist_ef max')
#distgfmax = np.amax(dist_gf_vals)
#sweep_distgfmax = sweep_vals[np.where(dist_gf_vals == distgfmax)]
#plt.plot(sweep_distgfmax,distgfmax,'*',label='dist_gf max')
##
##plt.xlabel('IQ angle (degrees)')
##plt.xlabel('readout duration (ns)')
#plt.xlabel('readout frequency (GHz)')
###plt.xlabel('JPA power (dBm)')
#plt.ylabel('SNR (blob distance)')
#plt.legend()
#    
    
#    trainx = np.concatenate((rec_readout_pi[0],rec_readout_nopi[0]))
#    trainy = np.concatenate((rec_readout_pi[1],rec_readout_nopi[1]))
#    y = IQ_cluster.kMeans(trainx,trainy,plot=True)
#    
#    IQ_cluster.svm_class(trainx,trainy, y, C = 0.05, plot=True)
    
### IQ CLUSTERING CODE

def cluster_model(ssm_ge,ssm_ef,pi_ge,pi_ef,ro_dur,IQ_angle):
    rec_g = [[]]*2;
    rec_e = [[]]*2;
    rec_f = [[]]*2;
    rec_read = [[]]*2;
    c1 = 0
    c2 = 0 
     
    
    for i in range(3):
        
        if i == 0:
            c1 = 0
            c2 = 0
        elif i == 1:
            c1 = 1
            c2 = 0
        elif i == 2:
            c1 = 1
            c2 = 1
        Fs.pi_nopi(coef=c1,coefpief=c2,pi_ge=pi_ge,pi_ef=pi_ef,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
        wx_programs.wx_set_and_amplitude_and_offset(amp=wx_amps)
        daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins,counts,rec_readout,rec_avg_all,rec_all = daq_programs.run_daq(3,reps,ro_dur,IQangle=IQ_angle);
           
    
        if i == 0:
    
            rec_g[0] = np.concatenate((rec_g[0], rec_readout[0]))
            rec_g[1] = np.concatenate((rec_g[1], rec_readout[1]))
            rec_read[0] = np.concatenate((rec_read[0], rec_g[0]))
            rec_read[1] = np.concatenate((rec_read[1], rec_g[1]))
    
        elif i == 1:
    
            rec_e[0] = np.concatenate((rec_e[0], rec_readout[0]))
            rec_e[1] = np.concatenate((rec_e[1], rec_readout[1]))
            rec_read[0] = np.concatenate((rec_read[0], rec_e[0]))
            rec_read[1] = np.concatenate((rec_read[1], rec_e[1]))
        else:
    
            rec_f[0] = np.concatenate((rec_f[0], rec_readout[0]))
            rec_f[1] = np.concatenate((rec_f[1], rec_readout[1]))
            rec_read[0] = np.concatenate((rec_read[0], rec_f[0]))
            rec_read[1] = np.concatenate((rec_read[1], rec_f[1]))
    
     
    
    y_train = []
    y_g = IQ_cluster.kMeans(rec_g, n_clusters = 1,plot=True); y_g[:] = +1
    y_e = IQ_cluster.kMeans(rec_e, n_clusters = 1,plot=True); y_e[:] = -1
    y_f = IQ_cluster.kMeans(rec_f, n_clusters = 1,plot=True); y_f[:] = 0
    
    y_train = np.concatenate((y_train,y_g)); y_train = np.concatenate((y_train,y_e)); y_train = np.concatenate((y_train,y_f))
    
    model,x, svm_coef,svm_intercept,Ax,By = IQ_cluster.svm_class(rec_read, y_train, rec_read, 0.05,plot=True)
    
    return model,svm_coef,svm_intercept,Ax,By
#    
#    
    
    
    
    
    
    
    
    
    
    
#####################################################################
#tempQ = np.zeros(101)
#tempI = np.zeros(101)
#avgs = 10
#
#for i in range(avgs): 
#    daq_params, rec_readout_vs_pats, p_vs_pats,I,Q= run_daq(101,10000);
#    tempI = tempI+I; 
#    tempQ=tempQ+Q
#
#
#I_avgs = tempI/avgs
#Q_avgs = tempQ/avgs
#
#plt.plot(I_avgs)
#plt.show()
#plt.plot(Q_avgs)

#tempQ = np.zeros(101)
#tempI = np.zeros(101)
    
#save_dir = r"Z:\candle qubit/"
#np.savetxt(save_dir+'7_4' + 'rabi_ef_without_pi_I_10k_10avgs',I_avgs)
#np.savetxt(save_dir+'7_4' + 'rabi_ef_without_pi_Q_10k_10avgs',Q_avgs)
#np.savetxt(save_dir+'5_28' + 'rabi_ef_with_pi_I_10k_20avgs',I_avgs)
#np.savetxt(save_dir+'5_28' + 'rabi_ef_with_pi_Q_10k_20avgs',Q_avgs)