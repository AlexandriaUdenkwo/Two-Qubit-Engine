# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 17:26:51 2022

@author: crow104
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May 27 15:14:48 2022

@author: crow104
"""

import sys
file_dir= r"C:\Users\crow104\Documents\Python Scripts\data_acquisition_scripts"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Users\crow104\Documents\Python Scripts\sequence_generator"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Data\2020\201021_thermal_readout\Ramsey vs thermal noise"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Users\crow104\Documents\Python Scripts\sequence_generator\py_sequences"
if file_dir not in sys.path: sys.path.append(file_dir)

from generator import *
import os
pi = np.pi
import pyvisa
plt.rcParams.update({'font.size': 13})
from tkinter import Tk
import tkinter as tki
from tkinter.filedialog import askopenfilename, asksaveasfilename, asksaveasfile
import scipy.io as io
import time
from IPython import get_ipython
from scipy import optimize as opt
from matplotlib.patches import Rectangle
import warnings
from scipy.optimize import fsolve


import numpy as np
import matplotlib.pyplot as plt
import FPJPA_sequence as Fs
import wx_programs as wx
import daq_programs_homo
import analysis
import chevron
import bnc
import NH_pi_nopi_RO as thr
from sklearn import svm
from sklearn.cluster import KMeans

instrument_path = r"C:\Users\Crow108\Documents\Python\instr\analyzer"
if instrument_path not in sys.path: sys.path.append(instrument_path )

analyzer_path = r"C:\Users\Crow108\Documents\Python\instr\python_interface\python_without_WX2184C"
if analyzer_path not in sys.path: sys.path.append(analyzer_path )
#save_dir = r"C:\Data\2023\two_qubit_tuna\chevron_sweeps"
save_dir = r"C:/Data/2023/two_qubit_tuna/heterodyne_qubit_freq_sweep/"

#readout bnc (BNC 2)
target_bnc_black_address = 'GPIB0::19::INSTR'

#qubit bnc (BNC 1)
target_bnc_address_6 = 'USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR' #835 with touchscreen 'USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR'


reps=300
ro_dur = 3000
## [RO,PC,Q1,Q2]
wx_amps =[0.7,1.1,.55,.5] # [0.7,1.1,.55,.5]
wx_offs = [0,0.01,0,0] # [0,0.2,0,0]
IQ_angle =250 #225

RO_LO = 6.77#6.77#6.85 #6.6
ROq1 = 6.872626 #6.872870#6.873020 
ROq2 = 6.68862
ssm_geq1 =-0.1985 #-0.0767 #-0.0767#-0.0564  #   #-0.0796#Hypres    #-0.0796#Keithley                  #Keithley#-0.0564#-.2 #-0.077737#-.090+0.00025#-.126#-0.3275  
ssm_geq2 = -0.1468#64#-0.1013-0.0012+0.0001 #-0.088
ssm_efq1 = -.2375
ssm_efq2 = -.249
#q1 = qubit(1)
#q1.ssbge
#q1.ROfre

which_qubit = 1 #qubit 1 = 1; qubit 2 = 0; both qubits = -1

if which_qubit ==1:
    pi_ge_time =  34 #Keithley  #41#Hypres        #Keithley#33#41 #40.6 #37 #38 #q1
if which_qubit == 0:
    pi_ge_time = 44#50#49#50 #q2
    
qubit_1_thr= [-1500,-400]#[-800,550]
qubit_2_thr=[-2500,-1100]#[500,1700]#
IQ_angle_q1 =  110 #130
IQ_angle_q2 = 30
#pi_ef_time = 40.46

# set 1 for which sequence to run
spec_ge= 1
spec_ef= 0
cluster_thr = 0
pis_nopi = 0
pi_2_nopi_2 = 0
no_pi_pm = 0
n_pi_2 = 0
run_rabi = 0
run_rabi_ge_2qubit = 0  
run_rabi_wx_sweep = 0
run_rabief = 0
run_rabi_J = 0
run_rabi_2qubit = 0
run_T1 = 0
run_T1ef = 0
run_T1_modified=0
run_ramsey = 0
run_ramsey_ef = 0
run_sweep = 0
run_J_sweep = 0
run_bnc_sweep = 0
run_bnc_sweep_2q = 0
run_bnc_power_sweep = 0
run_bnc_pow_freq_sweep = 0
run_ssbef_sweep = 0 #chevron for ssb_ef
run_x_measure = 0
run_y_measure = 0
teaching=0
ef_teach=0
parametric_coupling=0
parametric_coupling_time_domain=0
parametric_coupling_sweep = 0

if which_qubit == 1:
    scale_matrix = np.array([[0.48206667, 0.48806667, 0.48213333],
                             [0.1,         0.,         0.        ], 
                             [0.51793333, 0.51193333, 0.51786667]])
    
    scale_matrix = np.linalg.inv(scale_matrix)
    ssm_ge = ssm_geq1
    ssm_ef = -0.262
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
   # ssm_ef = ssm_ge - 0.15
    #RO_freq = ROq1
    RO_freq = RO_LO
    qubit_bnc =  4.3#4.312013#4.312113 #4.312013 #4.312024
    bnc.set_bnc_output(RO_freq, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc, power_dBm=7,bnc_addr=target_bnc_address_6)
if which_qubit == 0:
    scale_matrix = np.array([[9.66666667e-01, 3.73666667e-01, 1.31400000e-01],
                             [3.25333333e-02, 5.64333333e-01, 2.29000000e-01],
                             [8.00000000e-04, 6.20000000e-02, 6.39600000e-01]])
    scale_matrix = np.linalg.inv(scale_matrix)
    scale_matrix_ES = np.array([[0.99593333, 0.58973333, 0.18366667],
                                [0.00273333, 0.31706667, 0.16833333],
                                [0.00133333, 0.0932    , 0.648     ]])
    scale_matrix_ES = np.linalg.inv(scale_matrix_ES)
    ssm_ge = ssm_geq2
    ssm_ef = -0.25721#-0.25675#-0.2804#ssm_ge - 0.15# #-0.4091#-0.409 #ssm_ge - 0.15
    ROIF2 = ROq2 - RO_LO
    RO_freq = RO_LO
    #RO_freq = ROq2
    qubit_bnc =4.3#4.256
    bnc.set_bnc_output(RO_freq, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc, bnc_addr=target_bnc_address_6)
    
if which_qubit == -1:
    #Do this for heterodyne
    qubit_bnc = 4.3
#    ssm_geq1 = qubit_bnc - (4.312013 +-.090+0.00025)
#    ssm_geq2 = qubit_bnc - (4.256 +-0.1013-0.0012+0.0001)
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
    RO_freq = RO_LO
    bnc.set_bnc_output(RO_freq, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc, power_dBm=7,bnc_addr=target_bnc_address_6)


### TODO how should we save model? pickle?
if cluster_thr: #run to get "model" for blob clustering
    model,svm_coef,svm_intercept,Ax,By = thr.cluster_model(ssm_ge,ssm_ef,pi_ge_time,pi_ef_time,ro_dur,IQ_angle)
    Ax = []
    By = []
    num_lines = np.shape(svm_coef)[0]
    for i in range (num_lines):
        Ax.append(svm_coef[i][0]/svm_intercept[i])
        By.append(svm_coef[i][1]/svm_intercept[i])


if pis_nopi: #run to get scale_matrix = p 3x3 for g-e-f correction
    num_steps = 3
    reps = 15000
    Fs.pi_nopi(off = 0,coef=0,coefpief=0,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
    Fs.pi_nopi(off = 1,coef=1,coefpief=0,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
    Fs.pi_nopi(off = 2,coef=1,coefpief=1,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge = ssm_ge,ssm_ef=ssm_ef);

if pi_2_nopi_2: #run to correct two qubit rabi tomography
    num_steps = 3
    reps = 15000
    Fs.pi_nopi(off = 0,coef=0,coefpief=0,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge = ssm_ge,ssm_ef=ssm_ef); #|g>, keeps population in ground state
    Fs.pi_nopi(off = 1,coef=1,coefpief=0,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge = ssm_ge,ssm_ef=ssm_ef); #|e>, sends population to |e>
    #Fs.pi_nopi(off = 2,coef=1,coefpief=1,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge = ssm_ge,ssm_ef=ssm_ef); #|f>, send population to |f>
    
    
if n_pi_2: #do multiple pi/2s, use to tune pi/2 time
    num_steps = 3
    reps = 15000
    Fs.n_pi_2(off = 0,coef=0,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge = ssm_ge,q=which_qubit); #|g>
    Fs.n_pi_2(off = 1,coef=4,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge = ssm_ge,q=which_qubit);
    Fs.n_pi_2(off = 2,coef=8,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge = ssm_ge,q=which_qubit); #multiples of 4 to send to back to |g>



if no_pi_pm:
    amp_plus = 1-0.03
    amp_minus = 1+0.05
    phase_plus = 0
    phase_minus = 180
    num_steps = 3
    reps = 15000
    Fs.nopi_plus_minus(coef=0,coefES=0,coef_rev=1,off=0,amp_ef=0,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge=ssm_ge,ssm_ef=ssm_ef,es_phase=0) #|g> prepare
    Fs.nopi_plus_minus(coef=1,coefES=1,coef_rev=-1,off=1,amp_ef=amp_plus,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge=ssm_ge,ssm_ef=ssm_ef,es_phase=phase_plus) #|+> prepare and rotate to |e>
    Fs.nopi_plus_minus(coef=1,coefES=1,coef_rev=1,off=2,amp_ef=amp_minus,pi_ge=pi_ge_time,pi_ef=pi_ef_time,ssm_ge=ssm_ge,ssm_ef=ssm_ef,es_phase=phase_minus) #|-> prepare and rotate to |f>
    
if spec_ge:
    num_steps = 101
    f1=-.3 #-0.17
    f2= -.01#-0.13
    if which_qubit == 0:
        Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF=ROIF2,q=which_qubit)
    if which_qubit == 1:
        Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF=ROIF1,q=which_qubit)
        
#        num_steps=101,ssm_start=-.05,ssm_stop=-.15,spec_amp=.5,ROIF2 = 0,q=0
#        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
#    Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp= 1, q = which_qubit)
    

if spec_ef:
    num_steps = 51
    f1=-.2
    f2=-.25
    if which_qubit == 0:
#        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
        Fs.spectroscopy_ef(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF=ROIF2,q=which_qubit)
    if which_qubit == 1:
        Fs.spectroscopy_ef(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF=ROIF1,q=which_qubit)
        
#    Fs.spectroscopy_ef(num_steps,ssm_ge,ssm_start=f1,ssm_stop=f2,spec_amp=0.1)

if run_rabi:
    num_steps = 81
    sweep_time =300 #ns
    if which_qubit == 0:
        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
    if which_qubit == 1:
        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF1,which_qubit)
        
if run_rabi_ge_2qubit:
    num_steps = 81
    sweep_time =200 #ns
    Fs.rabi_ge_2qubit(num_steps,sweep_time,ssm_geq1,ssm_geq2, ROIF1,ROIF2)
    
if run_rabi_wx_sweep:
    num_steps = 81
    sweep_time =200 #ns
    Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)

if run_rabief:
    num_steps = 81
    sweep_time =200 #ns
    if which_qubit == 0:
#        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
        Fs.rabi_ef(num_steps,sweep_time,pi_ge_time,ssm_ge =ssm_geq2,ssm_ef=ssm_efq2,ef_amp=1,ROIF=ROIF2,q=which_qubit)
    if which_qubit == 1:
#        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF1,which_qubit)
        Fs.rabi_ef(num_steps,sweep_time,pi_ge_time,ssm_ge =ssm_geq1,ssm_ef=ssm_efq1,ef_amp=1,ROIF=ROIF1,q=which_qubit)
#    Fs.rabi_ef(num_steps,sweep_time,pi_ge_time,ssm_ge,ssm_ef,q=which_qubit)
    
if run_rabi_J: # 
    ef_amp = 0.01 #EP occurs at amp = 0.00227 
    num_steps = 81
    sweep_time =10000 #ns
    Fs.rabi_J(num_steps,sweep_time,ef_amp,pi_ge_time,pi_ef_time,ssm_ge,ssm_ef)
    
if run_rabi_2qubit:
    num_steps = 51
    sweep_time = 400
    Fs.rabi_ge_2qubit(num_steps,sweep_time,ssm_ge,ssm_geq1,ssm_geq2)
    
if run_T1:
    num_steps = 81
    sweep_time =100000#ns
    if which_qubit == 0:
        Fs.T1_ge(num_steps,sweep_time,ssm_ge,pi_ge_time,q=which_qubit,ifload = 1,ROIF=ROIF2)
    if which_qubit == 1:
        Fs.T1_ge(num_steps,sweep_time,ssm_ge,pi_ge_time,q=which_qubit,ifload = 1,ROIF=ROIF1)
    
if run_T1ef:
    num_steps = 81
    sweep_time =50000 #ns
    Fs.T1_ef(num_steps,sweep_time,ssm_ge,ssm_ef,pi_ge_time,pi_ef_time)
    
if run_T1_modified:
    num_steps = 51
    sweep_time =40000 #ns
    Fs.T1_ge_modified(num_steps,sweep_time)

if run_ramsey:
    piecho = 0
    osc_num = 4
    num_steps = 81
    sweep_time =20000 #ns
    if which_qubit == 0:
        Fs.ramsey(num_steps,sweep_time,piecho,osc_num,ssm_ge,pi_ge_time,ROIF2,q=which_qubit)
        #Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
    if which_qubit == 1:
        #Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF1,which_qubit)
        Fs.ramsey(num_steps,sweep_time,piecho,osc_num,ssm_ge,pi_ge_time,ROIF1,q=which_qubit)
    
    
if run_ramsey_ef:
    piecho = 0
    osc_num = 0
    num_steps = 81
    sweep_time =10000 #ns
    es_amp = 1
    Fs.ramsey_ef(num_steps,sweep_time,piecho,osc_num,ssm_ge,ssm_ef,pi_ge_time,pi_ef_time,es_amp,q=which_qubit)
    
if run_x_measure:
    num_steps = 81
    sweep_time =200 #ns
    Fs.x_measurement(num_steps,sweep_time,ssm_ge,which_qubit, pi_ge_time)
    
if run_y_measure:
    num_steps = 81
    sweep_time =200 #ns
    Fs.y_measurement(num_steps,sweep_time,ssm_ge,which_qubit, pi_ge_time)

#if teaching:
#    num_steps = 101
#    sweep_time=20000
#    Fs.teaching_test(num_steps=num_steps,sweep_time=sweep_time,pi_ge_time=pi_ge_time,ssm_ge=ssm_ge,ROIF2 = ROIF1,q=which_qubit)  

if ef_teach:
    num_steps = 81
    sweep_time =200 #ns
    Fs.ef_test(num_steps,sweep_time,pi_ge_time,ssm_ge,ssm_ef,q=which_qubit)
    
if parametric_coupling:
    num_steps = 101
    f1=-.025
    f2=-.04
    if which_qubit == 0:
#        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
        Fs.parametric_coupling(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_start=f1,ssm_stop=f2,spec_amp=1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
        #parametric_coupling(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_start=-.15,ssm_stop=-.25,spec_amp=.5,ROIF1=0,ROIF2=0,q=0)
#    if which_qubit == 1:
#        Fs.parametric_coupling(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF=ROIF1,q=which_qubit)
        
if parametric_coupling_time_domain:
    num_steps = 101
    f_parametric=-0.0344#-0.0337#-.0342
    sweep_time = 5000
    if which_qubit == 0:
#        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
        Fs.parametric_coupling_time_domain(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_para=f_parametric,spec_amp=1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit,sweep_time=sweep_time)
  


### no thresholding ###
#daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_pi,counts_pi,rec_readout,rec_avg_all,rec_all = daq_programs_homo.run_daq(num_steps,reps,ro_dur,IQangle=IQ_angle)
#plt.plot(Q);plt.show();plt.plot(I);plt.show()

#if run_rabi_wx_sweep:
#    #wx_amps = [2,.5,.55,.5]
#    wx_ch1amp_steps = np.arange(0.5, 2.0, .1)
#    for i in range(len(wx_ch1amp_steps)):
#        wx_amps[0] = wx_ch1amp_steps[i]
#        wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
#        
#        #rec_avg_all, rec_readout, rec_avg_vs_pats, rec_all_het, bins, counts = daq_programs_homo.run_daq_het(ssm_if=ROIF2, num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur, verbose=True)
#        I = rec_avg_vs_pats[0];plt.plot(I);plt.show()
#        Q = rec_avg_vs_pats[1];plt.plot(Q);plt.show()
#        
#        np.savetxt(save_dir+'\q2_amp'+str(wx_amps[0])+'_I_rabi_ROIF2_'+str(ROIF2)+'RO_LO_'+str(RO_LO),I)
#        np.savetxt(save_dir+'\q2_amp'+str(wx_amps[0])+'_Q_rabi_ROIF2_'+str(ROIF2)+'R0_LO_'+str(RO_LO),Q)
#        
    
#else:

##Heterodyne code
#if which_qubit == 0:
#    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
#
#    rec_avg_all, rec_readout, rec_avg_vs_pats, rec_all_het, bins, counts = daq_programs_homo.run_daq_het(ssm_if=ROIF2, num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur, verbose=True)
#    I = rec_avg_vs_pats[0];plt.plot(I);plt.show()
#    Q = rec_avg_vs_pats[1];plt.plot(Q);plt.show()
#
#if which_qubit == 1:
#    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
#
#    rec_avg_all, rec_readout, rec_avg_vs_pats, rec_all_het, bins, counts = daq_programs_homo.run_daq_het(ssm_if=ROIF1, num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur, verbose=True)
#    I = rec_avg_vs_pats[0];plt.plot(I);plt.show()
#    Q = rec_avg_vs_pats[1];plt.plot(Q);plt.show()

if parametric_coupling_sweep:
    
#    wx_amps =[0.7,.8,.55,.5]
#    wx_offs = [0,0,0,0]
    sweep_steps_amp = 21
    sweep_steps_off = 21
    wx_ch2amp_steps = np.linspace(0.5,1.1,sweep_steps_amp)#np.arange(0.5, 1.5, .1)
    wx_ch2offs_steps = np.linspace(0,0.2,sweep_steps_off)#np.arange(0, 1.0, .1)
    guess_vals = [0.6,0.05,0.15,90,0.5]
    
    coupling_freq = np.zeros((sweep_steps_off,sweep_steps_amp))
    parametric_freq = np.zeros(sweep_steps_off)
    #out_I, out_Q = np.zeros( (2, num_points, steps_in_seq))np.zeros()
    i=0
    j=0
    for i in range(len(wx_ch2offs_steps)):
        print('current amp:',wx_ch2amp_steps[j])
        print('current offset:',wx_ch2offs_steps[i])
        print("offset i:",i)
        print("amp j:",j)
        f1=-.025
        f2=-.04
        num_steps = 101
        
        Fs.parametric_coupling(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_start=f1,ssm_stop=f2,spec_amp=1,ROIF1=ROIF1,ROIF2=ROIF2,q=0)
       
        wx_amps[1] = wx_ch2amp_steps[j]
        wx_offs[1] = wx_ch2offs_steps[i]
        wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
#            
        rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
                                                                                                                                                                                                
        print("Qubit 1:")
        P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
        I_Q1 = rec_avg_vs_pats_1[0];
        Q_Q1 = rec_avg_vs_pats_1[1];
        plt.plot(I_Q1);plt.title('I Q1');plt.show()
        plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
        
        print("Qubit 2:")
        P_Q2 = prob_vs_pats_2[0];plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
        I_Q2 = rec_avg_vs_pats_2[0];
        Q_Q2 = rec_avg_vs_pats_2[1];
        plt.plot(I_Q2);plt.title('I Q2');plt.show()
        plt.plot(Q_Q2);plt.title('Q Q2');plt.show()
#            
        #Find frequency of parametric coupling
        freq = np.linspace(f1,f2,num_steps)

        Qrange = abs(np.max(Q_Q1)-np.min(Q_Q1))
        Irange = abs(np.max(I_Q1)-np.min(I_Q1))
        if Qrange>Irange:
            freq_index = np.where(Q_Q1 == np.amin(Q_Q1))  
            print("Q_Q1")
            plt.plot(freq,Q_Q1)
        if Irange>Qrange:
            freq_index = np.where(I_Q1 == np.amin(I_Q1))     
            print("I_Q1")
            plt.plot(freq,I_Q1)
        ssm_ge_para = freq[freq_index]
        print(ssm_ge_para)
#                
        parametric_freq[i] = ssm_ge_para
        f_parametric=ssm_ge_para#-0.0343
        for j in range(len(wx_ch2amp_steps)):
            print('current amp:',wx_ch2amp_steps[j])
            print('current offset:',wx_ch2offs_steps[i])
            print("offset i:",i)
            print("amp j:",j)
            
            sweep_time = 5000
            Fs.parametric_coupling_time_domain(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_para=f_parametric,spec_amp=1,ROIF1=ROIF1,ROIF2=ROIF2,q=0,sweep_time=sweep_time)
            
            wx_amps[1] = wx_ch2amp_steps[j]
            wx_offs[1] = wx_ch2offs_steps[i]
            wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
        
            rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
                                                                                                                                                                                                    
            print("Qubit 1:")
            P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
            I_Q1 = rec_avg_vs_pats_1[0];
            Q_Q1 = rec_avg_vs_pats_1[1];
            plt.plot(I_Q1);plt.title('I Q1');plt.show()
            plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
            
            print("Qubit 2:")
            P_Q2 = prob_vs_pats_2[0];plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
            I_Q2 = rec_avg_vs_pats_2[0];
            Q_Q2 = rec_avg_vs_pats_2[1];
            plt.plot(I_Q2);plt.title('I Q2');plt.show()
            plt.plot(Q_Q2);plt.title('Q Q2');plt.show()
            
#            guess_vals = fit_vals
          
            times =np.linspace(0,sweep_time/1000,num_steps);
            plt.plot(times,P_Q1,label='Q1');
            plt.plot(times,P_Q2,label='Q2');
            plt.legend();plt.xlabel('Time (\u03BCs)');
            plt.show();
            fit_vals,_,_,_ = analysis.fit_sine_decay(times,P_Q1,guess_vals);
#            guess_vals = fit_vals
            coupling_freq[i][j] = fit_vals[0]
            swap_time = abs((1/2/fit_vals[0])*1000)
            print("half swap time = {} ns".format(swap_time/2))
    
    plt.imshow(coupling_freq)
    np.savetxt(save_dir+'parametric_drive_sweep_0.5to1.1amp_0to0.2offset_coupling_freq',coupling_freq)
    np.savetxt(save_dir+'parametric_drive_sweep_0.5to1.1amp_0to0.2offset_parametric_freq',parametric_freq)
    
else:
    ##Do actual heterodyne
    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
    rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
                                                                                                                                                                                                    
    print("Qubit 1:")
    P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
    I_Q1 = rec_avg_vs_pats_1[0];
    Q_Q1 = rec_avg_vs_pats_1[1];
    plt.plot(I_Q1);plt.title('I Q1');plt.show()
    plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
    
    print("Qubit 2:")
    P_Q2 = prob_vs_pats_2[0];plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
    I_Q2 = rec_avg_vs_pats_2[0];
    Q_Q2 = rec_avg_vs_pats_2[1];
    plt.plot(I_Q2);plt.title('I Q2');plt.show()
    plt.plot(Q_Q2);plt.title('Q Q2');plt.show()   

### thresholding with I axis ###
#daq_params, t_histo, p,a,b,n_readout = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=m,num_patterns=num_steps, 
#                                                                                       num_records_per_pattern=reps,authr=0,fg=2,ro_dur=ro_dur,IQangle=IQ_angle)
#m = daq_params.threshold
#p_readout = p
#y=p_readout[1]

### thresholding with cluster blobs ###
#daq_params, t_histo, p,a,b,n_readout = daq_programs.run_daq_cluster_threshold(model,Ax=Ax,By=By,num_patterns=num_steps,
#                                                                                 num_records_per_pattern=reps,authr=0,fg=3,ro_dur=ro_dur,IQangle=IQ_angle)
#p_readout = p
#
#p_readout= np.matmul(scale_matrix,p)
#plt.plot(p_readout[0],label='|g>');plt.plot(p_readout[1],label='|e>');plt.plot(p_readout[2],label='|f>')
#plt.legend();plt.title('data scaled with matrix');plt.show()
####
#################### regular readout
#y=p_readout[2]
####
####
########## with post selection
#p_post = analysis.p_readout_postselected(p_readout)
#y = p_post[2]

## analysis code for each sequence ###
if which_qubit == 0:
    I = I_Q2; Q = Q_Q2
    y = P_Q2
if which_qubit == 1:
    I = I_Q1; Q = Q_Q1
    y = P_Q1
    
if parametric_coupling:
    freq = np.linspace(f1,f2,num_steps)
    if which_qubit == 0:
        Qrange = abs(np.max(Q_Q1)-np.min(Q_Q1))
        Irange = abs(np.max(I_Q1)-np.min(I_Q1))
        if Qrange>Irange:
            freq_index = np.where(Q_Q1 == np.amax(Q_Q1))  
            print("Q_Q1")
            plt.plot(freq,Q_Q1)
        if Irange>Qrange:
            freq_index = np.where(I_Q1 == np.amin(I_Q1))     
            print("I_Q1")
            plt.plot(freq,I_Q1)
        ssm_ge_para = freq[freq_index]
        print(ssm_ge_para)

if parametric_coupling_time_domain:
    times =np.linspace(0,sweep_time/1000,num_steps);
    plt.plot(times,P_Q1,label='Q1');
    plt.plot(times,P_Q2,label='Q2');
    plt.legend();plt.xlabel('Time (\u03BCs)');
    plt.show();
    fit_vals,_,_,_ = analysis.fit_sine_decay(times,P_Q1,[0.6,0.05,0.15,90,0.5]);
    swap_time = abs((1/2/fit_vals[0])*1000)
    print("half swap time = {} ns".format(swap_time/2))
    
if teaching:
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,I,guess_vals=[0.04,0.05,35,-90,490])
    pi_ge = abs((1/2/pi_ge_fit_vals[0])*1000)


if ef_teach:
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ef_fit_vals,_,_,_ = analysis.fit_sine_decay(times,Q,guess_vals=[12,0.5,0.14,120,-180])
    pi_ef = abs((1/2/pi_ge_fit_vals[0])*1000)
    print("\u03C0_ef time = {} ns".format(pi_ef))


if run_rabi:
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,y,guess_vals=[6,0.05,0.07,90,0.5])
    pi_ge = abs((1/2/pi_ge_fit_vals[0])*1000)
    print("\u03C0_ge time = {} ns".format(pi_ge))
#    
if run_rabief:
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ef_fit_vals,_,_,_ = analysis.fit_sine_decay(times,Q,guess_vals=[20,0.5,0.14,-600,-100])#[10,2,2,1,-179]
    pi_ef = abs((1/2/pi_ef_fit_vals[0])*1000)
    print("\u03C0_ef time = {} ns".format(pi_ef))
    
if run_x_measure:
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,I,guess_vals=[6,0.05,0.07,90,-46])
    pi_ge = abs((1/2/pi_ge_fit_vals[0])*1000)
#    np.savetxt(save_dir+'\q1_amp'+str(wx_amps[2])+'_I_xrabi_flux_1mA_'+str(qubit_bnc)+'_ghz',I)
#    np.savetxt(save_dir+'\q2_amp'+str(wx_amps[2])+'_Q_xrabi_flux_1mA_'+str(qubit_bnc)+'_ghz',Q)
    print("\u03C0_ge time = {} ns".format(pi_ge)) 
    
if run_y_measure:
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,I,guess_vals=[6,0.05,0.07,90,-46])
    pi_ge = abs((1/2/pi_ge_fit_vals[0])*1000)
#    np.savetxt(save_dir+'\q1_amp'+str(wx_amps[2])+'_I_yrabi_flux_1mA_'+str(qubit_bnc)+'_ghz',I)
#    np.savetxt(save_dir+'\q2_amp'+str(wx_amps[2])+'_Q_yrabi_flux_1mA_'+str(qubit_bnc)+'_ghz',Q)
    print("\u03C0_ge time = {} ns".format(pi_ge)) 

#    
if run_rabi_J:
    times = np.linspace(0,sweep_time/1000,num_steps)
    analysis.fit_sine_decay(times,y,guess_vals=[4,0.9,-0.4,-600,-100])#[10,2,2,1,-179]
    
if spec_ge:
    freq = np.linspace(f1,f2,num_steps)
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        freq_index = np.where(Q == np.amax(Q))
        print("Q")
        plt.plot(freq,Q)
    if Irange>Qrange:
        freq_index = np.where(I == np.amax(I))
        print("I")
        plt.plot(freq,I)
#    freq_index = np.where(y == np.max(y))
    ssm_ge = freq[freq_index]
    print(ssm_ge)
#    plt.plot(qubit_bnc+freq,Q)
    
#    analysis.fit_lorentzian(freq,y,guess_vals=[2.6457e-10,-2.949e-10,0.5,-0.0423])#[10,2,2,1,-179]

if spec_ef:
    freq = np.linspace(f1,f2,num_steps)
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        freq_index = np.where(Q == np.amax(Q))
    if Irange>Qrange:
        freq_index = np.where(I == np.amin(I))
#    freq_index = np.where(y == np.max(y))
    ssm_ef = freq[freq_index]
    print(ssm_ef)
         

if run_T1:
    times = np.linspace(0,sweep_time/1000,num_steps)
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        T1_ge_fit_vals,_,_,_ = analysis.fit_exp_decay(times,Q,guess_vals=[-400,0.001,-1700])
    if Irange>Qrange:
        T1_ge_fit_vals,_,_,_ = analysis.fit_exp_decay(times,I,guess_vals=[400,1,-1700])
#    T1_ge_fit_vals,_,_,_ = analysis.fit_exp_decay(times,P_Q2,guess_vals=[10,1,0.5])
    T1_ge = 1/T1_ge_fit_vals[1]
    print("T1_ge = {} \u03BCs".format(T1_ge))
    
if run_T1ef:
    times = np.linspace(0,sweep_time/1000,num_steps)
    analysis.fit_exp_decay(times,y,guess_vals=[-4,0.08,0.5])
    
if run_T1_modified:
    times = np.linspace(0,sweep_time/1000,num_steps)
    analysis.fit_exp_decay(times,I,guess_vals=[-4,0.08,0.5])

if run_ramsey or run_ramsey_ef:
    times = np.linspace(0,sweep_time/1000,num_steps)
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        T2_fit_vals,_,_,_ = analysis.fit_sine_decay(times,Q,guess_vals=[0.8,0.4,0.08,-260,-1600])
    if Irange>Qrange:
        T2_fit_vals,_,_,_ = analysis.fit_sine_decay(times,I,guess_vals=[0.8,0.4,0.08,-260,-1600])
    T2_fit_vals,_,_,_ = analysis.fit_sine_decay(times,y,guess_vals=[0.4,1,180,90,0.5])
    T2 = 1/T2_fit_vals[1]
    print("T2* = {} \u03BCs".format(T2))
#    analysis.fit_exp_decay(times,Q,guess_vals=[5,1,120])
    
#### run keithley sweep ####
if run_sweep:
#    center_cur = -0.1
    out_Q, out_I = chevron.sweep_keithley(-0.8,-0.6,21,num_steps,reps,ro_dur)
    np.savetxt(save_dir+'specq2_amp.1_-0.4,0ssb_I_keithleysweep_-0.8,-0.6mA,21steps',out_I)
    np.savetxt(save_dir+'specq2_amp.1_-0.4,0ssb_Q_keithleysweep_-0.8,-0.6mA,21steps',out_Q)

if run_J_sweep:
    num_steps = 81
    sweep_time = 1000
    out_y = chevron.sweep_rabi_J(ssm_ef-0.005,ssm_ef+0.005,51,model,scale_matrix,Ax,By,num_steps,reps,ro_dur,sweep_time,pi_ge_time,pi_ef_time,ssm_ge,ssm_ef,IQ_angle)
#    np.savetxt(save_dir+'sweep_JEP_-0.02to0.05_51steps_rabi_ef_4us_51steps_10kavg_scaledmatrix',out_y)
    np.savetxt(save_dir+'sweep_ssbef_chevron_10MHz_rabi_ef_0.1amp_1us_81steps_2kavg_scaledmatrix',out_y)

if run_ssbef_sweep: #chevron
    sweep_steps = 41
    num_steps = 51
    sweep_time = 1000
    start_freq = ssm_ef-0.002
    stop_freq = ssm_ef+0.002
    out_y = chevron.sweep_rabi_ef(start_freq, stop_freq,sweep_steps,model,scale_matrix,Ax,By,num_steps,reps,ro_dur,sweep_time,pi_ge_time,ssm_ge,ssm_ef,IQ_angle,which_qubit)
    np.savetxt(save_dir+'sweep_ssbef_chevron_4MHz_rabi_ef_0.08amp_1us_51steps_2kavg_scaledmatrix',out_y)

if run_bnc_sweep:
    sweep_steps = 101
    start_freq = qubit_bnc - 0.1
    stop_freq = qubit_bnc + 0.1
    out_Q, out_I = chevron.sweep_bnc_freq(start_freq, stop_freq,sweep_steps,num_steps,reps,ro_dur, IQ_angle)
    np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q1_I',out_I)
    np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q1_Q',out_Q)
    
if run_bnc_sweep_2q:
    sweep_steps = 401
    start_freq = qubit_bnc - 0.2
    stop_freq = qubit_bnc + 0.2
    out_Q1, out_I1,out_Q2, out_I2 = chevron.sweep_bnc_freq_2q(start_freq, stop_freq,sweep_steps,num_steps,reps,ro_dur, ROIF1,ROIF2, deg_1 = 0, deg_2 = 0,verbose=True)
    np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q1_I',out_I1)
    np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q1_Q',out_Q1)
    np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q2_I',out_I2)
    np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q2_Q',out_Q2)
    
#ROIF1,ROIF2, deg_1 = 0, deg_2 = 0,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur, verbose=True    

if run_bnc_power_sweep:    
    out_Q, out_I = chevron.sweep_bnc_power(-30, -10, 21 ,num_steps,reps,ro_dur)
    np.savetxt(save_dir+'-155.78mVsweep_BNCpower_to_filter4.2GHz_-30_-10dBm_10dB_Attn_ram15usoutI',out_I)
    np.savetxt(save_dir+'-155.78mVsweep_BNCpower_to_filter4.2GHz_-30_-10dBm_10dB_Attn_ram15usoutQ',out_Q)
#    out_y = chevron.sweep_bnc_power(0,19,21,num_steps,reps,ro_dur)
#    np.savetxt(save_dir+'-2640mVsweep_BNCpower_to_filter_0to19dBm_40dB_Attn_outy_10usram4osc_ssbsweep',out_y)
    
if run_bnc_pow_freq_sweep:
    #12.67
    out_bins_pi, out_counts_pi,out_bins_nopi, out_counts_nopi = chevron.sweep_bnc_freq_and_power(12.67-.5, 12.67+.5, 2,12.7239-0.005, 12.7239+0.005, 1,10000,2800)#pi,pf,psteps,fi,ff,fsteps,reps,ro_dur
    np.savetxt(save_dir+'sweep_JPA_BNCpow_freq_bins_pi_q1',out_bins_pi)
    np.savetxt(save_dir+'sweep_JPA_BNCpow_freq_counts_pi_q1',out_counts_pi)
    np.savetxt(save_dir+'sweep_JPA_BNCpow_freq_bins_nopi_q1',out_bins_nopi)
    np.savetxt(save_dir+'sweep_JPA_BNCpow_freq_counts_nopi_q1',out_counts_nopi)