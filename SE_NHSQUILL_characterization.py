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
save_dir = r"C:\Data\2023\entanglement_engine_data"

#readout bnc (BNC 2)
target_bnc_black_address = 'GPIB0::19::INSTR'

#qubit bnc (BNC 1)
target_bnc_address_6 = 'USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR' #835 with touchscreen 'USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR'


reps=1000
ro_dur = 5000
## [RO,QC,Q1,Q2]
coupler_off_value = .8
wx_amps =[.6,1.6,.5,.5]
wx_offs = [0,coupler_off_value,0,0]
       

import keithley2401
keithley2401.set_current(19.6, step_size_mA=0.1) #20 mhz detuning 6.5, 6 mhz detuning 6.85, 0 mhz detuning 6.95


qubit_bnc =  4.32
RO_LO = 6.77
ROq1 =  6.872982 #6.87270 
ROq2 =  6.68923 #6.68992
ssm_geq1 = -.153 
ssm_geq2 = -.1538
ssm_efq1 = -.2375
ssm_efq2 = -.249
pi_ge_time_q1 = 33
pi_ge_time_q2 = 45


which_qubit = 1 #qubit 1 = 1; qubit 2 = 0; both qubits = -1

if which_qubit ==1:
    pi_ge_time = pi_ge_time_q1
if which_qubit == 0:
    pi_ge_time = pi_ge_time_q2
    
qubit_1_thr= [-1000,570]
qubit_2_thr=[-2050,-95]
IQ_angle_q1 = 10
IQ_angle_q2 = 70#190 #increasing rotates CCW
#pi_ef_time = 40.46

#Scaling
scale_matrix1 = np.array([[0.95313333, 0.15633333, 0.0],
                          [0.04686667, 0.84366667, 0.0],
                          [0.        , 0.        , 1.0]])

scale_matrix1 = np.linalg.inv(scale_matrix1)


scale_matrix2 =np.array([[0.97113333, 0.15706667, 0.0],
                         [0.02886667, 0.84293333, 0.0],
                         [0.        , 0.        , 1.0]])
scale_matrix2 = np.linalg.inv(scale_matrix2)




# set 1 for which sequence to run
spec_ge=0
spec_ge_coupler_switch =0
spec_ef= 0
cluster_thr = 0
pis_nopi = 0
threshold_sweep = 0
pi_2_nopi_2 = 0
no_pi_pm = 0
n_pi_2 = 0
piA_piB_roB = 0
run_rabi =0
run_vacuum_rabi = 0
run_rabi_ge_2qubit = 0
run_rabi_wx_sweep = 0
run_rabief = 0
run_rabi_J = 0
run_rabi_2qubit = 0
run_T1 = 0
run_T1_ge_2q_RO = 0
run_T1_M_ge_2q_RO =0
run_T1ef = 0
run_T1_modified=0
run_ramsey = 0
run_ramsey_ef = 0
run_sweep =1
run_wx_sweep = 0
run_coupler_switch_sweep =0
run_wx_coupler_switch_sweep = 0
run_J_sweep = 0
run_bnc_sweep =0
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
sweep_ro_integration = 0



if which_qubit == 1:
    ssm_ge = ssm_geq1
    ssm_ef = -0.262
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
   # ssm_ef = ssm_ge - 0.15
    #RO_freq = ROq1
    RO_freq = RO_LO
    bnc.set_bnc_output(RO_freq, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc,power_dBm=13,bnc_addr=target_bnc_address_6)
if which_qubit == 0:
    ssm_ge = ssm_geq2
    ssm_ef = -0.25721
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
    RO_freq = RO_LO
    #RO_freq = ROq2
   # qubit_bnc =4.3#4.256
    bnc.set_bnc_output(RO_freq, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc,power_dBm=13, bnc_addr=target_bnc_address_6)
    
if which_qubit == -1:
    #Do this for heterodyne
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
    RO_freq = RO_LO
    bnc.set_bnc_output(RO_freq, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc, power_dBm=13,bnc_addr=target_bnc_address_6)

if pis_nopi: #run to correct two qubit rabi tomography
    num_steps = 3
    reps = 15000
    if which_qubit == 0:
        
        Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|g>, keeps population in ground state
        Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|e>, sends population to |e>
        #Fs.pi_nopi(off = 2,coef=1,coefpief=1,pi_ge=pi_ge_time,pi_ef=pi_ef_time,q= which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
        
    if which_qubit == 1:
        Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge); #|g>, keeps population in ground state
        Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,r=which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge); #|e>, sends population to |e>
        #Fs.pi_nopi(off = 2,coef=1,coefpief=1,pi_ge=pi_ge_time,pi_ef=pi_ef_time,q= which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge,ssm_ef=ssm_ef);

if piA_piB_roB: 
    num_steps = 3
    reps = 15000
    #pi pulse qubit 0 (A) readout qubit 1B
    Fs.pi_nopi(off = 0,coef=1,pi_ge=pi_ge_time_q2,q = 0,r=1,ROIF=ROIF1,ssm_ge = ssm_geq2); #|g>, keeps population in ground state
    #pi pulse qubit 1 (B) readout qubit 1B
    Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time_q1,q = 1,r=1,ROIF=ROIF1,ssm_ge = ssm_geq1); #|e>, sends population to |e>
        #Fs.pi_nopi(off = 2,coef=1,coefpief=1,pi_ge=pi_ge_time,pi_ef=pi_ef_time,q= which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
        

    
if spec_ge:
    num_steps = 101
    f1=-0.01
    f2= -0.31
    if which_qubit == 0:  #readout on both qubits, excited q2
        Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
    if which_qubit == 1: ##readout on both qubits, excited q1
        Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
        
if spec_ge_coupler_switch:
    num_steps = 101
    f1=-.2 #-0.17
    f2= -.1#-0.13
    time_coupler =2000 #ns
    if which_qubit == 0:  #readout on both qubits, excited q2
        Fs.spectroscopy_ge_coupler_switch(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit,time_coupler=time_coupler)
    if which_qubit == 1: ##readout on both qubits, excited q1
        Fs.spectroscopy_ge_coupler_switch(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit,time_coupler=time_coupler)
   


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
    num_steps = 101
    sweep_time =200 #ns
    if which_qubit == 0:
        Fs.rabi_ge(num_steps,sweep_time,ssm_geq2,ROIF1,ROIF2,which_qubit)
    if which_qubit == 1:
        Fs.rabi_ge(num_steps,sweep_time,ssm_geq1,ROIF1,ROIF2,which_qubit)

if run_vacuum_rabi:
    num_steps =101
    sweep_time =500 #ns #
    if which_qubit == 0:
        Fs.vacuum_rabi(num_steps,sweep_time,ssm_geq2,ROIF1,ROIF2,which_qubit,pi_ge_time)
    if which_qubit == 1:
        Fs.vacuum_rabi(num_steps,sweep_time,ssm_geq1,ROIF1,ROIF2,which_qubit,pi_ge_time)   
        
if run_rabi_ge_2qubit:
    num_steps = 81
    sweep_time =2000#200 #ns
    Fs.rabi_ge_2qubit(num_steps,sweep_time,ssm_geq1,ssm_geq2, ROIF1,ROIF2)
    

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
    

if run_rabi_2qubit:
    num_steps = 51
    sweep_time = 400
    Fs.rabi_ge_2qubit(num_steps,sweep_time,ssm_geq1,ssm_geq2,ROIF1, ROIF2)
    
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
    
if run_T1_ge_2q_RO:
    num_steps = 81
    sweep_time =80000#ns
    Fs.T1_ge_2q_RO(num_steps,sweep_time,ssm_ge,pi_ge_time,ROIF1, ROIF2,q=which_qubit,ifload = 1)
       
if run_T1_M_ge_2q_RO:
    num_steps = 81
    sweep_time =8000#ns
    Fs.T1_M_ge_2q_RO(num_steps,sweep_time,ssm_ge,pi_ge_time,ROIF1, ROIF2,q=which_qubit,ifload = 1)
 
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
    f1= -0.1 # -.025
    f2= -0.0 #-.04
    if which_qubit == 0:
#        Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
        Fs.parametric_coupling(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_start=f1,ssm_stop=f2,spec_amp=1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
        #parametric_coupling(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_start=-.15,ssm_stop=-.25,spec_amp=.5,ROIF1=0,ROIF2=0,q=0)
    if which_qubit == 1:
        Fs.parametric_coupling(num_steps,ssm_ge,pi_ge=pi_ge_time,ssm_start=f1,ssm_stop=f2,spec_amp=1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
        
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
        
if threshold_sweep:
    threshold_start = 500
    threshold_stop = 700
    threshold_step = 5
    num_threshold = int((threshold_stop- threshold_start)/threshold_step)
    
    prob_vs_pats_diag_array = np.zeros(num_threshold)
    for i in range(num_threshold): 
    
        num_steps = 3
        reps = 15000
        if which_qubit == 0:
            qubit_2_thr[1] = threshold_start + i*threshold_step
            Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|g>, keeps population in ground state
            Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,r=which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|e>, sends population to |e>
            #Fs.pi_nopi(off = 2,coef=1,coefpief=1,pi_ge=pi_ge_time,pi_ef=pi_ef_time,q= which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
            
        if which_qubit == 1:
            qubit_1_thr[1] = threshold_start + i*threshold_step
            Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,r=which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge); #|g>, keeps population in ground state
            Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,r=which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge); #|e>, sends population to |e>
            #Fs.pi_nopi(off = 2,coef=1,coefpief=1,pi_ge=pi_ge_time,pi_ef=pi_ef_time,q= which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge,ssm_ef=ssm_ef);
    
        wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
        
        if which_qubit == 1:
            prob_vs_pats_diag_array[i] = prob_vs_pats_1[0][0] + prob_vs_pats_1[1][1]
        if which_qubit == 0:
            prob_vs_pats_diag_array[i] = prob_vs_pats_2[0][0] + prob_vs_pats_2[1][1]
            
    thres_idx = (np.abs(prob_vs_pats_diag_array - 2)).argmin()
    print("done!")
    print("Threshold value closest to 2:", threshold_start + thres_idx*threshold_step)
            
    
        
    
if run_wx_sweep:
 
        sweep_steps_amp = 25#31
        
        wx_start = 0.5
        wx_stop = 2

        
        wx_ch2amp_steps = np.linspace(wx_start,wx_stop,sweep_steps_amp)
        #wx_ch1amp_steps = np.linspace(wx_start,wx_stop,sweep_steps_amp)
    
        if spec_ge_coupler_switch:
            xstart = f2
            xstop = f1
        if run_vacuum_rabi:
            xstart = 0
            xstop=sweep_time
        if pis_nopi:
            xstart = 0
            xstop=num_steps
        
            

        j=0
        out_I1, out_Q1, out_I2, out_Q2, out_P1, out_P2 = np.zeros( (6, sweep_steps_amp, num_steps))
        for j in range(len(wx_ch2amp_steps)):
            print('current amp:',wx_ch2amp_steps[j])
            print("amp j:",j)
            
            
            wx_amps[1] = wx_ch2amp_steps[j]
            wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
            
            n_vs_pats_1,n_vs_pats_2, rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
                                                                                                                                                                                                    

            out_Q1[j] = rec_avg_vs_pats_1[0];
            out_I1[j] = rec_avg_vs_pats_1[1];
            out_Q2[j] = rec_avg_vs_pats_2[0];
            out_I2[j] = rec_avg_vs_pats_2[1];
#            plt.plot(out_Q1[j]);plt.show()
#            plt.plot(out_I1[j]);plt.show()
#            plt.plot(out_Q2[j]);plt.show()
#            plt.plot(out_I2[j]);plt.show()
            
            p1_readout= np.matmul(scale_matrix1,prob_vs_pats_1)
            p2_readout= np.matmul(scale_matrix2,prob_vs_pats_2)
            
            plt.plot(rec_avg_vs_pats_1[0]);plt.title('Q Q1');plt.show()
            plt.plot(rec_avg_vs_pats_1[1]);plt.title('I Q1');plt.show()
            plt.plot(p1_readout[0]);plt.title('P Q1');plt.show()
            plt.plot(rec_avg_vs_pats_2[0]);plt.title('Q Q2');plt.show()
            plt.plot(rec_avg_vs_pats_2[1]);plt.title('I Q2');plt.show()
            plt.plot(p2_readout[0]);plt.title('P Q2');plt.show()

            
            out_P1[j] = p1_readout[0];
            out_P2[j] = p2_readout[0];
            
            print('current amp:',wx_ch2amp_steps[j])
            print("amp j:",j)
            print(j)
            print("probs vs pats for qubit 2")
            #print(prob_vs_pats_2)
              

           
            
            plt.imshow(out_Q1, extent=[xstart,xstop,wx_stop,wx_start],aspect='auto' );plt.show()
            plt.imshow(out_P1, extent=[xstart,xstop,wx_stop,wx_start],aspect='auto' );plt.show()
            plt.imshow(out_Q2, extent=[xstart,xstop,wx_stop,wx_start],aspect='auto' );plt.show()
            plt.imshow(out_P2, extent=[xstart,xstop,wx_stop,wx_start],aspect='auto' );plt.show()


        #save_basename = '\vacuumrabiq2_5000_ns_bnc_4.32_ssmq1ge_'+str(ssm_geq1)+'_ssmq2ge_'+str(ssm_geq2)+'_wxch2amp_'+str(wx_start)+'_to_'+str(wx_stop)+'steps_'+str(steps_in_seq)+'_wxch2off_coupler'+str(wx_offs[1])+'_.txt'
        save_basename = '\\wx_sweep4_vacuumrabi_500_detuning_13_mhz.txt'
        np.savetxt(save_dir+save_basename+"_I1",out_I1)
        np.savetxt(save_dir+save_basename+"_I2",out_I2)
        np.savetxt(save_dir+save_basename+"_Q1",out_Q1)
        np.savetxt(save_dir+save_basename+"_Q2",out_Q2)
        np.savetxt(save_dir+save_basename+"_P1",out_P1)
        np.savetxt(save_dir+save_basename+"_P2",out_P2)
        
#        delta_array = np.zeros(len(wx_ch2amp_steps))
#        for j in range(len(wx_ch2amp_steps)):
#            dummy_array_Q1 = out_Q1[j]
#            dummy_array_Q2 = out_Q2[j]
#            Q1_freq = wx_ch2amp_steps[np.abs(dummy_array_Q1).argmin()]
#            Q2_freq = wx_ch2amp_steps[np.abs(dummy_array_Q2).argmin()]
#            delta_array[j] = np.abs(Q1_freq - Q2_freq)
if sweep_ro_integration:
 
        sweep_steps_integration = 41#31
        integration_start = 100
        integration_stop = 4100
        
        integration_steps =np.linspace(integration_start,integration_stop,sweep_steps_integration)
    
        if spec_ge_coupler_switch:
            xstart = f2
            xstop = f1
        if run_vacuum_rabi:
            xstart = 0
            xstop=sweep_time
        if pis_nopi:
            xstart = 0
            xstop=num_steps
            
            
        j=0
        reps=15000
        out_rec_readout_vs_pats_1_A,out_rec_readout_vs_pats_1_B = np.zeros( (2, sweep_steps_integration, reps))
        for j in range(len(integration_steps)):
            print('current integration:',integration_steps[j])
            print("amp j:",j)
            
            ro_dur = np.int(integration_steps[j])
            wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
            
            n_vs_pats_1,n_vs_pats_2, rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
                                                                                                                                                                                                    
#            plt.plot(out_Q1[j]);plt.show()
#            plt.plot(out_I1[j]);plt.show()
#            plt.plot(out_Q2[j]);plt.show()
#            plt.plot(out_I2[j]);plt.show()
            
            import daq_processing
            daq_params = daq_programs_homo.get_daq_parameters(num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur)
            rec_readout_vs_pats_1 = daq_processing.record_vs_patterns(daq_params, rec_readout_1)
            out_rec_readout_vs_pats_1_A[j] = rec_readout_vs_pats_1[0,:,0]
            out_rec_readout_vs_pats_1_B[j] = rec_readout_vs_pats_1[0,:,1]
            

        save_basename = '\\integratpn_sweep5_0_mhz_detuning_piA_piB_ro_B.txt'
        np.savetxt(save_dir+save_basename+"A",out_rec_readout_vs_pats_1_A)
        np.savetxt(save_dir+save_basename+"B",out_rec_readout_vs_pats_1_B)





else:
    ##Do actual heterodyne
    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
    n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
    
#Correlation code
#calculate the readout correlation for <ZZ>
#    z1=np.ones([reps,num_steps])
#    z2=np.ones([reps,num_steps])
#    index = np.where(n_vs_pats_1[0,:,:]==0)
#    z1[index]=-1
#    index = np.where(n_vs_pats_2[0,:,:]==0)
#    z2[index]=-1
#    zz = np.mean(z1*z2,axis=0)
#    plt.plot(zz);plt.title('zz');plt.show()
    
    
    #rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = 0, deg_2 = 0,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur, verbose=True)
   
                                                                                                                                                                                                
    print("Qubit 1: P_Q1")
    P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
    I_Q1 = rec_avg_vs_pats_1[0];
    Q_Q1 = rec_avg_vs_pats_1[1];
    plt.plot(I_Q1);plt.title('I Q1');plt.show()
    plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
    
#Scaled Q1
    p1_readout= np.matmul(scale_matrix1,prob_vs_pats_1)
    plt.plot(p1_readout[0],label='|g>');plt.plot(p1_readout[1],label='|e>')
    plt.legend();plt.title('Q1 data scaled with matrix');plt.show()
    
    print("Qubit 2:")
    P_Q2 = prob_vs_pats_2[0];plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
    I_Q2 = rec_avg_vs_pats_2[0];
    Q_Q2 = rec_avg_vs_pats_2[1];
    plt.plot(I_Q2);plt.title('I Q2');plt.show()
    plt.plot(Q_Q2);plt.title('Q Q2');plt.show()  
    
#Scaled Q2
    p2_readout= np.matmul(scale_matrix2,prob_vs_pats_2)
    plt.plot(p2_readout[0],label='|g>');plt.plot(p2_readout[1],label='|e>')
    plt.legend();plt.title('Q2 data scaled with matrix');plt.show()
    

#    save_basename = "\\vacuumrabi_state01_3000reps.txt"
#    np.savetxt(save_dir+save_basename+"_P1",P_Q1)
#    np.savetxt(save_dir+save_basename+"_P2",P_Q2)
#    np.savetxt(save_dir+save_basename+"_P1_scaled_g",p1_readout[0])
#    np.savetxt(save_dir+save_basename+"_P1_scaled_e",p1_readout[1])
#    np.savetxt(save_dir+save_basename+"_P2_scaled_g",p2_readout[0])
#    np.savetxt(save_dir+save_basename+"_P2_scaled_e",p2_readout[1])
#    
#Entropy
    
#    entropy_meas =-p1_readout[0]*p2_readout[0]*np.log2(p1_readout[0]*p2_readout[0])-p1_readout[0]*p2_readout[1]*np.log2(p1_readout[0]*p2_readout[1])-p1_readout[1]*p2_readout[0]*np.log2(p1_readout[1]*p2_readout[0])-p1_readout[1]*p2_readout[1]*np.log2(p1_readout[1]*p2_readout[1])
#    plt.plot(entropy_meas)
#    plt.title('Joint Entropy');plt.show()
#    
#    entropy_meas_q1 =-p1_readout[0]*np.log2(p1_readout[0])-p1_readout[1]*np.log2(p1_readout[1])
#    plt.plot(entropy_meas_q1)
#    plt.title('Q1(B) Entropy');plt.show()
#    
#    entropy_meas_q2 =-p2_readout[0]*np.log2(p2_readout[0])-p2_readout[1]*np.log2(p2_readout[1])
#    plt.plot(entropy_meas_q2)
#    plt.title('Q2(A) Entropy');plt.show()
    

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
    

if run_rabi or run_vacuum_rabi:
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,y,guess_vals=[19,0.3,0.08,38,0.5])
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


    
if spec_ge:
    freq = np.linspace(f1,f2,num_steps)
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        freq_index = np.where(Q == np.min(Q))
        print("Q")
        plt.plot(freq,Q)
    if Irange>Qrange:
        freq_index = np.where(I == np.min(I))
        print("I")
        plt.plot(freq,I)
#    freq_index = np.where(y == np.max(y))
    ssm_ge = freq[freq_index]
    print(ssm_ge, freq_index)
#    plt.plot(qubit_bnc+freq,Q)
    
#    analysis.fit_lorentzian(freq,y,guess_vals=[2.6457e-10,-2.949e-10,0.5,-0.0423])#[10,2,2,1,-179]
if spec_ge_coupler_switch:
    freq = np.linspace(f1,f2,num_steps)
    Qrange = abs(np.max(Q)-np.min(Q))
    Irange = abs(np.max(I)-np.min(I))
    if Qrange>Irange:
        freq_index = np.where(Q == np.max(Q))
        print("Q")
        plt.plot(freq,Q)
    if Irange>Qrange:
        freq_index = np.where(I == np.max(I))
        print("I")
        plt.plot(freq,I)
#    freq_index = np.where(y == np.max(y))
    ssm_ge = freq[freq_index]
    print(ssm_ge, freq_index)

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

if run_T1_ge_2q_RO or run_T1_M_ge_2q_RO:
    times = np.linspace(0,sweep_time/1000,num_steps)
    Q1range = abs(np.max(Q_Q1)-np.min(Q_Q1))
    I1range = abs(np.max(I_Q1)-np.min(I_Q1))
    Q2range = abs(np.max(Q_Q2)-np.min(Q_Q2))
    I2range = abs(np.max(I_Q2)-np.min(I_Q2))
    #For Q1
    if Q1range>I1range:
        T1_ge_fit_vals_Q1,_,_,_ = analysis.fit_exp_decay(times,Q_Q1,guess_vals=[-400,0.001,-1700])
    if I1range>Q1range:
        T1_ge_fit_vals_Q1,_,_,_ = analysis.fit_exp_decay(times,I_Q1,guess_vals=[400,1,-1700])
    #Thresholding
#    T1_ge_fit_vals_Q1,_,_,_ = analysis.fit_exp_decay(times,P_Q1,guess_vals=[10,1,0.5])
    T1_ge_Q1 = 1/T1_ge_fit_vals_Q1[1]
    print("T1_ge_Q1 = {} \u03BCs".format(T1_ge_Q1))
    
    #For Q2
    if Q2range>I2range:
        T1_ge_fit_vals_Q2,_,_,_ = analysis.fit_exp_decay(times,Q_Q2,guess_vals=[-400,0.001,-1700])
    if I1range>Q1range:
        T1_ge_fit_vals_Q2,_,_,_ = analysis.fit_exp_decay(times,I_Q2,guess_vals=[400,1,-1700])
    #Thresholding
#    T1_ge_fit_vals_Q2,_,_,_ = analysis.fit_exp_decay(times,P_Q2,guess_vals=[10,1,0.5])
    T1_ge_Q2 = 1/T1_ge_fit_vals_Q2[1]
    print("T1_ge_Q2 = {} \u03BCs".format(T1_ge_Q2))
    
    
    
if piA_piB_roB:
    import daq_processing
    daq_params = daq_programs_homo.get_daq_parameters(num_patterns=num_steps, 
                                    num_records_per_pattern=reps,ro_dur=ro_dur)
    rec_readout_vs_pats_1 = daq_processing.record_vs_patterns(daq_params, rec_readout_1)
    
        
    
    
    #Histogram
    from scipy.stats import norm
    plt.figure(1, figsize=(8, 8))
    
    ax_hist1 = plt.axes()
    ax_hist1.tick_params(direction='in', labelbottom=False)
    binwidth =20
    
    pi_a_min = np.min([rec_readout_vs_pats_1[:,:,0][0]])
    pi_a_max = np.max([rec_readout_vs_pats_1[:,:,0][0]])
    bins_pi_a = np.arange(pi_a_min, pi_a_max, binwidth)
    
    hist_out_pi_a = ax_hist1.hist(rec_readout_vs_pats_1[:,:,0][0], bins=bins_pi_a, histtype='step', orientation='vertical', color = "firebrick",density = True)
    raw_data_a_y = np.sort(rec_readout_vs_pats_1[:,:,0][0])
    #note: for histtype, bar takes longer to run compared to step
    
    pi_b_min = np.min([rec_readout_vs_pats_1[:,:,1][0]])
    pi_b_max = np.max([rec_readout_vs_pats_1[:,:,1][0]])
    bins_pi_b = np.arange(pi_b_min, pi_b_max, binwidth)
    
    hist_out_pi_b= ax_hist1.hist(rec_readout_vs_pats_1[:,:,1][0], bins=bins_pi_b, histtype='step', orientation='vertical', color = "mediumpurple", density = True)
    raw_data_b_y = np.sort(rec_readout_vs_pats_1[:,:,1][0])
    
    mu_pi_a, std_pi_a = norm.fit(raw_data_a_y)
    gaussian_pi_a = norm.pdf(bins_pi_a, mu_pi_a,std_pi_a)
    
    mu_pi_b, std_pi_b = norm.fit(raw_data_b_y)
    gaussian_pi_b = norm.pdf(bins_pi_b, mu_pi_b,std_pi_b)
    ax_hist1.set_title('Pi_A Pi_B RO_B')
    
    plt.plot(bins_pi_a, gaussian_pi_a,'k', linewidth=2)
    plt.plot(bins_pi_b, gaussian_pi_b,'k', linewidth=2)
    plt.show()
    
    #is gaussian fit bad or good?
    print("is gaussian fit bad or good? input 'bad' or 'good' : ")
    fit_quality = input()
    
    if fit_quality == "good":
        print("cool")
    if fit_quality == "bad":
        #Find maxes of both histograms
        x_a_hist = np.where(hist_out_pi_a[0] == max(hist_out_pi_a[0]))[0]
        x_b_hist =np.where(hist_out_pi_b[0] == max(hist_out_pi_b[0]))[0]
        x_a_hist =x_a_hist[0]; x_b_hist=x_b_hist[0]
         
        #Finding which histogram has a smaller array
        if len(hist_out_pi_a[0]) != len(hist_out_pi_b[0]):
            smaller = min(len(hist_out_pi_a[0]),len(hist_out_pi_b[0]))
            bigger = max(len(hist_out_pi_a[0]),len(hist_out_pi_b[0]))
        if len(hist_out_pi_a[0]) < len(hist_out_pi_b[0]):
            print("smallest: hist_out_pi_a")
            small_array_freq = hist_out_pi_a[0]
            small_array_bins = hist_out_pi_a[1]
            big_array = hist_out_pi_b
        else:
            print("smallest: hist_out_pi_b")
            small_array_freq = hist_out_pi_b[0]
            small_array_bins = hist_out_pi_b[1]
            big_array = hist_out_pi_a
            
        #Add zeros to end of smaller array so we can compare both
        for i in range(bigger - smaller):
            small_array_freq = np.append(small_array_freq,0)
            small_array_bins = np.append(small_array_bins,small_array_bins[-1]+ binwidth)
            
        #Find index where intersection of histograms happens
        dummy_array = np.abs(small_array_freq[x_b_hist:x_a_hist] - big_array[0][x_b_hist:x_a_hist])
        idx = (np.abs(dummy_array - 0.0)).argmin()
        idx = idx + x_b_hist
        big_array_bin = big_array[1][idx]
        small_array_bin = small_array_bins[idx]
    
        if len(hist_out_pi_a[0]) == smaller:
            raw_idx_a = (np.abs(raw_data_a_y - small_array_bin)).argmin()
        else:
            raw_idx_a = (np.abs(raw_data_a_y - big_array_bin)).argmin()
        raw_data_a_y = raw_data_a_y[raw_idx_a:]
        mu_pi_a, std_pi_a = norm.fit(raw_data_a_y)
        print("mu_A = ",mu_pi_a )
        print("std_A = ",std_pi_a)
        gaussian_pi_a = norm.pdf(bins_pi_a,mu_pi_a,std_pi_a)
        
        if len(hist_out_pi_b[0]) == smaller:
            raw_idx_b = (np.abs(raw_data_b_y - small_array_bin)).argmin()
        else:
            raw_idx_b = (np.abs(raw_data_b_y - big_array_bin)).argmin()
        raw_data_b_y = raw_data_b_y[:raw_idx_b]
        mu_pi_b, std_pi_b = norm.fit(raw_data_b_y)
        print("mu_B = ",mu_pi_b )
        print("std_B = ",std_pi_b)
        gaussian_pi_b = norm.pdf(bins_pi_b, mu_pi_b,std_pi_b)
        
        
        #Plot histogram and Gaussian fits
        plt.figure(2, figsize=(8, 8))
        ax_hist2 = plt.axes()
        ax_hist2.tick_params(direction='in', labelbottom=False)
        hist_out_pi_a = ax_hist2.hist(rec_readout_vs_pats_1[:,:,0][0], bins=bins_pi_a, histtype='step', orientation='vertical', color = "firebrick",density = True)
        hist_out_pi_b= ax_hist2.hist(rec_readout_vs_pats_1[:,:,1][0], bins=bins_pi_b, histtype='step', orientation='vertical', color = "mediumpurple", density = True)
        ax_hist2.set_title('Pi_A Pi_B RO_B')
        plt.plot(bins_pi_a, gaussian_pi_a,'k', linewidth=2)
        plt.plot(bins_pi_b, gaussian_pi_b,'k', linewidth=2)
        plt.show()
     

            
        
        
        
    
    
    
    
    
#    save_basename = '\pi_a_pi_b_ro_b_binsize_' +str(binwidth)+'_RO_'+str(ro_dur)+'_.txt'
#    np.savetxt(save_dir+save_basename+"_pi_a_bins",hist_out_pi_a[1])
#    np.savetxt(save_dir+save_basename+"_pi_a_normalized_counts",hist_out_pi_a[0])
#    np.savetxt(save_dir+save_basename+"_pi_b_bins",hist_out_pi_b[1])
#    np.savetxt(save_dir+save_basename+"_pi_b_normalized_counts",hist_out_pi_b[0])
#    

    
    
    
    
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
    #num_steps=101#may need to remove eventually
    out_Q1, out_I1,  out_Q2, out_I2 = chevron.sweep_keithley(19.6,22,3,num_steps,reps,ro_dur,ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2,qubit_1_thr,qubit_2_thr)
    #save_basename = '\specq2_amp_4.32-.16,-.14_coupler_off'+str(wx_amps[1])+'coupler_off'+ str(wx_offs[1])+'keithley10.5,11.5'+'.txt'
    save_basename = '\calibration11.txt'
    np.savetxt(save_dir+save_basename+"_I1",out_I1)
    np.savetxt(save_dir+save_basename+"_I2",out_I2)
    np.savetxt(save_dir+save_basename+"_Q1",out_Q1)
    np.savetxt(save_dir+save_basename+"_Q2",out_Q2)
    #np.savetxt(save_dir+'\specq2_amp.1_-.2,-.05ssb_Q_keithleysweep_10,15mA,41steps_coupler0.1',out_Q)

#if run_J_sweep:
#    num_steps = 81
#    sweep_time = 1000
#    out_y = chevron.sweep_rabi_J(ssm_ef-0.005,ssm_ef+0.005,51,model,scale_matrix,Ax,By,num_steps,reps,ro_dur,sweep_time,pi_ge_time,pi_ef_time,ssm_ge,ssm_ef,IQ_angle)
##    np.savetxt(save_dir+'sweep_JEP_-0.02to0.05_51steps_rabi_ef_4us_51steps_10kavg_scaledmatrix',out_y)
#    np.savetxt(save_dir+'sweep_ssbef_chevron_10MHz_rabi_ef_0.1amp_1us_81steps_2kavg_scaledmatrix',out_y)

if run_ssbef_sweep: #chevron
    sweep_steps = 41
    num_steps = 51
    sweep_time = 1000
    start_freq = ssm_ef-0.002
    stop_freq = ssm_ef+0.002
    out_y = chevron.sweep_rabi_ef(start_freq, stop_freq,sweep_steps,model,scale_matrix,Ax,By,num_steps,reps,ro_dur,sweep_time,pi_ge_time,ssm_ge,ssm_ef,IQ_angle,which_qubit)
    np.savetxt(save_dir+'sweep_ssbef_chevron_4MHz_rabi_ef_0.08amp_1us_51steps_2kavg_scaledmatrix',out_y)

if run_bnc_sweep:
    sweep_steps = 11
    start_freq = qubit_bnc-.005
    stop_freq = qubit_bnc + 0.005
    out_Q, out_I = chevron.sweep_bnc_freq(start_freq, stop_freq,sweep_steps,num_steps,reps,ro_dur, IQ_angle)
    #np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q1_I',out_I)
    #np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q1_Q',out_Q)
    
if run_bnc_sweep_2q:
    sweep_steps = 51
    start_freq = 4.32- 0.005 #qubit_bnc - 0.005
    stop_freq = 4.32+ 0.005 #qubit_bnc + 0.005
    out_Q1, out_I1,out_Q2, out_I2 = chevron.sweep_bnc_freq_2q(start_freq, stop_freq,sweep_steps,num_steps,reps,ro_dur, ROIF1,ROIF2, deg_1 = 0, deg_2 = 0,verbose=True)
    np.savetxt(save_dir+'\coupler_p5_8usRabi_chev4'+str(start_freq)+'to'+str(stop_freq)+'_q1_I',out_I1)
    np.savetxt(save_dir+'\coupler_p5_8usRabi_chev4'+str(start_freq)+'to'+str(stop_freq)+'_q1_Q',out_Q1)
    np.savetxt(save_dir+'\coupler_p5_8usRabi_chev4'+str(start_freq)+'to'+str(stop_freq)+'_q2_I',out_I2)
    np.savetxt(save_dir+'\coupler_p5_8usRabi_chev4'+str(start_freq)+'to'+str(stop_freq)+'_q2_Q',out_Q2)
    
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