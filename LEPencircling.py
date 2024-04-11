# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:14:04 2020

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

#from generator import *
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
import copy_FPJPA_sequence as Fs
import wx_programs as wx
import daq_programs_homo
import analysis
import copy_chevron
import bnc
import NH_pi_nopi_RO as thr
from sklearn import svm
from sklearn.cluster import KMeans

instrument_path = r"C:\Users\Crow108\Documents\Python\instr\analyzer"
if instrument_path not in sys.path: sys.path.append(instrument_path )

analyzer_path = r"C:\Users\Crow108\Documents\Python\instr\python_interface\python_without_WX2184C"
if analyzer_path not in sys.path: sys.path.append(analyzer_path )
#save_dir = r"C:\Data\2023\two_qubit_tuna\chevron_sweeps"
save_dir = r"Z:\Weijian\tunable_dissipation/"

#readout bnc (BNC 2)
target_bnc_black_address = 'GPIB0::19::INSTR'

#qubit bnc (BNC 1)
target_bnc_address_6 = 'USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR' #835 with touchscreen 'USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR'

import datetime
from datetime import datetime

#import expt_parameters
#import seq_experiments
#import seq_programs
#import daq_programs
#import analysis
#import math
import black_nonHermitian_ma_V1
#import scipy as sp

#from Nop_class import Nop

#import wx_programs
#import keithley2401


    
reps = 500
ro_dur = 5000
## [RO,QC,Q1,Q2]
wx_amps =[1.99,0.9,.02,1.99] #0.52# [0.7,1.1,.55,.5]
wx_offs = [0,.5*2,0,0]
IQ_angle =160#250 #225

qubit_bnc =  4.32
RO_LO = 6.77#6.85 #6.6
ROq1 = 6.873182#6.873082#6.872626 
ROq2 =  6.68923#6.68913#6.68902
ssm_geq1 =-.1556#-.0054#-0.1528
ssm_geq2 = -0.1455-0.0004# -0.1538+0.00016 #-.162 #-0.1468
ssm_efq1 = -.2375
ssm_efq2 = -.249
#q1 = qubit(1)
#q1.ssbge
#q1.ROfre

which_qubit = 0 #qubit 1 = 1; qubit 2 = 0; both qubits = -1

if which_qubit ==1:
    pi_ge_time = 459#474#61
if which_qubit == 0:
    pi_ge_time = 44 #31#44#50#49#50 #q2
    
#qubit_1_thr= [-1500,300]#[-800,550]
#qubit_2_thr= [-1500,-500] #[-3000,-1200]#[500,1700]#
#IQ_angle_q1 = -90 #130
#IQ_angle_q2 = 240 - 60 +180 #increasing rotates CCW

qubit_1_thr= [-2500,600]
qubit_2_thr=[-2050,-250]
IQ_angle_q1 = 10
IQ_angle_q2 = 110#190 #increasing rotates CCW

#pi_ef_time = 40.46

# set 1 for which sequence to run
# qubit calibration
run_rabi = 0
run_ramsey = 0
run_T1 =  0
spec_ge=0
spec_ef= 0
run_tomo_calib = 0
run_rabi_raman = 0
run_rabi_raman_sweep_detun = 0
run_rabi_tomography = 0
run_T1_ge_tunable_dissipation = 0

# single qubit encircling
run_LEPencircling = 0
run_LEPencircling_raman_drive = 1
run_LEPencircling_tunable_raman_drive = 0

# sweep
run_bnc_sweep = 0
run_bnc_sweep_2q = 0
run_bnc_power_sweep = 0
run_bnc_pow_freq_sweep = 0
run_ssbef_sweep = 0 #chevron for ssb_ef

# two qubit
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
    bnc.set_bnc_output(RO_freq, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc,power_dBm=13,bnc_addr=target_bnc_address_6)
if which_qubit == 0:
#    scale_matrix = np.array([[9.66666667e-01, 3.73666667e-01, 1.31400000e-01],
#                             [3.25333333e-02, 5.64333333e-01, 2.29000000e-01],
#                             [8.00000000e-04, 6.20000000e-02, 6.39600000e-01]])
    
#    scale_matrix = np.array([[0.008   , 0.8884    , 0.],
#                          [0.992    , 0.1116    , 0.],
#                          [0.        , 0.        , 1        ]])
    scale_matrix = np.array([[0.9904    , 0.1262    , 0.],
                          [0.0096    , 0.8738    , 0.],
                          [0.        , 0.        , 1        ]])    

    scale_matrix = np.linalg.inv(scale_matrix)
    scale_matrix_ES = np.array([[0.99593333, 0.58973333, 0.18366667],
                                [0.00273333, 0.31706667, 0.16833333],
                                [0.00133333, 0.0932    , 0.648     ]])
    scale_matrix_ES = np.linalg.inv(scale_matrix_ES)
    ssm_ge = ssm_geq2
    ssm_ef = -0.25721#-0.25675#-0.2804#ssm_ge - 0.15# #-0.4091#-0.409 #ssm_ge - 0.15
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
    RO_freq = RO_LO
    #RO_freq = ROq2
   # qubit_bnc =4.3#4.256
    bnc.set_bnc_output(RO_freq, power_dBm=13, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc,power_dBm=13, bnc_addr=target_bnc_address_6)
    
if which_qubit == -1:
    #Do this for heterodyne
    #qubit_bnc = 4.3
#    ssm_geq1 = qubit_bnc - (4.312013 +-.090+0.00025)
#    ssm_geq2 = qubit_bnc - (4.256 +-0.1013-0.0012+0.0001)
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
    RO_freq = RO_LO
    bnc.set_bnc_output(RO_freq, bnc_addr=target_bnc_black_address)
    bnc.set_bnc_output(qubit_bnc, power_dBm=13,bnc_addr=target_bnc_address_6)
    
if run_rabi:
    num_steps = 51
    sweep_time =100 #ns
    if which_qubit == 0:
        Fs.rabi_ge(num_steps,sweep_time,ssm_geq2,ROIF1,ROIF2,which_qubit)
    if which_qubit == 1:
        Fs.rabi_ge(num_steps,sweep_time,ssm_geq1,ROIF1,ROIF2,which_qubit)
        
    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
    n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
                                                                                                                                                                                                    
    print("Qubit 1: P_Q1")
    P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
    I_Q1 = rec_avg_vs_pats_1[0];
    Q_Q1 = rec_avg_vs_pats_1[1];
    plt.plot(I_Q1);plt.title('I Q1');plt.show()
    plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
##    
    print("Qubit 2:")
    P_Q2 = prob_vs_pats_2[0];plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
    I_Q2 = rec_avg_vs_pats_2[0];
    Q_Q2 = rec_avg_vs_pats_2[1];
    plt.plot(I_Q2);plt.title('I Q2');plt.show()
    plt.plot(Q_Q2);plt.title('Q Q2');plt.show()   

if run_rabi_raman:
    
    num_steps = 51
    sweep_time = 4000
    num_records_per_pattern = 500

    times = np.linspace(0.0, sweep_time, num_steps)*1e-3

    # detuned drives
    detun = -0.03 #-0.015 
    ssm_coax = ROIF2 + detun
    amp_coax = 1
    ssm_Q2 = ssm_ge - detun - 0.0055  # excitation -detun; decay + detun
    amp_Q2 = 0.5 #0.54 
    
    black_nonHermitian_ma_V1.rabi_ge_withRamanDrive(num_steps=num_steps,sweep_time=sweep_time,ssm_ge=ssm_ge,ROIF1=ROIF1,ROIF2=ROIF2,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2, q=0)
    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
    
    P_Q2 = []
        
    n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=num_records_per_pattern,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
##    
#        print("Qubit 2:")
#        plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
    P_Q2 = prob_vs_pats_2[0]; 
#        I_Q2 = rec_avg_vs_pats_2[0];
#        Q_Q2 = rec_avg_vs_pats_2[1];
    
    p2_readout = np.matmul(scale_matrix,prob_vs_pats_2)    
    P_Q2 = p2_readout[1]        
    
    
    plt.plot(times,P_Q2,'b-')
    plt.xlabel("time (us)")
    plt.ylabel("pop_e")
    plt.ylim([0,1])
    plt.show()

if run_rabi_raman_sweep_detun:
    
    num_steps = 51
    sweep_time = 4000
    num_records_per_pattern = 500

    times = np.linspace(0.0, sweep_time, num_steps)*1e-3
    
    
    detun_add_start = 0.003 #-0.01# 0.026 #-0.02 #-0.005 #-0.003 # -0.1 #-0.125 #-0.2 # 0.15 #ssb_coax - 0.03#0.094 #ssb_coax - 0.03 
    detun_add_stop = 0.01#0.030 #0.094 #ssb_coax - 0.02 #0.0066
    sweep_vals = np.linspace(detun_add_start,detun_add_stop,sweep_steps)


    # detuned drives
    detun = -0.03 #-0.015 
    ssm_coax = ROIF2 + detun
    amp_coax = 1
    ssm_Q2 = ssm_ge + detun + 0.0025  # excitation -detun; decay + detun
    amp_Q2 = 0.56 #0.54 
    
    
    black_nonHermitian_ma_V1.rabi_ge_withRamanDrive(num_steps=num_steps,sweep_time=sweep_time,ssm_ge=ssm_ge,ROIF1=ROIF1,ROIF2=ROIF2,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2, q=0)
    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
    
    P_Q2 = []
        
    n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=num_records_per_pattern,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
##    
#        print("Qubit 2:")
#        plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
    P_Q2 = prob_vs_pats_2[0]; 
#        I_Q2 = rec_avg_vs_pats_2[0];
#        Q_Q2 = rec_avg_vs_pats_2[1];
    
    p2_readout = np.matmul(scale_matrix,prob_vs_pats_2)    
    P_Q2 = p2_readout[1]        
    
    
    plt.plot(times,P_Q2,'b-')
    plt.xlabel("time (us)")
    plt.ylabel("pop_e")
    plt.ylim([0,1])
    plt.show()

    
if run_rabi_tomography:


    num_patterns = 51 # 51 # 51 #101#201
    
    sweep_time = 100 #2000 #2000
    num_records_per_pattern = 500
    times = np.linspace(0.0, sweep_time, num_patterns)*1e-3
    num_avgs = 1
    
    amp_tom = 0.251 #0.51 # 0.513
    pi_tomo=pi_ge_time


    
    black_nonHermitian_ma_V1.rabi_tomography(num_steps=num_patterns,sweep_time=sweep_time,ssm_ge=ssm_ge,pi_ge=pi_ge_time, dur_pi2_ge=pi_ge_time/2,ro_phase=90,ROIF1=ROIF1,ROIF2=ROIF2,q=0,off_set=0*num_patterns)
    black_nonHermitian_ma_V1.rabi_tomography(num_steps=num_patterns,sweep_time=sweep_time,ssm_ge=ssm_ge,pi_ge=pi_ge_time, dur_pi2_ge=pi_ge_time/2,ro_phase=0,ROIF1=ROIF1,ROIF2=ROIF2,q=0,off_set=1*num_patterns)
    black_nonHermitian_ma_V1.rabi_tomography(num_steps=num_patterns,sweep_time=sweep_time,ssm_ge=ssm_ge,pi_ge=pi_ge_time, dur_pi2_ge=0,ro_phase=0,ROIF1=ROIF1,ROIF2=ROIF2,q=0,off_set=2*num_patterns)
 
    black_nonHermitian_ma_V1.loading(num_steps = 3*num_patterns)
    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
    
    p_readout = np.zeros(3*num_patterns)
    for k in range(num_avgs):
        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*3*seq.num_patterns,
#                                                                                   num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)

        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=3*num_patterns, num_records_per_pattern=num_records_per_pattern,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
           
                                                                                                                                                                                                    
#        print("Qubit 1: P_Q1")
#        P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
#        I_Q1 = rec_avg_vs_pats_1[0];
#        Q_Q1 = rec_avg_vs_pats_1[1];
#        plt.plot(I_Q1);plt.title('I Q1');plt.show()
#        plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
    ##    
#        print("Qubit 2:")
#        plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
        P_Q2 = prob_vs_pats_2[0]; 
#        I_Q2 = rec_avg_vs_pats_2[0];
#        Q_Q2 = rec_avg_vs_pats_2[1];
        

        if k is 0:
            p_readout = P_Q2
#                   
        else:
            p_readout += P_Q2

    p_readout = p_readout/num_avgs
#    p_post = analysis.p_readout_postselected(p_readout)
#    p_post=analysis.p_readout_postselected_pief(p_readout)

    N = num_patterns    
    x_p = 2*p_readout[0:N]-1 
    y_p = 1-2*p_readout[N:2*N]
    z_p = 1-2*p_readout[2*N:3*N]
    
    plt.plot(times,x_p,'b-')
    plt.plot(times,y_p,'r-')
    plt.plot(times,z_p,'k-')
    plt.xlabel("time (us)")
    plt.ylabel("x,y,z")
    plt.legend(["x","y","z"])
    plt.ylim([-1,1])
    plt.show()


if run_T1_ge_tunable_dissipation:
    
    detun = -0.03 #-0.01 #-0.005 # -0.015 
    ssb_coax = ROIF2 + detun #-0.249  #+ 0.0065
    amp_coax = 1
    ssb_Q2 = ssm_ge + detun  # excitation -detun; decay + detun
    amp_Q2 =  0.5 #0.5 #0.34 #0.17 # 0.54 

    num_steps = 51
    t1_time = 20000#50000
    sweep_steps = 1 # 41

    # sweep detuning
    detun_add_start = 0.001 #-0.01# 0.026 #-0.02 #-0.005 #-0.003 # -0.1 #-0.125 #-0.2 # 0.15 #ssb_coax - 0.03#0.094 #ssb_coax - 0.03 
    detun_add_stop = 0.01#0.030 #0.094 #ssb_coax - 0.02 #0.0066
    sweep_vals = np.linspace(detun_add_start,detun_add_stop,sweep_steps)
#
    reps = 5000

    # sweep ssb_coax 
#    ssb_start = ssb_Q2 - 0.05 #ssb_coax - 0.03#0.094 #ssb_coax - 0.03 
#    ssb_stop =  ssb_Q2 + 0.05 #0.094 #ssb_coax - 0.02 #0.0066
#    sweep_vals = np.linspace(ssb_start,ssb_stop,sweep_steps)
#    
##    sweep amp_Q2
#    amp_Q2_start = 0.62
#    amp_Q2_stop = 1
#    sweep_vals = np.linspace(amp_Q2_start,amp_Q2_stop,sweep_steps)

    
    T1_time = []
    out_Q = np.zeros( (sweep_steps, num_steps))
    pop_e = np.zeros( (sweep_steps, num_steps))
    p2_readout = []
    
    for i,val in enumerate(sweep_vals):
    
        print('current sweep val = ', val)
        print('current index = ', i)
        black_nonHermitian_ma_V1.T1_ge_tunable_dissipation_gaussian_pulse_vary_amp(num_steps= num_steps, sweeptime = t1_time,ssm_ge=ssm_ge,pi_ge_time=pi_ge_time,ssm_coax=ssb_coax, amp_coax=amp_coax, ssm_Q2=ssb_Q2 + val, amp_Q2=amp_Q2, ROIF=ROIF2, q=which_qubit, ifload = 1)
        
#        Fs.T1_ge_tunable_dissipation_gaussian_pulse(num_steps= num_steps, sweeptime = t1_time,ssm_ge=ssm_ge,pi_ge_time=pi_ge_time,ssm_coax=val, amp_coax=amp_coax, ssm_Q2=ssm_ge + freq_raman_bnc[i] - 4.5, amp_Q2=amp_Q2, q=which_qubit, ifload = 1)

        wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave



        # No thresholding
#        daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_pi,counts_pi,rec_readout,rec_avg_all,rec_all = daq_programs.run_daq(num_steps,reps,ro_dur,IQangle=100)
#        plt.plot(Q);plt.show();plt.plot(I);plt.show()
#        out_Q[i] = I
#        
#        times = np.linspace(0,t1_time/1000,num_steps)
#        T1_fit_vals,_,_,_ = analysis.fit_exp_decay(times,I,guess_vals=[50,0.6,-200])
#        T1 = 1/T1_fit_vals[1]
#        T1_time.append(T1)
#        print("T1 = {} \u03BCs".format(T1))
        
        # thresholding along I axis
        
#        daq_params, t_histo, p,a,b,n_readout = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=[-120, -83],num_patterns=num_steps, 
#                                                                                       num_records_per_pattern=reps,authr=0,fg=2,ro_dur=ro_dur,IQangle=100)
#        
#        m = daq_params.threshold
#        p_readout = p
#        y=p_readout[1] + p_readout[2]

        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
           
                                                                                                                                                                                                    
#        print("Qubit 1: P_Q1")
#        P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
#        I_Q1 = rec_avg_vs_pats_1[0];
#        Q_Q1 = rec_avg_vs_pats_1[1];
#        plt.plot(I_Q1);plt.title('I Q1');plt.show()
#        plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
    ##    
        print("Qubit 2:")
        P_Q2 = prob_vs_pats_2[0];plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
        I_Q2 = rec_avg_vs_pats_2[0];
        Q_Q2 = rec_avg_vs_pats_2[1];
        plt.plot(I_Q2);plt.title('I Q2');plt.show()
        plt.plot(Q_Q2);plt.title('Q Q2');plt.show()   
        
        p2_readout = np.matmul(scale_matrix,prob_vs_pats_2)    
        P_Q2 = p2_readout[1]

        
        pop_e[i] = P_Q2

        times = np.linspace(0,t1_time/1000,num_steps)
        try:
            T1_fit_vals,_,_,_ = analysis.fit_exp_decay(times,P_Q2,guess_vals=[1.1,0.5,0.2])
            T1 = 1/T1_fit_vals[1]
            T1_time.append(T1)
            print("T1 = {} \u03BCs".format(T1))
        except RuntimeError:
            print("error - curve fit failed")
            T1 = 1000
            T1_time.append(T1)
        
        # cluster thresholing
#        daq_params, t_histo, p,a,b,n_readout = daq_programs.run_daq_cluster_threshold(model,Ax=Ax,By=By,num_patterns=num_steps,
#                                                                                         num_records_per_pattern=reps,authr=0,fg=3,ro_dur=ro_dur,IQangle=IQ_angle)
#        p_readout = p
#        p_readout= np.matmul(scale_matrix,p)
#        plt.plot(p_readout[0],label='|g>');plt.plot(p_readout[1],label='|e>');plt.plot(p_readout[2],label='|f>')
#        plt.legend();plt.title('data scaled with matrix');plt.show()
#        y=p_readout[1]
#        plt.plot(y)
#
#        times = np.linspace(0,t1_time/1000,num_steps)
#        T1_fit_vals,_,_,_ = analysis.fit_exp_decay(times,y,guess_vals=[50,0.6,20])
#        T1 = 1/T1_fit_vals[1]
#        T1_time.append(T1)
#        print("T1 = {} \u03BCs".format(T1))
        
#        np.savetxt(save_dir+'ramsey_coax_drive_'+str(ssb_start)+','+str(ssb_stop)+'_coax_amp'+str(amp_coax)+'t2sweeptime'+str(t2_time)+'_fullsequence',out_I)
#        np.savetxt(save_dir+'ramsey_coax_drive_'+str(ssb_start)+','+str(ssb_stop)+'_coax_amp'+str(amp_coax)+'t2sweeptime'+str(t2_time)+'_T2_times',T2_time)
#        np.savetxt(save_dir+'ramsey_coax_drive_'+str(ssb_start)+','+str(ssb_stop)+'_coax_amp'+str(amp_coax)+'t2sweeptime'+str(t2_time)+'_freq_fit',freq_fit)
        
    plt.plot(sweep_vals,T1_time,'b-*',label='T1');
    plt.xlabel('additional detuning (GHz)')
    plt.ylim([0,40])
    plt.legend();plt.show()
    
#   plt.figure(figsize=(4, 3))
    plt.imshow(pop_e, extent=[0, t1_time/1000, detun_add_stop, detun_add_start],aspect='auto')
    plt.xlabel('time (us)')
    plt.ylabel('additional detun (GHz)')
    plt.colorbar()
    plt.show()

#    np.savetxt(save_dir+'t1_tunable_decay_4periods_detun-30MHz_ampCavity1.0_ampQubit0.5_qubitDetunAdd1MHz_pop_e_022324',pop_e)
#    np.savetxt(save_dir+'t1_tunable_decay_'+'detun_'+str(detun)+','+str(detun_add_start)+','+str(detun_add_stop)+'t1sweeptime'+str(t1_time)+'_ROamp_'+str(amp_coax)+'_Q2amp_'+str(amp_Q2)+'_pop_e'+'_020424',pop_e) #'_Q1_amp'+str(amp_Q2)+
#    np.savetxt(save_dir+'t1_coax_drive_'+str(ssb_start)+','+str(ssb_stop)+'_coax_amp'+str(amp_coax)+'t1sweeptime'+str(t1_time)+'_T1_times',T1_time)

# fit data for T1 decay
#    decay_rate = []; osc_freq = []
#    amplitude = []; offset = []
#
#    for i,val in enumerate(sweep_vals):
#        times = np.linspace(0,t1_time/1000,num_steps)
#        
#        if i < 151:
#        
#            try:
#                T1_fit_vals,_,_,_ = analysis.fit_exp_decay(times[1:-1],pop_e[i,1:-1],guess_vals=[1,0.1,0.05])
#                decay = T1_fit_vals[1]
#                amp_temp = T1_fit_vals[0]
#                offset_temp = T1_fit_vals[2]
#                freq = 0
#                decay_rate.append(decay)
#                osc_freq.append(freq)
#                amplitude.append(amp_temp)
#                offset.append(offset_temp)
#            except RuntimeError:
#                print("error - curve fit failed")
#                decay = 0
#                freq = 0
#                decay_rate.append(decay)
#                osc_freq.append(freq)
#                
#
#        else: 
#            
#            T1_fit_vals,_,_,_ = analysis.fit_sine_decay(times,pop_e[i],guess_vals=[0.8,0.2,2,1,0.5])
#            decay = T1_fit_vals[1];    
#            freq = T1_fit_vals[0];
#            decay_rate.append(decay); 
#            osc_freq.append(freq);
#        
##    plt.plot(sweep_vals,np.reciprocal(decay_rate),'b-s',label='T1 (us)');
##    plt.plot(sweep_vals,osc_freq,'r-o',label='freq (MHz)');
##    plt.xlabel('coax amp')
##    plt.ylim([0,15])
##    plt.legend();plt.show()
#
#    plt.plot(sweep_vals,decay_rate,'b-s',label='gamma');
#    plt.xlabel('detun')
#    plt.ylim([0,5])
#    plt.legend();plt.show()

    
    
if run_LEPencircling:

    # define the encircling loop (encircling type-I LEP) Two directions
    amp_min = 0 # 0.005 #0.01
    amp_max = 0.19 #  0.22 #0.12 #0.27 #0.22 #0.25
    freq_detun = 0.005

    num_patterns = 51 # 51 #101#201
    
    sweep_time = 10000 #2000 #2000
    num_records_per_pattern = 500
    times = np.linspace(0.0, sweep_time, num_patterns)*1e-3
    num_avgs = 10
    
    amp_tom = 0.251 #0.51 # 0.513
    pi_tomo=pi_ge_time
    
    # plus direction
    black_nonHermitian_ma_V1.encirclingLEP1(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,num_steps = num_patterns,phase_tomo=90,off_set=0*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
    
    black_nonHermitian_ma_V1.encirclingLEP1(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,num_steps = num_patterns,phase_tomo=0,off_set=1*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
        
    black_nonHermitian_ma_V1.encirclingLEP1(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,num_steps = num_patterns,phase_tomo=0,off_set=2*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)

    # minus direction
    freq_detun = -freq_detun
    black_nonHermitian_ma_V1.encirclingLEP1(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,num_steps = num_patterns,phase_tomo=90,off_set=3*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
    
    black_nonHermitian_ma_V1.encirclingLEP1(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,num_steps = num_patterns,phase_tomo=0,off_set=4*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
        
    black_nonHermitian_ma_V1.encirclingLEP1(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,num_steps = num_patterns,phase_tomo=0,off_set=5*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)


   
    black_nonHermitian_ma_V1.loading(num_steps = 2*3*num_patterns)
    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
    
    p_readout = np.zeros(2*3*num_patterns)
    for k in range(num_avgs):
        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*3*seq.num_patterns,
#                                                                                   num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)

        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=2*3*num_patterns, num_records_per_pattern=num_records_per_pattern,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
           
                                                                                                                                                                                                    
#        print("Qubit 1: P_Q1")
#        P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
#        I_Q1 = rec_avg_vs_pats_1[0];
#        Q_Q1 = rec_avg_vs_pats_1[1];
#        plt.plot(I_Q1);plt.title('I Q1');plt.show()
#        plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
    ##    
#        print("Qubit 2:")
#        plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
        P_Q2 = prob_vs_pats_2[0]; 
#        I_Q2 = rec_avg_vs_pats_2[0];
        Q_Q2 = rec_avg_vs_pats_2[1];
        

        if k is 0:
            p_readout = P_Q2
#                   
        else:
            p_readout += P_Q2

    p_readout = p_readout/num_avgs
#    p_post = analysis.p_readout_postselected(p_readout)
#    p_post=analysis.p_readout_postselected_pief(p_readout)

    N = num_patterns    
    x_p = 2*p_readout[0:N]-1 
    y_p = 1-2*p_readout[N:2*N]
    z_p = 1-2*p_readout[2*N:3*N]
    
    x_m = 2*p_readout[3*N:4*N]-1
    y_m = 1-2*p_readout[4*N:5*N]
    z_m = 1-2*p_readout[5*N:6*N]
    


    plt.plot(times,x_p,'b-')
    plt.plot(times,y_p,'r-')
    plt.plot(times,z_p,'k-')
    
    plt.plot(times,x_m,'b--')
    plt.plot(times,y_m,'r--')
    plt.plot(times,z_m,'k--')
    plt.xlabel("time (us)")
    plt.ylabel("x,y,z")
    plt.legend(["x(+)","y(+)","z(+)","x(-)","y(-)","z(-)"])
    plt.ylim([-1,1])
    plt.show()
    
#    # analysis: entropy (purity), chirality v.s. duration time T
#    R_p = np.sqrt(x_pief**2 + y_pief**2 + z_pief**2)
#    R_m = np.sqrt(x_pief_m**2 + y_pief_m**2 + z_pief_m**2)
#    entropy_p = -(1+R_p)/2*np.log2((1+R_p)/2)-(1-R_p)/2*np.log2((1-R_p)/2)  
#    entropy_m = -(1+R_m)/2*np.log2((1+R_m)/2)-(1-R_m)/2*np.log2((1-R_m)/2)  
#    
#    chirality = 1/2*np.sqrt((x_pief-x_pief_m)**2 + (y_pief-y_pief_m)**2 + (z_pief-z_pief_m)**2) # defined by trace distance
#    
#    plt.plot(seq.times,entropy_p,'b-')
#    plt.plot(seq.times,entropy_m,'r-')
#    plt.plot(seq.times,chirality,'k--')
#    plt.xlabel("duration time (us)")
#    plt.ylabel("entropy, chirality")
#    plt.legend(["entropy(+)","entropy(-)","chirality"])
#    plt.ylim([-0.1,1])
#    plt.show()



#    np.savetxt(save_dir+'tomography_encirclingLEPI_1us_ampmin0.0ampmax0.19freq5MHzTwoDirectionsPlusX_secondtry_10000repetitions_'+ add_t, p_readout)
#    np.savetxt(save_dir+'ztomo2usJ60to30radpusD5post', p_post)

if run_LEPencircling_raman_drive:
    
    

    # define the encircling loop (encircling type-I LEP) Two directions
    amp_min = 0 # 0.005 #0.01
    amp_max = 0.1 #  0.22 #0.12 #0.27 #0.22 #0.25
    freq_detun = 0.005 #0.005

    num_patterns = 5#51 # 51 #101#201
    
    sweep_time = 20000# 10000#10000 #2000 #2000
    num_records_per_pattern = 500
    times = np.linspace(0.0, sweep_time, num_patterns)*1e-3
    num_avgs = 1
    
    amp_tom = 0.251 #0.251 # 0.513
    pi_tomo=pi_ge_time
    
#    # detuned drives (cooling)
    detun = -0.03 
    ssm_coax = ROIF2 + detun
    amp_coax = 1
    ssm_Q2 = ssm_ge + detun # +0.001 for amp_Q2=0.5 # excitation -detun; decay + detun
    amp_Q2 = 0.5 #0.5 #0.54 

    # detuned drives (heating)
#    detun = -0.03  
#    ssm_coax = ROIF2 + detun     
#    amp_coax = 1 
#    ssm_Q2 = ssm_ge - detun - 0.0055  # excitation -detun; decay + detun
#    amp_Q2 =  0.6
    
    # plus direction
    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=90,off_set=0*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
        
    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=0,off_set=1*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
        
    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=0,off_set=2*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
#    # minus direction
    freq_detun = -freq_detun
    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=90,off_set=3*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=0,off_set=4*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
        
    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=0,off_set=5*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)


   
    black_nonHermitian_ma_V1.loading(num_steps = 2*3*num_patterns)
    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
    
    p_readout = np.zeros(2*3*num_patterns)
    for k in range(num_avgs):
        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*3*seq.num_patterns,
#                                                                                   num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)

        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=2*3*num_patterns, num_records_per_pattern=num_records_per_pattern,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
           
                                                                                                                                                                                                    
#        print("Qubit 1: P_Q1")
#        P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
#        I_Q1 = rec_avg_vs_pats_1[0];
#        Q_Q1 = rec_avg_vs_pats_1[1];
#        plt.plot(I_Q1);plt.title('I Q1');plt.show()
#        plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
    ##    

#        print("Qubit 2:")
#        plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
        P_Q2 = prob_vs_pats_2[0]; 
#        I_Q2 = rec_avg_vs_pats_2[0];
        Q_Q2 = rec_avg_vs_pats_2[1];

        p2_readout = np.matmul(scale_matrix,prob_vs_pats_2)    
        P_Q2 = p2_readout[1]
        
#    plt.plot(p2_readout[0],label='|g>');plt.plot(p2_readout[1],label='|e>')
#    plt.legend();plt.title('Q2 data scaled with matrix');plt.show()


        if k is 0:
            p_readout = P_Q2
#                   
        else:
            p_readout += P_Q2

    p_readout = p_readout/num_avgs
#    p_post = analysis.p_readout_postselected(p_readout)
#    p_post=analysis.p_readout_postselected_pief(p_readout)

#    N = num_patterns    
#    x_p = 2*p_readout[0:N]-1 
#    x_m = 2*p_readout[N:2*N]-1
#
#
#    plt.plot(times,x_p,'b-')    
#    plt.plot(times,x_m,'b--')
#    plt.xlabel("time ($\\mu$s)")
#    plt.legend(["x (+)","x (-)"])
#    plt.ylim([-1,1])
#    plt.show()


    N = num_patterns    
    x_p = 2*p_readout[0:N]-1 
    y_p = 1-2*p_readout[N:2*N]
    z_p = 1-2*p_readout[2*N:3*N]
    
    x_m = 2*p_readout[3*N:4*N]-1
    y_m = 1-2*p_readout[4*N:5*N]
    z_m = 1-2*p_readout[5*N:6*N]
    


    plt.plot(times,x_p,'b-')
    plt.plot(times,y_p,'r-')
    plt.plot(times,z_p,'k-')
    plt.xlabel("time (us)")
    plt.ylabel("x,y,z")
    plt.legend(["x(+)","y(+)","z(+)"])
    plt.ylim([-1,1])
    plt.show()

    
    plt.plot(times,x_m,'b-')
    plt.plot(times,y_m,'r-')
    plt.plot(times,z_m,'k-')
    plt.xlabel("time (us)")
    plt.ylabel("x,y,z")
    plt.legend(["x(-)","y(-)","z(-)"])
    plt.ylim([-1,1])
    plt.show()
    


#    # analysis: entropy (purity), chirality v.s. duration time T
#    R_p = np.sqrt(x_pief**2 + y_pief**2 + z_pief**2)
#    R_m = np.sqrt(x_pief_m**2 + y_pief_m**2 + z_pief_m**2)
#    entropy_p = -(1+R_p)/2*np.log2((1+R_p)/2)-(1-R_p)/2*np.log2((1-R_p)/2)  
#    entropy_m = -(1+R_m)/2*np.log2((1+R_m)/2)-(1-R_m)/2*np.log2((1-R_m)/2)  
#    
#    chirality = 1/2*np.sqrt((x_pief-x_pief_m)**2 + (y_pief-y_pief_m)**2 + (z_pief-z_pief_m)**2) # defined by trace distance
#    
#    plt.plot(seq.times,entropy_p,'b-')
#    plt.plot(seq.times,entropy_m,'r-')
#    plt.plot(seq.times,chirality,'k--')
#    plt.xlabel("duration time (us)")
#    plt.ylabel("entropy, chirality")
#    plt.legend(["entropy(+)","entropy(-)","chirality"])
#    plt.ylim([-0.1,1])
#    plt.show()



#    np.savetxt(save_dir+'tomography_encirclingLEPI_20us_ampmin0.0ampmax0.1freq5MHzTwoDirectionsPlusX_5000repetitions_022924', p_readout)
#    np.savetxt(save_dir+'ztomo2usJ60to30radpusD5post', p_post)

if run_LEPencircling_tunable_raman_drive:
    
    

    # define the encircling loop (encircling type-I LEP) Two directions
    amp_min = 0 # 0.005 #0.01
    amp_max = 0.19 #  0.22 #0.12 #0.27 #0.22 #0.25
    freq_detun = 0.005 #0.005

    num_patterns = 51#51 # 51 #101#201
    
    sweep_time = 10000#100000# 10000#10000 #2000 #2000
    num_records_per_pattern = 500
    times = np.linspace(0.0, sweep_time, num_patterns)*1e-3
    num_avgs = 1
    
    amp_tom = 0.251 #0.51 # 0.513
    pi_tomo=pi_ge_time
    
    # detuned drives
    
    detun = -0.03 
    ssm_coax = ROIF2 + detun
    amp_coax = 1
    ssm_Q2 = ssm_ge + detun + 0.001 # +0.001 for amp_Q2=0.5 # excitation -detun; decay + detun
    amp_Q2 = 0.5 #0.5 #0.54 
    
    # plus direction
    black_nonHermitian_ma_V1.encirclingLEP1_tunable_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=90,off_set=0*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
        
#    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
#                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=0,off_set=1*num_patterns,
#                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#        
#    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
#                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=0,off_set=2*num_patterns,
#                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
#    # minus direction
    freq_detun = -freq_detun
    black_nonHermitian_ma_V1.encirclingLEP1_tunable_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=90,off_set=1*num_patterns,
                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
#                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=0,off_set=4*num_patterns,
#                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=pi_ge_time/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#        
#    black_nonHermitian_ma_V1.encirclingLEP1_raman_drive(rabi_time = sweep_time,pi_ge=pi_ge_time,ssm_ge =ssm_ge, ROIF=ROIF2
#                                       ,q=0,ssm_coax=ssm_coax, amp_coax=amp_coax, ssm_Q2=ssm_Q2, amp_Q2=amp_Q2,num_steps = num_patterns,phase_tomo=0,off_set=5*num_patterns,
#                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)


   
    black_nonHermitian_ma_V1.loading(num_steps = 2*1*num_patterns)
    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
    
    p_readout = np.zeros(2*1*num_patterns)
    for k in range(num_avgs):
        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*3*seq.num_patterns,
#                                                                                   num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)

        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=2*1*num_patterns, num_records_per_pattern=num_records_per_pattern,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
           
                                                                                                                                                                                                    
#        print("Qubit 1: P_Q1")
#        P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
#        I_Q1 = rec_avg_vs_pats_1[0];
#        Q_Q1 = rec_avg_vs_pats_1[1];
#        plt.plot(I_Q1);plt.title('I Q1');plt.show()
#        plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
    ##    

#        print("Qubit 2:")
#        plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
        P_Q2 = prob_vs_pats_2[0]; 
#        I_Q2 = rec_avg_vs_pats_2[0];
        Q_Q2 = rec_avg_vs_pats_2[1];

        p2_readout = np.matmul(scale_matrix,prob_vs_pats_2)    
        P_Q2 = p2_readout[1]
        
#    plt.plot(p2_readout[0],label='|g>');plt.plot(p2_readout[1],label='|e>')
#    plt.legend();plt.title('Q2 data scaled with matrix');plt.show()


        if k is 0:
            p_readout = P_Q2
#                   
        else:
            p_readout += P_Q2

    p_readout = p_readout/num_avgs
#    p_post = analysis.p_readout_postselected(p_readout)
#    p_post=analysis.p_readout_postselected_pief(p_readout)

    N = num_patterns    
    x_p = 2*p_readout[0:N]-1 
    x_m = 2*p_readout[N:2*N]-1


    plt.plot(times,x_p,'b-')    
    plt.plot(times,x_m,'b--')
    plt.xlabel("time ($\\mu$s)")
    plt.legend(["x (+)","x (-)"])
    plt.ylim([-1,1])
    plt.show()


#    N = num_patterns    
#    x_p = 2*p_readout[0:N]-1 
#    y_p = 1-2*p_readout[N:2*N]
#    z_p = 1-2*p_readout[2*N:3*N]
#    
#    x_m = 2*p_readout[3*N:4*N]-1
#    y_m = 1-2*p_readout[4*N:5*N]
#    z_m = 1-2*p_readout[5*N:6*N]
#    
#
#
#    plt.plot(times,x_p,'b-')
#    plt.plot(times,y_p,'r-')
#    plt.plot(times,z_p,'k-')
#    plt.xlabel("time (us)")
#    plt.ylabel("x,y,z")
#    plt.legend(["x(+)","y(+)","z(+)"])
#    plt.ylim([-1,1])
#    plt.show()
#
#    
#    plt.plot(times,x_m,'b-')
#    plt.plot(times,y_m,'r-')
#    plt.plot(times,z_m,'k-')
#    plt.xlabel("time (us)")
#    plt.ylabel("x,y,z")
#    plt.legend(["x(-)","y(-)","z(-)"])
#    plt.ylim([-1,1])
#    plt.show()
    


#    # analysis: entropy (purity), chirality v.s. duration time T
#    R_p = np.sqrt(x_pief**2 + y_pief**2 + z_pief**2)
#    R_m = np.sqrt(x_pief_m**2 + y_pief_m**2 + z_pief_m**2)
#    entropy_p = -(1+R_p)/2*np.log2((1+R_p)/2)-(1-R_p)/2*np.log2((1-R_p)/2)  
#    entropy_m = -(1+R_m)/2*np.log2((1+R_m)/2)-(1-R_m)/2*np.log2((1-R_m)/2)  
#    
#    chirality = 1/2*np.sqrt((x_pief-x_pief_m)**2 + (y_pief-y_pief_m)**2 + (z_pief-z_pief_m)**2) # defined by trace distance
#    
#    plt.plot(seq.times,entropy_p,'b-')
#    plt.plot(seq.times,entropy_m,'r-')
#    plt.plot(seq.times,chirality,'k--')
#    plt.xlabel("duration time (us)")
#    plt.ylabel("entropy, chirality")
#    plt.legend(["entropy(+)","entropy(-)","chirality"])
#    plt.ylim([-0.1,1])
#    plt.show()



#    np.savetxt(save_dir+'tomography_encirclingLEPI_1us_ampmin0.0ampmax0.19freq5MHzTwoDirectionsPlusX_secondtry_10000repetitions_'+ add_t, p_readout)
#    np.savetxt(save_dir+'ztomo2usJ60to30radpusD5post', p_post)




#   #########################################################

#    ssb_ge = 0.4045+0.0001 + 0.0047 # left side of the bath engineered spectrum
#    ssb_ef = 0.1101+0.00025 +0.0032
#    p_ge = 39
#    p_ef = 32

    ###########################################################
#    # define the encircling loop (encircling NHH EP) Two directions
#    # EP encircling tomography with error bands 08/19/2021
#
#    
#    save_dir = r"Z:\Weijian\NHH_Jumps\slow_driving/"
#    amp_min = 0.01 #0.01
#    amp_min = 1
#    amp_max = 0.27 #0.25 #0.25
#    freq_detun = -0.015
#
#    seq.comment = "measurement"
#    seq.num_patterns = 101#201#51
#    
##    0.3
#    seq.sweep_time = 4000 #2000
#    seq.num_records_per_pattern = 2000
#    seq.times = np.linspace(0.0, seq.sweep_time, seq.num_patterns)*1e-3
#    seq.num_avgs= 5
#    
#    sweep.comment = "rabi ssb freq near EP"
#    sweep.steps = 10
#    sweep.vals = np.linspace(0,0.3, sweep.steps)*1e-3 - 0.1202 # ge ramsey 0.3925
##    sweep.vals = np.linspace(-0.5,0.5, sweep.steps)*1e-3 + 0.0856  #+ 0.0915 # ef
##    sweep.vals = np.linspace(0.37,0.38,sweep.steps) # wx voltage
##    sweep.vals = np.linspace(90,90,sweep.steps) # mixer nonorthogonality
##    sweep.vals = np.linspace(0.043,0.045,sweep.steps)
##    sweep.vals = np.linspace(0.38,0.40,sweep.steps) # search eigenstate
##    sweep.vals = np.linspace(0,1,sweep.steps)
#
#    
#    N = seq.num_patterns    
#
#    amp_tom = 0.514
#    ph_ang = 0
#
#    pop_3state = []
#    pop_post = []
#    x_tom = []
#    y_tom = []
#    z_tom = []
#
#
#    pi_tomo=p_ef
#    
#    for idx, _ in enumerate(sweep.vals):
#        
#        print(idx)
#        threshold = [-166.05195052355808, -141.55195052355808];
#    
##        # plus direction
##        black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                           ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=0*seq.num_patterns,
##                                           amp_tomo=amp_tom,pi_time=0,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
##        
##        black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                           ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=1*seq.num_patterns,
##                                           amp_tomo=amp_tom,pi_time=0,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
##        
##        black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                           ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=2*seq.num_patterns,
##                                           amp_tomo=amp_tom,pi_time=0,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
##        
##        black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                           ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=3*seq.num_patterns,
##                                           amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
##        
##        black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                           ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=4*seq.num_patterns,
##                                           amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
##        
##        black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                           ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=5*seq.num_patterns,
##                                           amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
##        
##        
##        #    # minus direction
##        #    freq_detun = -freq_detun
##        #    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##        #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=3*seq.num_patterns,
##        #                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##        #    
##        #    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##        #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=4*seq.num_patterns,
##        #                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##        #
##        #    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##        #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=5*seq.num_patterns,
##        #                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##        
##        #    pi_tomo=0
##        #    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##        #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,
##        #                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
##        #
##        #    pi_tomo=p_ef
##        #    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##        #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=seq.num_patterns,
##        #                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
##           
##        black_nonHermitian_ma_V1.loading(num_steps = 2*3*seq.num_patterns)
##        
#        
#        
##        [-0.076, 0.025 , 0.09, -0.086]
#        wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5], offset=[-0.076, 0.025 , 0.09, -0.086])  # .033,.01, 0.033, -0.077; offset=[-0.062, 0.016 , 0.09, -0.086]
#        
#        #    threshold = [-165.44885667419695, -140.69885667419695];    
#        
#        for k in range(seq.num_avgs):
#            
#            daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*3*seq.num_patterns,
#                                                                                       num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#            
#        
#            if k is 0:
#                p_readout = p
#        #                   
#            else:
#                p_readout += p
#        
#        p_readout = p_readout/seq.num_avgs
#        #    p_post = analysis.p_readout_postselected(p_readout)
#        p_post=analysis.p_readout_postselected_pief(p_readout)
#        
#        
#        x_pief = p_post[2][0:N]-p_post[1][0:N] # with pi_ef pusle at the end
#        y_pief = p_post[1][N:2*N]-p_post[2][N:2*N]
#        z_pief = p_post[1][2*N:3*N]-p_post[2][2*N:3*N]
#        
#        
#        pop_3state.append(p_readout)
#        pop_post.append(p_post)
#        x_tom.append(x_pief)
#        y_tom.append(y_pief)
#        z_tom.append(z_pief)    
#
#    
##    x_pief_m = 2*(p_readout[2][3*N:4*N]+p_readout[1][3*N:4*N])-1 # with pi_ef pusle at the end
##    y_pief_m = 1-2*(p_readout[2][4*N:5*N]+p_readout[1][4*N:5*N])
##    z_pief_m = 1-2*(p_readout[2][5*N:6*N]+p_readout[1][5*N:6*N])
#    
#    
#        plt.semilogy(seq.times,1-p_readout[0,0:N]) # submanifold population
#        plt.xlabel("time (us)")
#        plt.ylabel("pop_sub")
#        plt.show()
#
#
#        plt.plot(seq.times,x_pief,'b-') # plus direction (+)
#        plt.plot(seq.times,y_pief,'r-')
#        plt.plot(seq.times,z_pief,'k-')
#        
#    #    plt.plot(seq.times,x_pief_m,'b--') # minus  direction (-)
#    #    plt.plot(seq.times,y_pief_m,'r--')
#    #    plt.plot(seq.times,z_pief_m,'k--')
#        plt.xlabel("time (us)")
#        plt.ylabel("x,y,z")
#    #    plt.legend(["x(+)","y(+)","z(+)"])
#        plt.ylim([-1,1])
#        plt.show()
#        
#        # analysis: entropy (purity), chirality v.s. duration time T
#        R_p = np.sqrt(x_pief**2 + y_pief**2 + z_pief**2)
#    #    R_m = np.sqrt(x_pief_m**2 + y_pief_m**2 + z_pief_m**2)
#        entropy_p = -(1+R_p)/2*np.log2((1+R_p)/2)-(1-R_p)/2*np.log2((1-R_p)/2)  
#    #    entropy_m = -(1+R_m)/2*np.log2((1+R_m)/2)-(1-R_m)/2*np.log2((1-R_m)/2)  
#        
#    #    chirality = 1/2*np.sqrt((x_pief-x_pief_m)**2 + (y_pief-y_pief_m)**2 + (z_pief-z_pief_m)**2) # defined by trace distance
#        
#        plt.plot(seq.times,entropy_p,'b-')
#    #    plt.plot(seq.times,entropy_m,'r-')
#    #    plt.plot(seq.times,chirality,'k--')
#        plt.xlabel("duration time (us)")
#        plt.ylabel("entropy")
#    #    plt.legend(["entropy(+)","entropy(-)","chirality"])
#        plt.ylim([-0.1,1])
#        plt.show()
#
#
#
#    pop_3state_stack=np.dstack(pop_3state) # dimension: 3 (g,e,f) x number of patterns x sweep steps
#    pop_post_stack = np.dstack(pop_post)
#    x_tom_stack = np.stack(x_tom)
#    y_tom_stack = np.stack(y_tom)
#    z_tom_stack = np.stack(z_tom)
#    
##    pop_post1_stack = np.dstack(pop_post1)
##    msmt_fit_stack = np.stack(msmt_fit)
#    pop_3state_stack_reshape = pop_3state_stack.reshape(seq.num_patterns*6*3,sweep.steps)
#    pop_post_stack_reshape = pop_post_stack.reshape(int(seq.num_patterns*6*3/2),sweep.steps)
##    pop3state_f=popdstack.reshape(seq.num_patterns*3,steps)
##    np.savetxt('pop_3state_stack', pop_3state_stack)
##    np.savetxt('pop_post_stack', pop_post_stack)
#    
#    fig = plt.figure(figsize=(6,4))
#    plt.imshow(np.transpose(pop_post_stack[1,:,:]-pop_post_stack[2,:,:]),extent=[min(seq.times),max(seq.times),max(sweep.vals),min(sweep.vals)],cmap='bwr',aspect='auto')
##    plt.imshow(np.transpose(pop_post_stack[1,:,:]-pop_post_stack[2,:,:]),cmap='bwr',aspect='auto')
#    plt.xlabel("time")
#    plt.ylabel("ssm_ef")
#    plt.colorbar()
#    plt.show()
#    
#    
#    x_tom_mean = np.mean(x_tom_stack, axis=0)
#    y_tom_mean = np.mean(y_tom_stack, axis=0)
#    z_tom_mean = np.mean(z_tom_stack, axis=0)
#
#    x_tom_std = np.std(x_tom_stack, axis=0)
#    y_tom_std = np.std(y_tom_stack, axis=0)
#    z_tom_std = np.std(z_tom_stack, axis=0)
#
#    plt.plot(seq.times, z_tom_mean,'k-')
#    plt.plot(seq.times, x_tom_mean,'b-')
#    plt.plot(seq.times, y_tom_mean,'r-')
#    plt.plot(seq.times, z_tom_mean + z_tom_std,'k--')
#    plt.plot(seq.times, x_tom_mean + x_tom_std,'b--')
#    plt.plot(seq.times, y_tom_mean + y_tom_std,'r--')
#    plt.plot(seq.times, z_tom_mean - z_tom_std,'k--')
#    plt.plot(seq.times, x_tom_mean - x_tom_std,'b--')
#    plt.plot(seq.times, y_tom_mean - y_tom_std,'r--')
#    plt.xlabel('time (us)')
#    plt.ylabel('x, y, z')
#    plt.ylim([-1,1])
#    plt.show()
#
#    plt.plot(seq.times, z_tom_std,'k-')
#    plt.plot(seq.times, x_tom_std,'b-')
#    plt.plot(seq.times, y_tom_std,'r-')
#    plt.xlabel('time (us)')
#    plt.ylabel('std of x,y,z')
#    plt.show()
#
#
#
#
##    plt.plot(sigma[0:101])
##    plt.plot(sigma[101:202])
##    plt.plot(sigma[202:303])
##    plt.shows()
#    
##    for ii in range(7):
##        msmt.popt, msmt.perr, _, _ = analysis.fit_sine_decay(seq.times, pop_post_stack[2,:,ii], guess_vals=[0.35,3,0.12,30,0.75])
##        msmt_fit.append(msmt.popt)
#
##    plt.plot(sweep.vals,abs(msmt_fit_stack[:,0]),'b-s')
##    plt.plot(sweep.vals,msmt_fit_stack[:,1],'r-o')
##    plt.xlabel("ssm_ef")
##    plt.ylabel("frequency (MHz), decay rate (us-1)")
##    plt.show()
#
###    
##    fig = plt.figure()
##    plt.plot(sweep.vals,abs(msmt_fit_stack[:,0]),'ro')
###    plt.plot(sweep.vals,3.629*sweep.vals-0.003)
##    analysis.fit_parabola(sweep.vals, abs(msmt_fit_stack[:,0]), guess_vals=None)
###    plt.xlabel("amp")
##    plt.xlabel("amp")
##    plt.ylabel("Frequency (MHz)")
##    plt.show
#    
#    np.savetxt(save_dir+'pop3state_slow_driving_tomography_4us_J30_detun15MHz_minus direction_ssm_ef -0.1202_10000repetitions_10times_amp0.27_08192021',pop_3state_stack_reshape)
#    np.savetxt(save_dir+'pop_post_slow_driving_tomography_4us_J30_detun15MHz_minus direction_ssm_ef -0.1202_10000repetitions_10times_amp_0.27_08192021',pop_post_stack_reshape)
########    np.savetxt(save_dir+'pop_post fit_pop f_2us_WX 1.5_amp 0.05_ssm_ef 0.0895 0.0935_07252020',msmt_fit_stack)   
######           
#    np.savetxt(save_dir+'slow_driving_tomographyX_10000_repetitions_10times_ssm_ef -0.1202_amp0.27_08192021',x_tom_stack)   
#    np.savetxt(save_dir+'slow_driving_tomographyY_10000_repetitions_10times_ssm_ef -0.1202_amp0.27_08192021',y_tom_stack)   
#    np.savetxt(save_dir+'slow_driving_tomographyZ_10000_repetitions_10times_ssm_ef -0.1202_amp0.27_08192021',z_tom_stack)   
##
#
#    ###########################################################
#    # define the encircling loop (encircling NHH EP) Two directions
#    # EP encircling tomography with error bands 08/19/2021
#
#    
#    save_dir = r"Z:\Weijian\NHH_Jumps\slow_driving/"
#    amp_min = 0.01 #0.01
#    amp_min = 1
#    amp_max = 0.28 #0.25 #0.25
#    freq_detun = -0.015
#
#    seq.comment = "measurement"
#    seq.num_patterns = 101#201#51
#    
##    0.3
#    seq.sweep_time = 4000 #2000
#    seq.num_records_per_pattern = 2000
#    seq.times = np.linspace(0.0, seq.sweep_time, seq.num_patterns)*1e-3
#    seq.num_avgs= 5
#    
#    sweep.comment = "rabi ssb freq near EP"
#    sweep.steps = 10
#    sweep.vals = np.linspace(0,0.3, sweep.steps)*1e-3 - 0.1202 # ge ramsey 0.3925
##    sweep.vals = np.linspace(-0.5,0.5, sweep.steps)*1e-3 + 0.0856  #+ 0.0915 # ef
##    sweep.vals = np.linspace(0.37,0.38,sweep.steps) # wx voltage
##    sweep.vals = np.linspace(90,90,sweep.steps) # mixer nonorthogonality
##    sweep.vals = np.linspace(0.043,0.045,sweep.steps)
##    sweep.vals = np.linspace(0.38,0.40,sweep.steps) # search eigenstate
##    sweep.vals = np.linspace(0,1,sweep.steps)
#
#    
#    N = seq.num_patterns    
#
#    amp_tom = 0.514
#    ph_ang = 0
#
#    pop_3state = []
#    pop_post = []
#    x_tom = []
#    y_tom = []
#    z_tom = []
#
#
#    pi_tomo=p_ef
#    
#    
#    # plus direction
#    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=0*seq.num_patterns,
#                                       amp_tomo=amp_tom,pi_time=0,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
#    
#    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=1*seq.num_patterns,
#                                       amp_tomo=amp_tom,pi_time=0,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
#    
#    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=2*seq.num_patterns,
#                                       amp_tomo=amp_tom,pi_time=0,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
#    
#    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=3*seq.num_patterns,
#                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
#    
#    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=4*seq.num_patterns,
#                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
#    
#    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=5*seq.num_patterns,
#                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=ph_ang,coef_in=1)
#    
#    
#    #    # minus direction
#    #    freq_detun = -freq_detun
#    #    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#    #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=3*seq.num_patterns,
#    #                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    #    
#    #    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#    #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=4*seq.num_patterns,
#    #                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ef/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    #
#    #    black_nonHermitian_ma_V1.encirclingNHHEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#    #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=5*seq.num_patterns,
#    #                                       amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#    #    pi_tomo=0
#    #    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#    #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,
#    #                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
#    #
#    #    pi_tomo=p_ef
#    #    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#    #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=seq.num_patterns,
#    #                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
#       
#    black_nonHermitian_ma_V1.loading(num_steps = 2*3*seq.num_patterns)
#        
#
#    
#    
#    for idx, _ in enumerate(sweep.vals):
#        
#        print(idx)
#        threshold = [-166.05195052355808, -141.55195052355808];
#    
#        
#        
##        [-0.076, 0.025 , 0.09, -0.086]
#        wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5], offset=[-0.076, 0.025 , 0.09, -0.086])  # .033,.01, 0.033, -0.077; offset=[-0.062, 0.016 , 0.09, -0.086]
#        
#        #    threshold = [-165.44885667419695, -140.69885667419695];    
#        
#        for k in range(seq.num_avgs):
#            
#            daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*3*seq.num_patterns,
#                                                                                       num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#            
#        
#            if k is 0:
#                p_readout = p
#        #                   
#            else:
#                p_readout += p
#        
#        p_readout = p_readout/seq.num_avgs
#        #    p_post = analysis.p_readout_postselected(p_readout)
#        p_post=analysis.p_readout_postselected_pief(p_readout)
#        
#        
#        x_pief = p_post[2][0:N]-p_post[1][0:N] # with pi_ef pusle at the end
#        y_pief = p_post[1][N:2*N]-p_post[2][N:2*N]
#        z_pief = p_post[1][2*N:3*N]-p_post[2][2*N:3*N]
#        
#        
#        pop_3state.append(p_readout)
#        pop_post.append(p_post)
#        x_tom.append(x_pief)
#        y_tom.append(y_pief)
#        z_tom.append(z_pief)    
#
#    
##    x_pief_m = 2*(p_readout[2][3*N:4*N]+p_readout[1][3*N:4*N])-1 # with pi_ef pusle at the end
##    y_pief_m = 1-2*(p_readout[2][4*N:5*N]+p_readout[1][4*N:5*N])
##    z_pief_m = 1-2*(p_readout[2][5*N:6*N]+p_readout[1][5*N:6*N])
#    
#    
#        plt.semilogy(seq.times,1-p_readout[0,0:N]) # submanifold population
#        plt.xlabel("time (us)")
#        plt.ylabel("pop_sub")
#        plt.show()
#
#
#        plt.plot(seq.times,x_pief,'b-') # plus direction (+)
#        plt.plot(seq.times,y_pief,'r-')
#        plt.plot(seq.times,z_pief,'k-')
#        
#    #    plt.plot(seq.times,x_pief_m,'b--') # minus  direction (-)
#    #    plt.plot(seq.times,y_pief_m,'r--')
#    #    plt.plot(seq.times,z_pief_m,'k--')
#        plt.xlabel("time (us)")
#        plt.ylabel("x,y,z")
#    #    plt.legend(["x(+)","y(+)","z(+)"])
#        plt.ylim([-1,1])
#        plt.show()
#        
#        # analysis: entropy (purity), chirality v.s. duration time T
#        R_p = np.sqrt(x_pief**2 + y_pief**2 + z_pief**2)
#    #    R_m = np.sqrt(x_pief_m**2 + y_pief_m**2 + z_pief_m**2)
#        entropy_p = -(1+R_p)/2*np.log2((1+R_p)/2)-(1-R_p)/2*np.log2((1-R_p)/2)  
#    #    entropy_m = -(1+R_m)/2*np.log2((1+R_m)/2)-(1-R_m)/2*np.log2((1-R_m)/2)  
#        
#    #    chirality = 1/2*np.sqrt((x_pief-x_pief_m)**2 + (y_pief-y_pief_m)**2 + (z_pief-z_pief_m)**2) # defined by trace distance
#        
#        plt.plot(seq.times,entropy_p,'b-')
#    #    plt.plot(seq.times,entropy_m,'r-')
#    #    plt.plot(seq.times,chirality,'k--')
#        plt.xlabel("duration time (us)")
#        plt.ylabel("entropy")
#    #    plt.legend(["entropy(+)","entropy(-)","chirality"])
#        plt.ylim([-0.1,1])
#        plt.show()
#
#
#
#    pop_3state_stack=np.dstack(pop_3state) # dimension: 3 (g,e,f) x number of patterns x sweep steps
#    pop_post_stack = np.dstack(pop_post)
#    x_tom_stack = np.stack(x_tom)
#    y_tom_stack = np.stack(y_tom)
#    z_tom_stack = np.stack(z_tom)
#    
##    pop_post1_stack = np.dstack(pop_post1)
##    msmt_fit_stack = np.stack(msmt_fit)
#    pop_3state_stack_reshape = pop_3state_stack.reshape(seq.num_patterns*6*3,sweep.steps)
#    pop_post_stack_reshape = pop_post_stack.reshape(int(seq.num_patterns*6*3/2),sweep.steps)
##    pop3state_f=popdstack.reshape(seq.num_patterns*3,steps)
##    np.savetxt('pop_3state_stack', pop_3state_stack)
##    np.savetxt('pop_post_stack', pop_post_stack)
#    
#    fig = plt.figure(figsize=(6,4))
#    plt.imshow(np.transpose(pop_post_stack[1,:,:]-pop_post_stack[2,:,:]),extent=[min(seq.times),max(seq.times),max(sweep.vals),min(sweep.vals)],cmap='bwr',aspect='auto')
##    plt.imshow(np.transpose(pop_post_stack[1,:,:]-pop_post_stack[2,:,:]),cmap='bwr',aspect='auto')
#    plt.xlabel("time")
#    plt.ylabel("ssm_ef")
#    plt.colorbar()
#    plt.show()
#    
#    
#    x_tom_mean = np.mean(x_tom_stack, axis=0)
#    y_tom_mean = np.mean(y_tom_stack, axis=0)
#    z_tom_mean = np.mean(z_tom_stack, axis=0)
#
#    x_tom_std = np.std(x_tom_stack, axis=0)
#    y_tom_std = np.std(y_tom_stack, axis=0)
#    z_tom_std = np.std(z_tom_stack, axis=0)
#
#    plt.plot(seq.times, z_tom_mean,'k-')
#    plt.plot(seq.times, x_tom_mean,'b-')
#    plt.plot(seq.times, y_tom_mean,'r-')
#    plt.plot(seq.times, z_tom_mean + z_tom_std,'k--')
#    plt.plot(seq.times, x_tom_mean + x_tom_std,'b--')
#    plt.plot(seq.times, y_tom_mean + y_tom_std,'r--')
#    plt.plot(seq.times, z_tom_mean - z_tom_std,'k--')
#    plt.plot(seq.times, x_tom_mean - x_tom_std,'b--')
#    plt.plot(seq.times, y_tom_mean - y_tom_std,'r--')
#    plt.xlabel('time (us)')
#    plt.ylabel('x, y, z')
#    plt.ylim([-1,1])
#    plt.show()
#
#    plt.plot(seq.times, z_tom_std,'k-')
#    plt.plot(seq.times, x_tom_std,'b-')
#    plt.plot(seq.times, y_tom_std,'r-')
#    plt.xlabel('time (us)')
#    plt.ylabel('std of x,y,z')
#    plt.show()
#
#
#
#
##    plt.plot(sigma[0:101])
##    plt.plot(sigma[101:202])
##    plt.plot(sigma[202:303])
##    plt.shows()
#    
##    for ii in range(7):
##        msmt.popt, msmt.perr, _, _ = analysis.fit_sine_decay(seq.times, pop_post_stack[2,:,ii], guess_vals=[0.35,3,0.12,30,0.75])
##        msmt_fit.append(msmt.popt)
#
##    plt.plot(sweep.vals,abs(msmt_fit_stack[:,0]),'b-s')
##    plt.plot(sweep.vals,msmt_fit_stack[:,1],'r-o')
##    plt.xlabel("ssm_ef")
##    plt.ylabel("frequency (MHz), decay rate (us-1)")
##    plt.show()
#
###    
##    fig = plt.figure()
##    plt.plot(sweep.vals,abs(msmt_fit_stack[:,0]),'ro')
###    plt.plot(sweep.vals,3.629*sweep.vals-0.003)
##    analysis.fit_parabola(sweep.vals, abs(msmt_fit_stack[:,0]), guess_vals=None)
###    plt.xlabel("amp")
##    plt.xlabel("amp")
##    plt.ylabel("Frequency (MHz)")
##    plt.show
#    
#    np.savetxt(save_dir+'pop3state_slow_driving_tomography_4us_J30_detun15MHz_minus direction_ssm_ef -0.1202_10000repetitions_10times_amp0.28_08202021',pop_3state_stack_reshape)
#    np.savetxt(save_dir+'pop_post_slow_driving_tomography_4us_J30_detun15MHz_minus direction_ssm_ef -0.1202_10000repetitions_10times_amp_0.28_08202021',pop_post_stack_reshape)
########    np.savetxt(save_dir+'pop_post fit_pop f_2us_WX 1.5_amp 0.05_ssm_ef 0.0895 0.0935_07252020',msmt_fit_stack)   
######           
#    np.savetxt(save_dir+'slow_driving_tomographyX_10000_repetitions_10times_ssm_ef -0.1202_amp0.28_08202021',x_tom_stack)   
#    np.savetxt(save_dir+'slow_driving_tomographyY_10000_repetitions_10times_ssm_ef -0.1202_amp0.28_08202021',y_tom_stack)   
#    np.savetxt(save_dir+'slow_driving_tomographyZ_10000_repetitions_10times_ssm_ef -0.1202_amp0.28_08202021',z_tom_stack)   
##

    
    ###########################################################

    
## define the encircling loop (encircling type-I LEP) Purity or chirality (distance of final states of two directions) vs Jmin
#    amp_min = 0 # 0.005 #0.01
#    amp_max = 0.19 #0.12 #0.22 #0.25
##    freq_detun = 0.015
#
#    sweep.comment = "sweep encircling duration"
#    sweep.steps = 201 #201
##    sweep.vals = np.linspace(1,3001, sweep.steps) # sweep encircling duration  
#    sweep.vals = np.linspace(-1,1, sweep.steps) # sweep amp_min
##    sweep.vals = np.linspace(0, 0.02, sweep.steps) # sweep frequency detuning
#
#    seq.comment = "measurement"
#    seq.num_patterns = sweep.steps
##    seq.times = sweep.vals*1e-3
#    seq.times = sweep.vals
#
#    seq.sweep_time = 2000 #2000
#    seq.num_records_per_pattern = 500
##    seq.times = np.linspace(0.0, seq.sweep_time, seq.num_patterns)*1e-3
#    seq.num_avgs= 20
#    
##    amp_tom = 0.525
##    amp_tom = 0.5
##    amp_tom = 0.49 #0.53 #0.522
#    amp_tom = 0.51 #0.52
#    
#    pi_tomo=p_ef
#    
#    for idx, amp_min in enumerate(sweep.vals): # sweep encirlcing duration
#
##    for idx, temp_freq in enumerate(sweep.vals): # sweeep amp_min
#        print(idx)
#        
#        # plus direction
##        freq_detun = temp_freq #0.005
#        freq_detun = 0.005
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=90,off_set=0+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=0,off_set=1+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=0,off_set=2+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
#        # minus direction
#        freq_detun = -freq_detun
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=90,off_set=3+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=180+0,off_set=4+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=0,off_set=5+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
##    pi_tomo=0
##    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,
##                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
##
##    pi_tomo=p_ef
##    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=seq.num_patterns,
##                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
#   
#    black_nonHermitian_ma_V1.loading(num_steps = 2*3*seq.num_patterns) #.033,.01, 0.033, -0.077
#        
#        
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])  # offset=[-0.062, 0.016 , 0.09, -0.086] -0.096, 0.028 , 0.09, -0.086
#    threshold = [-147.0968567120443, -134.5968567120443];    
#    
#    for k in range(seq.num_avgs):
#        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*3*seq.num_patterns,
#                                                                                   num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        
#
#        if k is 0:
#            p_readout = p
##                   
#        else:
#            p_readout += p
#
#    p_readout = p_readout/seq.num_avgs
##    p_post = analysis.p_readout_postselected(p_readout)
##    p_post=analysis.p_readout_postselected_pief(p_readout)
#
#    N = seq.num_patterns    
#    x_pief = 2*(p_readout[2][0::6]+p_readout[1][0::6])-1 # with pi_ef pusle at the end
#    y_pief = 1-2*(p_readout[2][1::6]+p_readout[1][1::6])
#    z_pief = 1-2*(p_readout[2][2::6]+p_readout[1][2::6])
#    
#    x_pief_m = 2*(p_readout[2][3::6]+p_readout[1][3::6])-1 # with pi_ef pusle at the end
#    y_pief_m = 2*(p_readout[2][4::6]+p_readout[1][4::6])-1
#    z_pief_m = 1-2*(p_readout[2][5::6]+p_readout[1][5::6])
#    
#
#
#    plt.plot(seq.times,x_pief,'b-')
#    plt.plot(seq.times,y_pief,'r-')
#    plt.plot(seq.times,z_pief,'k-')
#    
#    plt.plot(seq.times,x_pief_m,'b--')
#    plt.plot(seq.times,y_pief_m,'r--')
#    plt.plot(seq.times,z_pief_m,'k--')
##    plt.xlabel("duration time (us)")
#    plt.xlabel("freq detun (GHz)")
#    plt.ylabel("x,y,z")
#    plt.legend(["x(+)","y(+)","z(+)","x(-)","y(-)","z(-)"])
#    plt.ylim([-1,1])
#    plt.show()
#    
#    
#    # analysis: entropy (purity), chirality v.s. duration time T
#    R_p = np.sqrt(x_pief**2 + y_pief**2 + z_pief**2)
#    R_m = np.sqrt(x_pief_m**2 + y_pief_m**2 + z_pief_m**2)
#    entropy_p = -(1+R_p)/2*np.log2((1+R_p)/2)-(1-R_p)/2*np.log2((1-R_p)/2)  
#    entropy_m = -(1+R_m)/2*np.log2((1+R_m)/2)-(1-R_m)/2*np.log2((1-R_m)/2)  
#    
#    chirality = 1/2*np.sqrt((x_pief-x_pief_m)**2 + (y_pief-y_pief_m)**2 + (z_pief-z_pief_m)**2) # defined by trace distance
#    
#    plt.plot(seq.times,entropy_p,'b-')
#    plt.plot(seq.times,entropy_m,'r-')
#    plt.plot(seq.times,chirality,'k--')
##    plt.xlabel("duration time (us)")
#    plt.xlabel("freq detun (GHz)")
#    plt.ylabel("entropy, chirality")
#    plt.legend(["entropy(+)","entropy(-)","chirality"])
#    plt.ylim([-0.1,1])
#    plt.show()
#
#    np.savetxt(save_dir+'tomography_encirclingLEPI_sweepampmin-1to1_T2us_ampmax0.19freq5MHzTwoDirections_PlusX_NoteOppositeYforMinusDirection_' + add_t, p_readout)
##    np.savetxt(save_dir+'tomography_encirclingLEPI_sweepFreqDetun0to20MHz_T2us_ampmax0.12ampmin0.005TwoDirections_MinusX_NoteOppositeYforMinusDirection_second try_09162021', p_readout)
####    
##    


#   #########################################################
    
## define the encircling loop (encircling type-I LEP) Purity or chirality (distance of final states of two directions) vs duration T (or Jmin)
#    amp_min = 0.005 #0.01
#    amp_max = 0.12 #0.22 #0.25
##    freq_detun = 0.015
#
#    sweep.comment = "sweep encircling duration"
#    sweep.steps = 201
#    sweep.vals = np.linspace(1,3001, sweep.steps) # sweep encircling duration  
##    sweep.vals = np.linspace(-0.2,0.2, sweep.steps) # sweep amp_min
##    sweep.vals = np.linspace(0, 0.02, sweep.steps) # sweep frequency detuning
#
#    seq.comment = "measurement"
#    seq.num_patterns = sweep.steps
##    seq.times = sweep.vals*1e-3
#    seq.times = sweep.vals
#
#    seq.sweep_time = 2000 #2000
#    seq.num_records_per_pattern = 1000
##    seq.times = np.linspace(0.0, seq.sweep_time, seq.num_patterns)*1e-3
#    seq.num_avgs= 10
#    
##    amp_tom = 0.525
##    amp_tom = 0.5
##    amp_tom = 0.49 #0.53 #0.522
#    amp_tom = 0.5 #0.52
#    
#    pi_tomo=p_ef
#    
#    for idx, seq.sweep_time in enumerate(sweep.vals): # sweep encirlcing duration
#
##    for idx, temp_freq in enumerate(sweep.vals): # sweeep amp_min
#        print(idx)
#        
#        # plus direction
##        freq_detun = temp_freq #0.005
#        freq_detun = 0.005
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=90,off_set=0+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=0,off_set=1+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=0,off_set=2+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
#        # minus direction
#        freq_detun = -freq_detun
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=90,off_set=3+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=180+0,off_set=4+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
#        black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                ,ssm_ef = ssb_ef,phase_tomo=0,off_set=5+idx*6,
#                                                amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#
##    pi_tomo=0
##    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,
##                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
##
##    pi_tomo=p_ef
##    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=seq.num_patterns,
##                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
#   
#    black_nonHermitian_ma_V1.loading(num_steps = 2*3*seq.num_patterns) #.033,.01, 0.033, -0.077
#        
#        
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])  # offset=[-0.062, 0.016 , 0.09, -0.086] -0.096, 0.028 , 0.09, -0.086
#    threshold = [-171.9742140308683, -148];    
#    
#    for k in range(seq.num_avgs):
#        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*3*seq.num_patterns,
#                                                                                   num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        
#
#        if k is 0:
#            p_readout = p
##                   
#        else:
#            p_readout += p
#
#    p_readout = p_readout/seq.num_avgs
##    p_post = analysis.p_readout_postselected(p_readout)
##    p_post=analysis.p_readout_postselected_pief(p_readout)
#
#    N = seq.num_patterns    
#    x_pief = 2*(p_readout[2][0::6]+p_readout[1][0::6])-1 # with pi_ef pusle at the end
#    y_pief = 1-2*(p_readout[2][1::6]+p_readout[1][1::6])
#    z_pief = 1-2*(p_readout[2][2::6]+p_readout[1][2::6])
#    
#    x_pief_m = 2*(p_readout[2][3::6]+p_readout[1][3::6])-1 # with pi_ef pusle at the end
#    y_pief_m = 2*(p_readout[2][4::6]+p_readout[1][4::6])-1
#    z_pief_m = 1-2*(p_readout[2][5::6]+p_readout[1][5::6])
#    
#
#
#    plt.plot(seq.times,x_pief,'b-')
#    plt.plot(seq.times,y_pief,'r-')
#    plt.plot(seq.times,z_pief,'k-')
#    
#    plt.plot(seq.times,x_pief_m,'b--')
#    plt.plot(seq.times,y_pief_m,'r--')
#    plt.plot(seq.times,z_pief_m,'k--')
##    plt.xlabel("duration time (us)")
#    plt.xlabel("freq detun (GHz)")
#    plt.ylabel("x,y,z")
#    plt.legend(["x(+)","y(+)","z(+)","x(-)","y(-)","z(-)"])
#    plt.ylim([-1,1])
#    plt.show()
#    
#    
#    # analysis: entropy (purity), chirality v.s. duration time T
#    R_p = np.sqrt(x_pief**2 + y_pief**2 + z_pief**2)
#    R_m = np.sqrt(x_pief_m**2 + y_pief_m**2 + z_pief_m**2)
#    entropy_p = -(1+R_p)/2*np.log2((1+R_p)/2)-(1-R_p)/2*np.log2((1-R_p)/2)  
#    entropy_m = -(1+R_m)/2*np.log2((1+R_m)/2)-(1-R_m)/2*np.log2((1-R_m)/2)  
#    
#    chirality = 1/2*np.sqrt((x_pief-x_pief_m)**2 + (y_pief-y_pief_m)**2 + (z_pief-z_pief_m)**2) # defined by trace distance
#    
#    plt.plot(seq.times,entropy_p,'b-')
#    plt.plot(seq.times,entropy_m,'r-')
#    plt.plot(seq.times,chirality,'k--')
##    plt.xlabel("duration time (us)")
#    plt.xlabel("freq detun (GHz)")
#    plt.ylabel("entropy, chirality")
#    plt.legend(["entropy(+)","entropy(-)","chirality"])
#    plt.ylim([-0.1,1])
#    plt.show()
#
##    np.savetxt(save_dir+'tomography_encirclingLEPI_sweepT0usto3us_ampmin0.005ampmax0.12freq5MHzTwoDirections_MinusX_NoteOppositeYforMinusDirection_09162021', p_readout)
##    np.savetxt(save_dir+'tomography_encirclingLEPI_sweepFreqDetun0to20MHz_T2us_ampmax0.12ampmin0.005TwoDirections_MinusX_NoteOppositeYforMinusDirection_second try_09162021', p_readout)
####    
##    
    ###########################################################
    
    
## define the encircling loop (encircling type-I LEP) Purity or chirality (distance of final states of two directions) vs duration T (or Jmin)
## 2D sweep (sweep duration T and Jmin)
#
#
#    amp_min = 0.005 #0.01
#    amp_max = 0.22 # 0.2 #0.25
##    freq_detun = 0.015
#
#    sweep.comment = "sweep encircling duration"
#    sweep.steps = 51 #201
#    sweep.steps1 = 60
#    sweep.vals = np.linspace(-0.22,0.22, sweep.steps) # sweep amp_min
#    sweep.vals1 = np.linspace(50,3000, sweep.steps1) # sweep encircling duration  
#
#    seq.comment = "measurement"
#    seq.num_patterns = sweep.steps
##    seq.times = sweep.vals*1e-3
#    seq.times = sweep.vals
#
#    seq.sweep_time = 2000 #2000
#    seq.num_records_per_pattern = 1000
##    seq.times = np.linspace(0.0, seq.sweep_time, seq.num_patterns)*1e-3
#    seq.num_avgs=10
#    
#    amp_tom = 0.525 #0.522
##    amp_tom = 0.5
#
#    pi_tomo=p_ef
#    
#    
#    pop_3state = []
#    x_tom_p = []
#    y_tom_p = []
#    z_tom_p = []
#    x_tom_m = []
#    y_tom_m = []
#    z_tom_m = []
#    
#    for idx, seq.sweep_time in enumerate(sweep.vals1): # sweep encirlcing duration
#        
#        print(idx/sweep.steps1)
#        
#        for idx, amp_min in enumerate(sweep.vals): # sweeep amp_min
#            print(idx)
#            
#            # plus direction
#            freq_detun = 0.005
#            black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                    ,ssm_ef = ssb_ef,phase_tomo=90,off_set=0+idx*6,
#                                                    amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#        
#            black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                    ,ssm_ef = ssb_ef,phase_tomo=0,off_set=1+idx*6,
#                                                    amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#            black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                    ,ssm_ef = ssb_ef,phase_tomo=0,off_set=2+idx*6,
#                                                    amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#            # minus direction
#            freq_detun = -freq_detun
#            black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                    ,ssm_ef = ssb_ef,phase_tomo=90,off_set=3+idx*6,
#                                                    amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#        
#            black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                    ,ssm_ef = ssb_ef,phase_tomo=0,off_set=4+idx*6,
#                                                    amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=p_ge/2,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#            black_nonHermitian_ma_V1.encirclingLEP1_tVar(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                    ,ssm_ef = ssb_ef,phase_tomo=0,off_set=5+idx*6,
#                                                    amp_tomo=amp_tom,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#    #    pi_tomo=0
#    #    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#    #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,
#    #                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
#    #
#    #    pi_tomo=p_ef
#    #    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#    #                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=seq.num_patterns,
#    #                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
#       
#        black_nonHermitian_ma_V1.loading(num_steps = 2*3*seq.num_patterns)
#        wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[-0.076, 0.025 , 0.09, -0.086])  # offset=[-0.062, 0.016 , 0.09, -0.086]
#        threshold = [-169.06348890613324, -147.56348890613324];    
#        
#        for k in range(seq.num_avgs):
#            
#            daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*3*seq.num_patterns,
#                                                                                       num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#            
#    
#            if k is 0:
#                p_readout = p
#    #                   
#            else:
#                p_readout += p
#    
#        p_readout = p_readout/seq.num_avgs
#    #    p_post = analysis.p_readout_postselected(p_readout)
#    #    p_post=analysis.p_readout_postselected_pief(p_readout)
#    
#        N = seq.num_patterns    
#        x_pief = 2*(p_readout[2][0::6]+p_readout[1][0::6])-1 # with pi_ef pusle at the end
#        y_pief = 1-2*(p_readout[2][1::6]+p_readout[1][1::6])
#        z_pief = 1-2*(p_readout[2][2::6]+p_readout[1][2::6])
#        
#        x_pief_m = 2*(p_readout[2][3::6]+p_readout[1][3::6])-1 # with pi_ef pusle at the end
#        y_pief_m = 1-2*(p_readout[2][4::6]+p_readout[1][4::6])
#        z_pief_m = 1-2*(p_readout[2][5::6]+p_readout[1][5::6])
#        
#        pop_3state.append(p_readout)
#        x_tom_p.append(x_pief)  # p: pluse direction
#        y_tom_p.append(y_pief)
#        z_tom_p.append(z_pief)
#        x_tom_m.append(x_pief_m) # m: minus direction 
#        y_tom_m.append(y_pief_m)
#        z_tom_m.append(z_pief_m)
#    
#        plt.plot(seq.times,x_pief,'b-')
#        plt.plot(seq.times,y_pief,'r-')
#        plt.plot(seq.times,z_pief,'k-')
#        
#        plt.plot(seq.times,x_pief_m,'b--')
#        plt.plot(seq.times,y_pief_m,'r--')
#        plt.plot(seq.times,z_pief_m,'k--')
#    #    plt.xlabel("duration time (us)")
#        plt.xlabel("amp_min")
#        plt.ylabel("x,y,z")
#        plt.legend(["x(+)","y(+)","z(+)","x(-)","y(-)","z(-)"])
#        plt.ylim([-1,1])
#        plt.show()
#        
#        
#        # analysis: entropy (purity), chirality v.s. duration time T
#        R_p = np.sqrt(x_pief**2 + y_pief**2 + z_pief**2)
#        R_m = np.sqrt(x_pief_m**2 + y_pief_m**2 + z_pief_m**2)
#        entropy_p = np.real(-(1+R_p)/2*np.log2((1+R_p)/2)-(1-R_p)/2*np.log2((1-R_p)/2))  
#        entropy_m = np.real(-(1+R_m)/2*np.log2((1+R_m)/2)-(1-R_m)/2*np.log2((1-R_m)/2))  
#        
#        chirality = 1/2*np.sqrt((x_pief-x_pief_m)**2 + (y_pief-y_pief_m)**2 + (z_pief-z_pief_m)**2) # defined by trace distance
#        
#        plt.plot(seq.times,entropy_p,'b-')
#        plt.plot(seq.times,entropy_m,'r-')
#        plt.plot(seq.times,chirality,'k--')
#    #    plt.xlabel("duration time (us)")
#        plt.xlabel("amp_min")
#        plt.ylabel("entropy, chirality")
#        plt.legend(["entropy(+)","entropy(-)","chirality"])
#        plt.ylim([-0.1,1])
#        plt.show()
#
#
#    pop_3state_stack=np.dstack(pop_3state) # dimension: 3 (g,e,f) x number of patterns x sweep steps
#    pop_3state_stack_reshape = pop_3state_stack.reshape(sweep.steps1*3*2*3,sweep.steps)
#    
#    x_tom_p_stack = np.stack(x_tom_p)
#    y_tom_p_stack = np.stack(y_tom_p)
#    z_tom_p_stack = np.stack(z_tom_p)
#    x_tom_m_stack = np.stack(x_tom_m)
#    y_tom_m_stack = np.stack(y_tom_m)
#    z_tom_m_stack = np.stack(z_tom_m)
#    
#    
#    fig = plt.figure(figsize=(6,4))
#    plt.imshow(np.transpose(y_tom_p_stack[1:-1,:]),extent=[min(sweep.vals1),max(sweep.vals1),min(seq.times),max(seq.times)],cmap='bwr',aspect='auto')
##    plt.imshow(np.transpose(pop_post_stack[1,:,:]-pop_post_stack[2,:,:]),cmap='bwr',aspect='auto')
#    plt.xlabel("time (ns)")
#    plt.ylabel("Jmin")
#    plt.colorbar()
#    plt.show()
#
#    fig = plt.figure(figsize=(6,4))
#    plt.pcolor(np.transpose(y_tom_p_stack[1:-1,:]), vmin=-1, vmax=1, cmap='bwr')
#    plt.xlabel("time(ns)")
#    plt.ylabel("Jmin")
#    plt.colorbar()
#    plt.show()
#
#                        
#
##    np.savetxt(save_dir+'tomography_encirclingLEPI_sweepT0.5usto3us_sweep_ampmin-0.22to0.22_ampmax0.22freq5MHzTwoDirections_MinusX_10000repetitions_09072021', pop_3state_stack_reshape)
##    np.savetxt(save_dir+'tomography_X_sweep_Jmin_sweep_T_10000repetitions_plus_direction_09072021', x_tom_p_stack)
##    np.savetxt(save_dir+'tomography_Y_sweep_Jmin_sweep_T_10000repetitions_plus_direction_09072021', y_tom_p_stack)
##    np.savetxt(save_dir+'tomography_Z_sweep_Jmin_sweep_T_10000repetitions_plus_direction_09072021', z_tom_p_stack)
##    np.savetxt(save_dir+'tomography_X_sweep_Jmin_sweep_T_10000repetitions_minus_direction_09072021', x_tom_m_stack)
##    np.savetxt(save_dir+'tomography_Y_sweep_Jmin_sweep_T_10000repetitions_minus_direction_09072021', y_tom_m_stack)
##    np.savetxt(save_dir+'tomography_Z_sweep_Jmin_sweep_T_10000repetitions_minus_direction_09072021', z_tom_m_stack)    
#    
    ###########################################################
    
## define the encircling loop (encircling type-I LEP) Purity or chirality (distance of final states of two directions) vs duration T
#    amp_min = 0.01
#    amp_max = 0.25
#    freq_detun = -0.015
#
#    seq.comment = "measurement"
#    seq.num_patterns = 151#51
#    
##    0.3
#    seq.sweep_time = 500
#    seq.num_records_per_pattern = 500
#    seq.times = np.linspace(0.0, seq.sweep_time, seq.num_patterns)*1e-3
#    seq.num_avgs=2
#
#   
##    pi_tomo=0
##    black_nonHermitian_ma_V1.encirclingLEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,
##                                       amp_tomo=0.5,pi_time=pi_tomo,pi2_time=p_ge,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##    
##    black_nonHermitian_ma_V1.encirclingLEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=1*seq.num_patterns,
##                                       amp_tomo=0.5,pi_time=pi_tomo,pi2_time=p_ge,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##
##    black_nonHermitian_ma_V1.encirclingLEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=2*seq.num_patterns,
##                                       amp_tomo=0.5,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##    
##    black_nonHermitian_ma_V1.encirclingLEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=3*seq.num_patterns,
##                                       amp_tomo=0.5,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##
##
##    pi_tomo=p_ef
##    black_nonHermitian_ma_V1.encirclingLEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=4*seq.num_patterns,
##                                       amp_tomo=0.5,pi_time=pi_tomo,pi2_time=p_ge,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##    
##    black_nonHermitian_ma_V1.encirclingLEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=5*seq.num_patterns,
##                                       amp_tomo=0.5,pi_time=pi_tomo,pi2_time=p_ge,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##
##    black_nonHermitian_ma_V1.encirclingLEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=6*seq.num_patterns,
##                                       amp_tomo=0.5,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
##    
##    black_nonHermitian_ma_V1.encirclingLEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=90,off_set=7*seq.num_patterns,
##                                       amp_tomo=0.5,pi_time=pi_tomo,pi2_time=0,min_amp=amp_min,detun=freq_detun,amp_enc=amp_max,ph=0,coef_in=1)
#    
#    pi_tomo=0
#    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,
#                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
#
#    pi_tomo=p_ef
#    black_nonHermitian_ma_V1.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=seq.num_patterns,
#                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=10,coef_in=1)
#    
#    black_nonHermitian_ma_V1.loading(num_steps = 2*seq.num_patterns)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])  # offset=[-0.062, 0.016 , 0.09, -0.086]
#    threshold = [-160.96600532200222, -155.71600532200222];    
#    
#    for k in range(seq.num_avgs):
#        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=threshold,num_patterns=2*seq.num_patterns,
#                                                                                   num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        
#
#        if k is 0:
#            p_readout = p
##                   
#        else:
#            p_readout += p
#
##    p_post = analysis.p_readout_postselected(p_readout)
##    p_post=analysis.p_readout_postselected_pief(p_readout)
#
#    N = seq.num_patterns    
#    r_gf_real_pief = 1/2*(p_readout[2][4*N:5*N]-p_readout[2][0:N]) # with pi_ef pusle at the end
#    r_gf_imag_pief = 1/2*(p_readout[2][5*N:6*N]-p_readout[2][N:2*N])
#    r_ef_real_pief = 1/2*(p_readout[2][6*N:7*N]-p_readout[2][2*N:3*N])
#    r_ef_imag_pief = 1/2*(p_readout[2][7*N:8*N]-p_readout[2][3*N:4*N])
#    
#
#    x=seq.times
##    y=p_readout[1]
##    y=p_post[1]
#    plt.plot(x,r_gf_real_pief,'b-')
#    plt.plot(x,r_gf_imag_pief,'b--')
#    plt.plot(x,r_ef_real_pief,'r-')
#    plt.plot(x,r_ef_imag_pief,'r--')
#    
##    ratio = (np.square(r_gf_real_pief) + np.square(r_gf_imag_pief))/(np.square(r_ef_real_pief) + np.square(r_ef_imag_pief))
##    plt.plot(x,ratio)    
##    plt.ylim([0,15])
##    popt,peer, y_vals, _ = analysis.fit_sine_decay(x,y,guess_vals=None) 
##    popt,peer, y_vals, _ = analysis.fit_exp_decay(x,y,guess_vals=None)
##    popt,peer, y_vals, _,maxf = analysis.fit_sine(x,y,guess_vals=[14,.4,90,.5])  
##    popt, y_vals = analysis.fit_sine(x,y,guess_vals=None)
##    ind=np.argmax(y_vals)
##    print(x[ind])
#                 
##        p_readout = p_readout/seq.num_avgs
##        p_post = analysis.p_readout_postselected(p_readout)
##        plt.plot(p_post[1])
##        plt.show()
#####
##
#
##    np.savetxt(save_dir+'ztomo2usJ60to30radpusD5', p_readout)
##    np.savetxt(save_dir+'ztomo2usJ60to30radpusD5post', p_post)

    ###################################################################################


    
    
#    pi_tomo=0
#    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = rtime,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,amp_tomo=0.5
#                                                 ,pi_time=pi_tomo,min_amp=-1,max_amp=0,detun=0.005,amp_enc=.25,ph=0,pl=1)
    
#    black_nonHermitian_ma.ch12_leakage(pi_ge=p_ge,pi_ef=p_ef,ssm_ge = ssb_hf,ssm_ef = ssb_hf,x=90)    
#    black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,t1_time = seq.sweep_time,num_steps =seq.num_patterns,pi_ge=p_ge,pi_ef=p_ef,pi_echo=0)
#    black_nonHermitian_ma.ramsey(ssm_ge = ssb_ge,off_set=0,num_steps = seq.num_patterns,t1_time = seq.sweep_time,pi_ge=p_ge,amp=.5,pi_echo=0)     

#    black_nonHermitian_ma.t1_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,t1_time = seq.sweep_time,t_min=0,pi_ge=p_ge,pi_ef=p_ef,amp_g=.505)   
#    black_nonHermitian_ma.t1(ssm_ge = ssb_ge,off_set=0,num_steps=seq.num_patterns,t1_time = seq.sweep_time,t_min=0,pi_ge=p_ge,amp_g=.5)
#    black_nonHermitian_ma.t1_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,ssm_hf = ssb_hf,pi_ge=p_ge,pi_ef=p_ef,pi_hf=p_hf,amp=.5,t1_time = seq.sweep_time,num_steps = seq.num_patterns,off_set=0)
#
#    black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,
#                                  rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.5,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
#    black_nonHermitian_ma.rabi(ssm_ge = ssb_ge, num_steps = seq.num_patterns,amp = .5,rabi_time = seq.sweep_time)  
#    black_nonHermitian_ma.phase_meas_esdrive_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=180)
#    black_nonHermitian_ma.xdrive(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=188,ort=-5)
#    black_nonHermitian_ma.ramsey_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,t2_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,xp=0,p2ef=p_ef*0)
#    black_nonHermitian_ma.rabi_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.786,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf)
#    black_nonHermitian_ma.phase_meas_encircling(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
##                                                ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=0,detun=0.005,jmin=0.01)
#    black_nonHermitian_ma.check_phase_enc(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                    ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=180
###                                                    ,detun=-0.005,min_amp=0,amp_enc=.25,amp_h=.5)

    
####    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
###                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=180,off_set=seq.num_patterns,amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,max_amp=.5,detun=-0.005)
###    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
###                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=180,off_set=seq.num_patterns,amp_tomo=0.5
###                                                 ,pi_time=pi_tomo,min_amp=.01,max_amp=.5,detun=-0.003,amp_enc=.17,ph=15)
###    black_nonHermitian_ma.encirclingEPJvary(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
###                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=201,amp_tomo=0.5,pi_time=pi_tomo)
####
#    black_nonHermitian_ma.phase_meas_encircling(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                ,off_set=seq.num_patterns*idx,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=180
#                                                ,detun=0.005,min_amp=.01,amp_enc=.25,amp_h=.5)
#    black_nonHermitian_ma.check_phase_enc(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                    ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=180,esph=180,exph=15
#                                                    ,detun=0.005,min_amp=0,amp_enc=.25,amp_h=.5,coef_in=.95) 

#    black_nonHermitian_ma.test_phase_meas_encircling(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                    ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=180,esph=0,detun=-0.005,min_amp=.01)   
#    black_nonHermitian_ma.phase_meas_encircling_jmin(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                     ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                     ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=180,esph=0
##                                                     ,detun=-0.005,amp_enc=.25,amp_h=.5,min_amp=0,max_amp=1)
#    black_nonHermitian_ma.no_pi_pi_pipi(off_set=0,pi_ge=p_ge,pi_ef=p_ef,ssm_ge = ssb_ge,ssm_ef= ssb_ef,p1=0,p2=0)
#    black_nonHermitian_ma.no_pi_pi_pipi(off_set=1,pi_ge=p_ge,pi_ef=p_ef,ssm_ge = ssb_ge,ssm_ef= ssb_ef,p1=p_ge,p2=0)
#    black_nonHermitian_ma.no_pi_pi_pipi(off_set=2,pi_ge=p_ge,pi_ef=p_ef,ssm_ge = ssb_ge,ssm_ef= ssb_ef,p1=p_ge,p2=p_ef,p3=0,ssm_hf=ssb_ef)
##    black_nonHermitian_ma.no_pi_pi_pipi(off_set=2,pi_ge=p_ge,pi_ef=p_ef,ssm_ge = ssb_ge,ssm_ef= ssb_ef,p1=p_ge,p2=p_ef,p3=p_hf,ssm_hf=ssb_ef)


###########
####################################    
#################
#################
#################
#################
################################
#################################    daq_params, t_histo, p_readout,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
#################################            num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0)
################################
#    for k in range(seq.num_avgs):
#        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=m,num_patterns=2*seq.num_patterns,
#                                                                                   num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        
#
#        if k is 0:
#            p_readout = p
##                   
#        else:
#            p_readout += p
#
#    p_post = analysis.p_readout_postselected(p_readout)
#    p_post=analysis.p_readout_postselected_pief(p_readout)
#    
#    r_gf_real_pief = 1/2*(p_readout[2][N:2*N]-p_readout[2][0:N]) # with pi_ef pusle at the end
#
#    x=seq.times
#    y=p_readout[1]
#    y=p_post[1]     
#    plt.plot(x,r_gf_real_pief)
#    plt.ylim([-1,1])
#    popt,peer, y_vals, _ = analysis.fit_sine_decay(x,y,guess_vals=None) 
#    popt,peer, y_vals, _ = analysis.fit_exp_decay(x,y,guess_vals=None)
#    popt,peer, y_vals, _,maxf = analysis.fit_sine(x,y,guess_vals=[14,.4,90,.5])  
#    popt, y_vals = analysis.fit_sine(x,y,guess_vals=None)
#    ind=np.argmax(y_vals)
#    print(x[ind])
                 
#        p_readout = p_readout/seq.num_avgs
#        p_post = analysis.p_readout_postselected(p_readout)
#        plt.plot(p_post[1])
#        plt.show()
####
#

#    np.savetxt(save_dir+'ztomo2usJ60to30radpusD5', p_readout)
#    np.savetxt(save_dir+'ztomo2usJ60to30radpusD5post', p_post)
##    
    
#    black_nonHermitian_ma.t1(ssm_ge = ssb_ge,off_set=0,num_steps=seq.num_patterns,t1_time = 1000,t_min=0,amp_g=.505)
#    black_nonHermitian_ma.t1_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns1,t1_time = seq.sweep_time,t_min=0,pi_ge=p_ge,pi_ef=p_ef,amp_g=.505)
#    black_nonHermitian_ma.ramsey(ssm_ge = 0.3885,off_set=0)
#    black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,t1_time = 6000,num_steps =seq.num_patterns)
#    black_nonHermitian_ma.loading(num_steps = seq.num_patterns)
#    black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = 100)
#    black_nonHermitian_ma.ramsey(ssm_ge = ssb_ge,off_set=0,num_steps = seq.num_patterns,t1_time = 3000,pi_ge=p_ge)
#################################################################
#    black_nonHermitian_ma.check_phase_enc(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                    ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=180,esph=180
#                                                    ,detun=-0.005,min_amp=0,amp_enc=.25,amp_h=.5)    
#    black_nonHermitian_ma.loading(num_steps = 1*seq.num_patterns)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
#    prev_threshold = [152.385192421315, 155.385192421315];
#    for k in range(seq.num_avgs):
#        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
#            num_patterns=1*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        
#
#        if k is 0:
#            p_readout = p
##                   
#        else:
#            p_readout += p
##
#    p_post = analysis.p_readout_postselected(p_readout)
##    p_post=analysis.p_readout_postselected_pief(p_readout)
#    x=seq.times
#    y=p_readout[1]
#    y1=p_post[2]      
#    plt.plot(x,y1)
##################################################################
#    black_nonHermitian_ma.check_phase_enc(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                    ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=-180
#                                                    ,detun=-0.005,min_amp=0,amp_enc=.25,amp_h=.5)    
#    black_nonHermitian_ma.loading(num_steps = 1*seq.num_patterns)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
#    prev_threshold = [152.385192421315, 155.385192421315];
#    for k in range(seq.num_avgs):
#        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
#            num_patterns=1*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        
#
#        if k is 0:
#            p_readout = p
##                   
#        else:
#            p_readout += p
##
#    p_post = analysis.p_readout_postselected(p_readout)
##    p_post=analysis.p_readout_postselected_pief(p_readout)
#    x=seq.times
#    y=p_readout[1]
#    y2=p_post[2]      
#    plt.plot(x,y2)    
##################################################################
#    black_nonHermitian_ma.check_phase_enc(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                    ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=180
#                                                    ,detun=0.005,min_amp=0,amp_enc=.25,amp_h=.5)    
#    black_nonHermitian_ma.loading(num_steps = 1*seq.num_patterns)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
#    prev_threshold = [152.385192421315, 155.385192421315];
#    for k in range(seq.num_avgs):
#        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
#            num_patterns=1*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        
#
#        if k is 0:
#            p_readout = p
##                   
#        else:
#            p_readout += p
##
#    p_post = analysis.p_readout_postselected(p_readout)
##    p_post=analysis.p_readout_postselected_pief(p_readout)
#    x=seq.times
#    y=p_readout[1]
#    y3=p_post[2]      
#    plt.plot(x,y3)
##################################################################
#    black_nonHermitian_ma.check_phase_enc(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                    ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=180
#                                                    ,detun=0.005,min_amp=0,amp_enc=.25,amp_h=.5)    
#    black_nonHermitian_ma.loading(num_steps = 1*seq.num_patterns)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
#    prev_threshold = [152.385192421315, 155.385192421315];
#    for k in range(seq.num_avgs):
#        
#        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
#            num_patterns=1*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        
#
#        if k is 0:
#            p_readout = p
##                   
#        else:
#            p_readout += p
##
#    p_post = analysis.p_readout_postselected(p_readout)
##    p_post=analysis.p_readout_postselected_pief(p_readout)
#    x=seq.times
#    y=p_readout[1]
#    y4=p_post[2]      
#    plt.plot(x,y4)
##################################################################
#    st=21
##    sweep.vals=np.linspace(.407,.413,st)
##    sweep.vals=np.linspace(.108,.114,st)
#    sweep.vals=np.linspace(.218,.224,st)
#
#    for ink,ss in enumerate(sweep.vals):
#            print(ink)
##            black_nonHermitian_ma.ramsey(ssm_ge = ss,off_set=seq.num_patterns*ink,num_steps = seq.num_patterns,t1_time = seq.sweep_time,pi_ge=p_ge,amp=.5) 
##            black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ss,off_set=seq.num_patterns*ink,t1_time = seq.sweep_time,num_steps =seq.num_patterns,pi_ge=p_ge,pi_ef=p_ef)
#            black_nonHermitian_ma.ramsey_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=seq.num_patterns*ink,num_steps = seq.num_patterns,amp=.5,t2_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ss,pi_hf=p_hf,xp=0)
#
#    black_nonHermitian_ma.loading(num_steps = st*seq.num_patterns)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
##################################################################
#    black_nonHermitian_ma.t1(ssm_ge = ssb_ge,off_set=0,num_steps=seq.num_patterns,t1_time = 2000,t_min=0,pi_ge=p_ge,amp_g=.5)
#    black_nonHermitian_ma.t1_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=seq.num_patterns,num_steps = seq.num_patterns,t1_time = 6000,t_min=0,pi_ge=p_ge,pi_ef=p_ef,amp_g=.505)   
#    black_nonHermitian_ma.t1_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,ssm_hf = ssb_hf,pi_ge=p_ge,pi_ef=p_ef,pi_hf=p_hf,amp=.5,t1_time = 6000,num_steps = seq.num_patterns,off_set=2*seq.num_patterns)
#    black_nonHermitian_ma.loading(num_steps = 3*seq.num_patterns)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
###################################################################
#    st=21
#    sweep.vals=np.linspace(.48,.52,st)
#    sweep.vals=np.linspace(.78,.8,st)
#
#    for ink,ss in enumerate(sweep.vals):
#            print(ink)
##            black_nonHermitian_ma.rabi(ssm_ge = ssb_ge, num_steps = seq.num_patterns,amp = ss,rabi_time = seq.sweep_time,off_set=seq.num_patterns*ink)  
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=seq.num_patterns*ink,num_steps = seq.num_patterns,amp=ss,
##                                          rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.504,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
#            black_nonHermitian_ma.rabi_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=seq.num_patterns*ink,num_steps = seq.num_patterns,amp=ss,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf)
####
#    black_nonHermitian_ma.loading(num_steps = st*seq.num_patterns)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
#################################################################
#    sweep.vals=np.linspace(7.34,7.35,21)
#    sweep.vals=np.linspace(.165,.18,11)
####
###
####    #    sweep.vals=np.linspace(.219,.221,21)
######    sweep.vals=np.linspace(.76,.78,11 )
#####    sweep.vals=np.linspace(.48,.52,21)
#########    
###########################    sweep.vals=np.linspace(0,20,11)
#####    sweep.vals=np.linspace(.112,.114,21)
######################## 
#####################    sweep.vals=np.linspace(750,850,11)
#    freq=[]
#    dec=[]
#    fpes=[]
#    pi_pe=[]
#    
#    for ink,ssb_ge in enumerate(sweep.vals):
#            print(ink)
##            pi_tomo=0
##            black_nonHermitian_ma.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=0,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=ss,coef_in=1)
##            keithley2401.set_current(ss, step_size_mA=0.001)
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,
##                            rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.5,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
#            black_nonHermitian_ma.ramsey(ssm_ge = ssb_ge,off_set=0,num_steps = seq.num_patterns,t1_time = seq.sweep_time,pi_ge=p_ge,amp=.5) 
##            black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ss,off_set=0,t1_time = seq.sweep_time,num_steps =seq.num_patterns,pi_ge=p_ge,pi_ef=p_ef)
##         black_nonHermitian_ma.t1(ssm_ge = ssb_ge,off_set=0,num_steps=seq.num_patterns,t1_time = seq.sweep_time,t_min=0,pi_ge=p_ge,amp_g=.5)
##            pi_tomo=0
##            black_nonHermitian_ma.encirclingEP(rabi_time = ss,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=0,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=0.005,amp_enc=.25,ph=10,coef_in=1)
#        ########    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#        ########                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=180,off_set=0,amp_tomo=0.5
#        ########                                                 ,pi_time=pi_tomo,min_amp=.01,max_amp=.5,detun=0.003,amp_enc=.17,ph=15)
#        #########
##            pi_tomo=p_ef
##            black_nonHermitian_ma.encirclingEP(rabi_time = ss,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=seq.num_patterns,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=0.005,amp_enc=.25,ph=10,coef_in=1)
##            black_nonHermitian_ma.rabi(ssm_ge = ssb_ge, num_steps = seq.num_patterns,amp = .5,rabi_time = seq.sweep_time)  
##            black_nonHermitian_ma.ramsey_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,t2_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ss,pi_hf=p_hf,xp=0)
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=ss,
##                                          rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.498,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
##            black_nonHermitian_ma.rabi_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=ss,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf)
###
#            black_nonHermitian_ma.loading(num_steps = 1*seq.num_patterns)
#            wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
#            m = [0, 152.5];    
#
#            m= [0, 155.62557200692163];     
#            for k in range(seq.num_avgs):
#                
#                daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=m,
#                                    num_patterns=1*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=2)
#        #        
#
#        #    
#                prev_threshold=daq_params.threshold
#
#                if k is 0:
#                    p_readout = p
#        #                   
#                else:
#                    p_readout += p
##            p_post=analysis.p_readout_postselected_pief(p_readout)
#        
#            p_post = analysis.p_readout_postselected(p_readout)
#
#            x=seq.times
##            y=p_post[2]
#            y=p_readout[1]
#            fpes.append(y)   
##            plt.plot(x,y)
#            popt,peer, y_vals, _ = analysis.fit_sine_decay(x,y,guess_vals=None) 
#            popt,peer, y_vals, _ = analysis.fit_exp_decay(x,y,guess_vals=popt)
#
#            freq.append(popt[0])
#            pi_pe.append(1/2/popt[0])
#            dec.append(popt[1])

###             
#
##################################################################
#    sweep.vals=np.linspace(7.38,7.42,21)
##    sweep.vals=np.linspace(.408,.410,21)
#    sweep.vals=np.linspace(.219,.221,21)
##    sweep.vals=np.linspace(.76,.78,11 )
#    sweep.vals=np.linspace(.112,.114,21)
##    sweep.vals=np.linspace(.5,.54,21)
######    
########################    sweep.vals=np.linspace(0,20,11)
##
###################### 
###################    sweep.vals=np.linspace(750,850,11)
#    freq=[]
#    dec=[]
#    fpfs=[]
#    pi_pf=[]
#    
#    m = [147.985192421315, 152.9385192421315];      
#    for ink,ss in enumerate(sweep.vals):
#            print(ink)
##            pi_tomo=0
##            black_nonHermitian_ma.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=0,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=ss,coef_in=1)
##            keithley2401.set_current(ss, step_size_mA=0.001)
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,
##                            rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.5,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
##            black_nonHermitian_ma.ramsey(ssm_ge = ss,off_set=0,num_steps = seq.num_patterns,t1_time = seq.sweep_time,pi_ge=p_ge,amp=.5) 
#            black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ss,off_set=0,t1_time = seq.sweep_time,num_steps =seq.num_patterns,pi_ge=p_ge,pi_ef=p_ef)
##         black_nonHermitian_ma.t1(ssm_ge = ssb_ge,off_set=0,num_steps=seq.num_patterns,t1_time = seq.sweep_time,t_min=0,pi_ge=p_ge,amp_g=.5)
##            pi_tomo=0
##            black_nonHermitian_ma.encirclingEP(rabi_time = ss,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=0,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=0.005,amp_enc=.25,ph=10,coef_in=1)
##        ########    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##        ########                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=180,off_set=0,amp_tomo=0.5
##        ########                                                 ,pi_time=pi_tomo,min_amp=.01,max_amp=.5,detun=0.003,amp_enc=.17,ph=15)
##        #########
##            pi_tomo=p_ef
##            black_nonHermitian_ma.encirclingEP(rabi_time = ss,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=seq.num_patterns,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=0.005,amp_enc=.25,ph=10,coef_in=1)
##            black_nonHermitian_ma.rabi(ssm_ge = ssb_ge, num_steps = seq.num_patterns,amp = ss,rabi_time = seq.sweep_time)  
##            black_nonHermitian_ma.ramsey_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,t2_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ss,pi_hf=p_hf,xp=0)
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=ss,
##                                          rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.504,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
##            black_nonHermitian_ma.rabi_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=ss,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf)
###
#            black_nonHermitian_ma.loading(num_steps = 1*seq.num_patterns)
#            wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
#            m = [147.985192421315, 152.9385192421315];  
##            [145.12402489015983, 152.87402489015983];     
#            for k in range(seq.num_avgs):
#                
#                daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=m,
#                                    num_patterns=1*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=1,fg=3)
#        #        
#
#        #    
#                prev_threshold=daq_params.threshold
#
#                if k is 0:
#                    p_readout = p
#        #                   
#                else:
#                    p_readout += p
##            p_post=analysis.p_readout_postselected_pief(p_readout)
#        
#            p_post = analysis.p_readout_postselected(p_readout)
#
#            x=seq.times
#            y=p_post[2]
##            y=p_readout[1]
#            fpfs.append(y)   
#            plt.plot(x,y)
#            popt,peer, y_vals, _ = analysis.fit_sine_decay(x,y,guess_vals=None) 
##            popt,peer, y_vals, _ = analysis.fit_exp_decay(x,y,guess_vals=popt)
##
#            freq.append(popt[0])
#            pi_pf.append(1/2/popt[0])
#            dec.append(popt[1])
####
######             
#####################################################################
##    sweep.vals=np.linspace(7.38,7.42,21)
##    sweep.vals=np.linspace(.408,.410,21)
#    sweep.vals=np.linspace(.219,.221,21)
##    seq.sweep_time=2000
##    sweep.vals=np.linspace(.74,.78,21 )
##    sweep.vals=np.linspace(.48,.5,11)
#####    
#######################    sweep.vals=np.linspace(0,20,11)
###    sweep.vals=np.linspace(.112,.114,21)
##################### 
##################    sweep.vals=np.linspace(750,850,11)
#    freq=[]
#    dec=[]
#    fphs=[]
#    pi_ph=[]
#    
#    for ink,ss in enumerate(sweep.vals):
#            print(ink)
##            pi_tomo=0
##            black_nonHermitian_ma.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=0,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=ss,coef_in=1)
##            keithley2401.set_current(ss, step_size_mA=0.001)
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,
##                            rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.5,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
##            black_nonHermitian_ma.ramsey(ssm_ge = ss,off_set=0,num_steps = seq.num_patterns,t1_time = seq.sweep_time,pi_ge=p_ge,amp=.5) 
##            black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ss,off_set=0,t1_time = seq.sweep_time,num_steps =seq.num_patterns,pi_ge=p_ge,pi_ef=p_ef)
##         black_nonHermitian_ma.t1(ssm_ge = ssb_ge,off_set=0,num_steps=seq.num_patterns,t1_time = seq.sweep_time,t_min=0,pi_ge=p_ge,amp_g=.5)
##            pi_tomo=0
##            black_nonHermitian_ma.encirclingEP(rabi_time = ss,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=0,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=0.005,amp_enc=.25,ph=10,coef_in=1)
##        ########    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##        ########                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=180,off_set=0,amp_tomo=0.5
##        ########                                                 ,pi_time=pi_tomo,min_amp=.01,max_amp=.5,detun=0.003,amp_enc=.17,ph=15)
##        #########
##            pi_tomo=p_ef
##            black_nonHermitian_ma.encirclingEP(rabi_time = ss,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=seq.num_patterns,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=0.005,amp_enc=.25,ph=10,coef_in=1)
##            black_nonHermitian_ma.rabi(ssm_ge = ssb_ge, num_steps = seq.num_patterns,amp = ss,rabi_time = seq.sweep_time)  
#            black_nonHermitian_ma.ramsey_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,t2_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ss,pi_hf=p_hf,xp=0)
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=ss,
##                                          rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.498,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
##            black_nonHermitian_ma.rabi_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=ss,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf)
###
#            black_nonHermitian_ma.loading(num_steps = 1*seq.num_patterns)
#            wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[.033,.01, 0.033, -0.077])
#            m = [148.385192421315, 155.385192421315];     
##            [145.12402489015983, 152.87402489015983];     
#            for k in range(seq.num_avgs):
#                
#                daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=m,
#                                    num_patterns=1*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        #        
#
#        #    
#                prev_threshold=daq_params.threshold
#
#                if k is 0:
#                    p_readout = p
#        #                   
#                else:
#                    p_readout += p
##            p_post=analysis.p_readout_postselected_pief(p_readout)
#        
#            p_post = analysis.p_readout_postselected(p_readout)
#
#            x=seq.times
#            y=p_post[2]
##            y=p_readout[1]
#            fphs.append(y)   
#            plt.plot(x,y)
#            popt,peer, y_vals, _ = analysis.fit_sine_decay(x,y,guess_vals=None) 
##            popt,peer, y_vals, _ = analysis.fit_exp_decay(x,y,guess_vals=popt)
##
#            freq.append(popt[0])
#            pi_ph.append(1/2/popt[0])
#            dec.append(popt[1])

##             
#################################################################
#    sweep.vals=np.linspace(.383,.385,9)
#    sweep.vals=np.linspace(.244,.246,9)
#    sweep.vals=np.linspace(0.45,.55,11)
#    sweep.vals=np.linspace(.95,1.05,21)
##    sweep.vals=np.linspace(0,20,11)
#
#####    sweep.vals=np.linspace(.086,.089,10)
#    freq=[]
#    fp=[]
#    for ink,ss in enumerate(sweep.vals):
#            print(ink)
#            
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.25,
##                                          rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.5,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1,ef_dur=p_ef,coef_f=ss,phas=12)
#
#
##            black_nonHermitian_ma.rabi_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,amp_hf=ss,hf_dur=p_hf/2)
#
##            pi_tomo=0
#
##            black_nonHermitian_ma.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
##                                               ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=10,off_set=0,
##                                               amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=-0.005,amp_enc=.25,ph=ss,coef_in=)
#            keithley2401.set_current(ss, step_size_mA=0.001)
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ss,off_set=0,num_steps = seq.num_patterns,amp=.5,
##                            rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.5,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
##            black_nonHermitian_ma.ramsey(ssm_ge = ss,off_set=0,num_steps = seq.num_patterns,t1_time = seq.sweep_time,pi_ge=p_ge,amp=.5) 
##            black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ss,off_set=0,t1_time = seq.sweep_time,num_steps =seq.num_patterns,pi_ge=p_ge,pi_ef=p_ef)
##            black_nonHermitian_ma.ramsey_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,t2_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ss,pi_hf=p_hf)
##            black_nonHermitian_ma.loading(num_steps = 1*seq.num_patterns)
##            wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5])
#            [152.12402489015983, 156.87402489015983];     
#            for k in range(seq.num_avgs):
#                
#                daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
#                    num_patterns=1*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
#        #        
#
#        #    
#                prev_threshold=daq_params.threshold
#
#                if k is 0:
#                    p_readout = p
#        #                   
#                else:
#                    p_readout += p
#        
#            p_post = analysis.p_readout_postselected(p_readout)
#
#            x=seq.times
##            y=p_post[2]
#            y=p_readout[1]
#            fp.append(y)   
#            plt.plot(x,y)

################################################################################









# -*- coding: utf-8 -*-
#"""
#Created on Tue Mar 10 14:14:04 2020
#
#@author: crow104
#"""
#import numpy as np
#import matplotlib.pyplot as plt
#import pickle
#import datetime
#    
#import expt_parameters
#import seq_experiments
#import seq_programs
#import daq_programs
#import analysis
#import black_nonHermitian_ma
#
#from Nop_class import Nop
#
#import wx_programs
#import keithley2401
#
#
#if __name__ == '__main__':
#    expt = expt_parameters.expt_parameters()
#    seq = Nop("seq")
#    msmt = Nop("measurement")
#    sweep = Nop("sweep")
#    save_dir = r"C:\Data\2020\200309_nonHermitian\Ering data\loop/"
##    ssb_ge=.3865
##    ssb_ef=.0912
#    ssb_ge=.3906
#    ssb_ef=0.0945
#    p_ge=36
#    p_ef=28
    
    
    

#    seq.comment = "ramsey_ge calibration"
#    seq.num_patterns = 51
#    seq.sweep_time = 1000
#    seq.num_records_per_pattern = 500
#    seq.rabi_amp = 0.05
#    seq.times = np.linspace(0., seq.sweep_time, seq.num_patterns)*1e-3
#    steps=13
#    
#    sweep = Nop("sweep")
#    sweep.comment = "rabi ssb freq near EP"
#    sweep.vals = np.linspace(-3, 3, steps)*1e-3 + 0.092 # ge ramsey
##    sweep.vals = np.linspace(-0.5, 0.5, 3)*1e-3 + expt_cal.ssm.ef
#    pop_3state = []
#    fpop=[]
#
#    for idx, ssm_ef in enumerate(sweep.vals):
#        print(idx)
#        
#        black_nonHermitian_ma.ramsey_ef(ssm_ef)
#        wx_programs.wx_set_and_amplitude_and_offset()
#
##        seq_experiments.rabi_ef_prep_f(seq.num_patterns, seq.sweep_time, seq.rabi_amp, rabi_ssb_freq)
#        daq_params, rec_readout_vs_pats, p_readout = daq_programs.run_daq(
#                num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern)
#        p_post = analysis.p_readout_postselected(p_readout)
#        pop_3state.append(p_readout)
##        msmt.popt, msmt.perr, _, _ = analysis.fit_sine_decay(x, y, guess_vals=None)
#        x = seq.times
#        y = p_post[1]
##        analysis.fit_sine_decay(x,y,)
#        fpop.append(p_post[1])
##        save_by_pickle((expt, seq, daq, msmt))
#        plt.plot(seq.times, p_post[1])
#        plt.ylim([0,4])
#        plt.show()
##    pop_3state = np.stack(pop_3state)
#    popdstack=np.dstack(pop_3state)
#
#    pop3state_f=popdstack.reshape(seq.num_patterns*3,steps)
#    np.savetxt(save_dir+'popdstack', pop3state_f)
#    np.savetxt(save_dir+'fpop', fpop)
#    fig = plt.figure()
#    plt.imshow(fpop)
##    plt.imshow(pop_3state[:,1,:])
#    plt.show
       ##########################################
#   
#    seq.comment = "ep"
#    seq.num_patterns = 51
#    seq.sweep_time = 2000
#    seq.num_records_per_pattern = 3000
#    seq.rabi_amp = 0.05
#    seq.times = np.linspace(0., seq.sweep_time, seq.num_patterns)*1e-3
#    steps=50
#    
#    sweep = Nop("sweep")
#    sweep.comment = "rabi ssb freq near EP"
#    
##    sweep.vals = np.linspace(-1, 1, steps)*1e-3 + 0.091 # ge ramsey
#    sweep.vals = np.linspace(.2, .05, steps) # ge ramsey
##    sweep.vals = np.linspace(-0.5, 0.5, 3)*1e-3 + expt_cal.ssm.ef
#    pop_3state = []
#    f_pop=[]
#    g_pop=[]
#    e_pop=[]
#    ep_pop=[]
#    fp_pop=[]
#    dummy=[]
#    freq_ef=[]
#    popt1=[  0.87047485,   4.57162552,  -0.11505249, -99.82764558,  0.74976633]
#    
#
#    for idx, amp in enumerate(sweep.vals):
#        print(idx)
#        
##  
##        black_nonHermitian_ma.ramsey_ef(ssm_ge=ssb_ge, ssm_ef = ssm_efs)
#
#        wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, amp, amp])
#
##        seq_experiments.rabi_ef_prep_f(seq.num_patterns, seq.sweep_time, seq.rabi_amp, rabi_ssb_freq)
#        daq_params, rec_readout_vs_pats, p_readout = daq_programs.run_daq(
#                num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern)
#        p_post = analysis.p_readout_postselected(p_readout)
#        g_pop.append(p_readout[0])
#        e_pop.append(p_readout[1])
#        f_pop.append(p_readout[2])
#        x=seq.times
#        y=p_post[1]
#        popt1, perr1, _, _ = analysis.fit_sine_decay(x, y, guess_vals=popt1)
#        freq_ef.append(abs(popt1[0]))
#
#
##        analysis.fit_sine_decay(x,y,)
#
#        
#        ep_pop.append(p_post[1])
#        fp_pop.append(p_post[2])
##        print(np.shape(ep_pop))
##        save_by_pickle((expt, seq, daq, msmt))
##        plt.plot(seq.times, p_post[1])
##        plt.ylim([0,1])
##        plt.show()
##    pop_3state = np.stack(pop_3state)
##    popdstack=np.dstack(pop_3state)
##
##    pop3state_f=popdstack.reshape(seq.num_patterns*3,steps)
##    np.savetxt(save_dir+'popdstack', pop3state_f)
##    np.savetxt(save_dir+'fpop', fpop)
#    fig = plt.figure()
#    plt.imshow(fp_pop)
#    fig = plt.figure()
#    plt.plot(sweep.vals,freq_ef,'o')
#    #        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
#
##    plt.imshow(pop_3state[:,1,:])
#    plt.show
   ###########################################
#   
#    seq.comment = "rabi_ge calibration of drive amplitude"
#    seq.num_patterns = 51
#    seq.sweep_time = 4000
#    seq.num_records_per_pattern = 500
#    seq.rabi_amp = 0.05
#    seq.times = np.linspace(0., seq.sweep_time, seq.num_patterns)*1e-3
#    steps=9
#    
#    sweep = Nop("sweep")
#    sweep.comment = "rabi ssb freq near EP"
#    
#    sweep.vals = np.linspace(-1, 1, steps)*1e-3 + 0.091 # ge ramsey
##    sweep.vals = np.linspace(-0.5, 0.5, 3)*1e-3 + expt_cal.ssm.ef
#    pop_3state = []
#    f_pop=[]
#    g_pop=[]
#    e_pop=[]
#    ep_pop=[]
#    fp_pop=[]
#    dummy=[]
#    freq_ge=[]
#
#    for idx, ssm_efs in enumerate(sweep.vals):
#        print(idx)
#        
##  
#        black_nonHermitian_ma.ramsey_ef(ssm_ge=ssb_ge, ssm_ef = ssm_efs)
#
#        wx_programs.wx_set_and_amplitude_and_offset()
#
##        seq_experiments.rabi_ef_prep_f(seq.num_patterns, seq.sweep_time, seq.rabi_amp, rabi_ssb_freq)
#        daq_params, rec_readout_vs_pats, p_readout = daq_programs.run_daq(
#                num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern)
#        p_post = analysis.p_readout_postselected(p_readout)
#        g_pop.append(p_readout[0])
#        e_pop.append(p_readout[1])
#        f_pop.append(p_readout[2])
#        x=seq.times
#        y=p_post[1]
#        popt1, perr1, _, _ = analysis.fit_sine_decay(x, y, guess_vals=None)
#        freq_ge.append(abs(popt1[0]))
#
#
##        analysis.fit_sine_decay(x,y,)
#
#        
#        ep_pop.append(p_post[1])
#        fp_pop.append(p_post[2])
##        print(np.shape(ep_pop))
##        save_by_pickle((expt, seq, daq, msmt))
#        plt.plot(seq.times, p_post[1])
#        plt.ylim([0,1])
#        plt.show()
##    pop_3state = np.stack(pop_3state)
##    popdstack=np.dstack(pop_3state)
##
##    pop3state_f=popdstack.reshape(seq.num_patterns*3,steps)
##    np.savetxt(save_dir+'popdstack', pop3state_f)
##    np.savetxt(save_dir+'fpop', fpop)
#    fig = plt.figure()
#    plt.imshow(g_pop)
#    fig = plt.figure()
#    plt.plot(sweep.vals,freq_ge)
#    #        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
#
##    plt.imshow(pop_3state[:,1,:])
#    plt.show
   ##########################################
########    parameters
#    steps=61
##    sweep.vals = np.linspace(-2, 2, steps)*1e-3  # ge ramsey
#    sweep.vals = np.linspace(0, .1, steps)
#
#
##
#    for idx, amps in enumerate(sweep.vals):
#        print(idx)
#        off_sets=idx*51
#
##        
#        black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge, ssm_ef = ssb_ef,off_set=off_sets,amp=amps)  
##        black_nonHermitian_ma.ramsey(ssm_ge = ssb+ssm,off_set=off_sets)
##        black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge, ssm_ef = ssb_ef+ssm,off_set=off_sets)
#    black_nonHermitian_ma.loading(num_steps = 61*51)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 0.2, 0.2])
   
   ##########################################
######    parameters
#    steps=51
#    sweep.vals = np.linspace(0, .15, steps)  # ge ramsey
#
#
##
#    for idx, amps in enumerate(sweep.vals):
#        print(idx)
#        off_sets=idx*51
#
#        
##        black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=off_sets,amp=amps)
#        black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=off_sets,amp=amps)
#        
#        
   ##########################################
###    parameters
#    seq.comment = "One point measurement of Ering tVar"
#    steps=101
#    tmax=8000
#    tmin=100
#    sweep.vals = np.linspace(tmin, tmax, steps)  # ge ramsey
#    ampe=.05
#
#
#    print('x')
#    for idx, rabi_t in enumerate(sweep.vals):
#        print(idx)
#        off_sets=idx
##        black_nonHermitian_ma.xtransport_ef_tVar(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=off_sets,ph=90,ph2=0,amplz=.5,epamp=amp,rabi_time=rabi_t)
#        black_nonHermitian_ma.xtransport_ef_tVar(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=off_sets,ph=90,ph2=0,amplz=.5,epamp=ampe,rabi_time=rabi_t,pi_ge=p_ge,pi_ef=p_ef)
#    print('y')
#    for idx, rabi_t in enumerate(sweep.vals):
#        print(idx)
#        off_sets=idx
##        black_nonHermitian_ma.xtransport_ef_tVar(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=51+off_sets,ph=90,ph2=90,amplz=.5,epamp=amp,rabi_time=rabi_t)
#        black_nonHermitian_ma.xtransport_ef_tVar(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=steps+off_sets,ph=90,ph2=90,amplz=.5,epamp=ampe,rabi_time=rabi_t,pi_ge=p_ge,pi_ef=p_ef)
#    print('z')
#    for idx, rabi_t in enumerate(sweep.vals):
#        print(idx)
#        off_sets=idx
##        black_nonHermitian_ma.xtransport_ef_tVar_z(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=102+off_sets,ph=0,ph2=0,amplz=.5,epamp=amp,rabi_time=rabi_t)
#        black_nonHermitian_ma.xtransport_ef_tVar_z(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=2*steps+off_sets,ph=0,ph2=0,amplz=.5,epamp=ampe,rabi_time=rabi_t,pi_ge=p_ge,pi_ef=p_ef)
#   
#   
#
#    black_nonHermitian_ma.t1(ssm_ge = ssb_ge,off_set=steps*3,num_steps = steps,t1_time = 4000,t_min=0,pi_ge=p_ge)
#    black_nonHermitian_ma.t1_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=steps*4,num_steps = steps,t1_time = 4000,t_min=0,pi_ge=p_ge,pi_ef=p_ef)
#    black_nonHermitian_ma.ramsey(ssm_ge = ssb_ge,off_set=steps*5,num_steps = steps,t1_time = 1000,pi_ge=p_ge)
#    black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=steps*6,t1_time = 1000,num_steps = steps,pi_ge=p_ge,pi_ef=p_ef)
#    
##    black_nonHermitian_ma.t1(ssm_ge = 0.3865,off_set=51*3)
##    black_nonHermitian_ma.t1_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=51*4)
#    black_nonHermitian_ma.loading(num_steps = steps*7)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 0.2, 0.2])
###        
####    
#   
#  ###################################################  
#    seq.comment = "measuring Ring"
#    seq.num_patterns = 101*7
#    seq.sweep_time = 1000
#    seq.num_records_per_pattern = 1000
#    seq.num_avgs=5
#    
#    
#
#    
#
#    
#    
##    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, .2, .2])    
#    for k in range(seq.num_avgs):
#    #
#        daq_params, _, p = daq_programs.run_daq(
#            num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern)
#        if k is 0:
#            p_readout = p
#        else:
#            p_readout += p
#            
#    p_readout = p_readout/seq.num_avgs
#    p_post = analysis.p_readout_postselected(p_readout)
#    
#
#    np.savetxt(save_dir+"test2"+str(0.2), p_readout)
#    np.savetxt(save_dir+"p2test"+str(0.2), p_post)
####################################
##    
#    seq.comment = "ES_unbroken"
#    seq.num_patterns = 51*51
#    seq.sweep_time = 1000
#    seq.num_records_per_pattern = 50
#    seq.num_avgs=50
#    
#    
#
#    
#
#    
#    
#    
#    for k in range(seq.num_avgs):
#    #
#        daq_params, _, p = daq_programs.run_daq(
#            num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern)
#        if k is 0:
#            p_readout = p
#        else:
#            p_readout += p
#            
#    p_readout = p_readout/seq.num_avgs
#    p_post = analysis.p_readout_postselected(p_readout)
#    np.savetxt(save_dir+'p_readout', p_readout)
#    np.savetxt(save_dir+'p_post', p_post)
##
##        
#  ################################
#    seq.comment = "quantum state tomography of the Ering"
#    seq.num_patterns = 51*3
#    seq.sweep_time = 4000
#    seq.num_records_per_pattern = 5000
#    seq.times = np.linspace(0., seq.sweep_time, seq.num_patterns)*1e-3
#    amp=.03
#    
#    
#
#    
#
#    
#    
##    
##    black_nonHermitian_ma.xtransport_ef_z(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,ph=90,ph2=90,amplz=0.5,epamp=amp)
##    black_nonHermitian_ma.xtransport_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=51,ph=90,ph2=90,amplz=0.5,epamp=amp)
##    black_nonHermitian_ma.xtransport_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=102,ph=90,ph2=0,amplz=0.5,epamp=amp)
#        
#    black_nonHermitian_ma.xtransport_ef_z(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,ph=90,ph2=90,amplz=0.5,epamp=amp,rabi_time=seq.sweep_time)
#    black_nonHermitian_ma.xtransport_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=51,ph=90,ph2=90,amplz=0.5,epamp=amp,rabi_time=seq.sweep_time)
#    black_nonHermitian_ma.xtransport_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=102,ph=90,ph2=0,amplz=0.5,epamp=amp,rabi_time=seq.sweep_time)
#    black_nonHermitian_ma.loading(num_steps = 51*3)
#     
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 0.2, 0.2])
#    #
#    daq_params, t_histo, p_readout, = daq_programs.run_daq(
#            num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern)
#    p_post = analysis.p_readout_postselected(p_readout)
#    x=seq.times
#    y=p_post[2]  
##    analysis.fit_exp_decay
#    popt,peer, _, _ = analysis.fit_sine_decay(x,y,guess_vals=None)      
#    plt.plot(x,y)
#
#
#    np.savetxt(save_dir+'r_tomozxy_ring_amp03t4us', p_readout)
#    np.savetxt(save_dir+'p_tomozxy_ring_amp03t4us', p_post)
   #  ################################
#    seq.comment = "measurement"
#    seq.num_patterns = 101
#    seq.sweep_time = 3000
#    seq.num_records_per_pattern = 1000
#    seq.times = np.linspace(0., seq.sweep_time, seq.num_patterns)*1e-3
#    amp=.13
#    
#    
#
#    
#
##    
##    
##    
##    black_nonHermitian_ma.t1(ssm_ge = ssb_ge,off_set=0,num_steps=seq.num_patterns,t1_time = 1000,t_min=0)
## #   black_nonHermitian_ma.t1_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=51)
###    black_nonHermitian_ma.ramsey(ssm_ge = 0.3885,off_set=0)
###    black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,t1_time = 6000,num_steps =seq.num_patterns)
###    black_nonHermitian_ma.loading(num_steps = seq.num_patterns)
##    black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = 100)
#    black_nonHermitian_ma.ramsey(ssm_ge = ssb_ge,off_set=0,num_steps = seq.num_patterns,t1_time = 3000,pi_ge=p_ge)
#    black_nonHermitian_ma.loading(num_steps = seq.num_patterns)
#    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, .2, .2])
#    #
#    daq_params, t_histo, p_readout, = daq_programs.run_daq(
#            num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern)
##    p_post = analysis.p_readout_postselected(p_readout)
#    x=seq.times
#    y=p_readout[0]  
###    analysis.fit_exp_decay
#    popt,peer, _, _ = analysis.fit_sine_decay(x,y,guess_vals=None)      
##    plt.plot(x,y)
##
##
#    np.savetxt(save_dir+'tp_readout', p_readout)
#    np.savetxt(save_dir+'tp_post', p_post)
#######################################

#    seq.num_patterns = 101*5
#    seq.sweep_time = 1000
#    seq.num_records_per_pattern = 1000
#    seq.times = np.linspace(0., seq.sweep_time, seq.num_patterns)*1e-3
# 
#    steps=3
#    seq.num_avgs=2
#    
#
#    sweep.comment = "average of amp change TVar"
#    
#    sweep.vals = np.linspace(.1, .3, steps)
#
#    pop_3state = []
#    f_pop=[]
#    g_pop=[]
#    e_pop=[]
#    ep_pop=[]
#    fp_pop=[]
#    dummy=[]
#    freq_ge=[]
#
#    for idx, amp in enumerate(sweep.vals):
#        print(idx)
#        amppi=amp
#        
##
#
#        wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, amp, amp])
#
#   
#    
#        for k in range(seq.num_avgs):
#        #
#            daq_params, _, p = daq_programs.run_daq(
#                num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern)
#            if k is 0:
#                p_readout = p
#            else:
#                p_readout += p
#                
#        p_readout = p_readout/seq.num_avgs
#        p_post = analysis.p_readout_postselected(p_readout)
#        np.savetxt(save_dir+"r_tomoxyz_tVar8_.1to5us_d6t"+str(idk), p_readout)
#        np.savetxt(save_dir+"p_tomoxyz_tVar8_.1to5us_d6t"+str(idx), p_post)
#
#
#        g_pop.append(p_readout[0])
#        e_pop.append(p_readout[1])
#        f_pop.append(p_readout[2])
#
#
#
#        ep_pop.append(p_post[1])
#        fp_pop.append(p_post[2])
#
#        plt.plot(seq.times, p_post[1])
#        plt.ylim([0,1])
#        plt.show()
#        
#    fig = plt.figure()
#    plt.imshow(g_pop)
#    np.savetxt(save_dir+'fp_pop',fp_pop)
#    np.savetxt(save_dir+'ep_pop',ep_pop)
#    np.savetxt(save_dir+'f_pop',f_pop)
#    np.savetxt(save_dir+'g_pop',g_pop)
#    np.savetxt(save_dir+'e_pop',e_pop)
#    
#    
#    
#
#    plt.show