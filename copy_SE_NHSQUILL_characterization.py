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
ssm_geq2 = -0.1455-0.00045 #-0.1538 + 0.00016#-.162 #-0.1468
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


#Scaling
scale_matrix1 = np.array([[9.34000000e-01, 1.44800000e-01, 0],
                          [6.58000000e-02, 8.55000000e-01, 0],
                          [0, 0, 1]])

scale_matrix1 = np.linalg.inv(scale_matrix1)
#ge_matrix1 = np.array([[9.40866667e-01, 1.92466667e-01],
#                       [5.90000000e-02, 8.07333333e-01]])

#scale_matrix2 = np.array([[0.9904    , 0.1262    , 0.],
#                          [0.0096    , 0.8738    , 0.],
#                          [0.        , 0.        , 1        ]])

#scale_matrix2 = np.array([[0.008   , 0.8884    , 0.],
#                          [0.992    , 0.1116    , 0.],
#                          [0.        , 0.        , 1        ]])
scale_matrix2 = np.array([[0.9904    , 0.1262    , 0.],
                          [0.0096    , 0.8738    , 0.],
                          [0.        , 0.        , 1        ]])    

scale_matrix2 = np.linalg.inv(scale_matrix2)

    
#qubit_1_thr= [-1500,300]#[-800,550]
#qubit_2_thr=[-3000,-1200]#[500,1700]#
#IQ_angle_q1 = -90 #130
#IQ_angle_q2 = 240 #increasing rotates CCW
#pi_ef_time = 40.46

# set 1 for which sequence to run
spec_ge=0
spec_ef= 0
cluster_thr = 0
pis_nopi = 0
pi_2_nopi_2 = 0
no_pi_pm = 0
n_pi_2 = 0
run_rabi =0
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
run_sweep =0
run_J_sweep = 0
run_bnc_sweep = 0
run_bnc_sweep_2q = 0
run_bnc_power_sweep = 0
run_bnc_pow_freq_sweep = 0
run_ssbef_sweep = 0 #chevron for ssb_ef
run_x_measure = 0
run_y_measure = 0
run_tomo_calib = 0
teaching=0
ef_teach=0
parametric_coupling=0
parametric_coupling_time_domain=0
parametric_coupling_sweep = 0
run_ramsey_coax_drive = 0
run_T1_ge_raman_cooling = 0
#run_T1_ge_tunable_dissipation = 0
run_qubit_spectroscopy_vs_cavity_drive = 0
run_spec_ROcavity = 0
run_spec_ROcavity_powerSweep = 0
run_raman_cavityPowerSweep = 1
run_raman_cavityDetuningSweep = 0
run_sweep_keithley = 0
run_ramsey_raman = 0


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
    ROIF1 = ROq1 - RO_LO
    ROIF2 = ROq2 - RO_LO
    RO_freq = RO_LO
    #RO_freq = ROq2
   # qubit_bnc =4.3#4.256
    bnc.set_bnc_output(RO_freq, bnc_addr=target_bnc_black_address)
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

if pis_nopi: #run to correct two qubit rabi tomography
    num_steps = 3
    reps = 15000
    if which_qubit == 0:
        
        Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|g>, keeps population in ground state
        Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,ROIF=ROIF2,ssm_ge = ssm_ge); #|e>, sends population to |e>
        
    if which_qubit == 1:
        Fs.pi_nopi(off = 0,coef=0,pi_ge=pi_ge_time,q = which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge); #|g>, keeps population in ground state
        Fs.pi_nopi(off = 1,coef=1,pi_ge=pi_ge_time,q= which_qubit,ROIF=ROIF1,ssm_ge = ssm_ge); #|e>, sends population to |e>
        

if spec_ge:
    num_steps = 101
    f1=-.170 #-0.17
    f2= -.140#-0.13
    if which_qubit == 0:  #readout on both qubits, excited q2
        Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.01,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
    if which_qubit == 1: ##readout on both qubits, excited q1
        Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.1,ROIF1=ROIF1,ROIF2=ROIF2,q=which_qubit)
   


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
    num_steps = 51
    sweep_time =100 #ns
    if which_qubit == 0:
        Fs.rabi_ge(num_steps,sweep_time,ssm_geq2,ROIF1,ROIF2,which_qubit)
    if which_qubit == 1:
        Fs.rabi_ge(num_steps,sweep_time,ssm_geq1,ROIF1,ROIF2,which_qubit)
        
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
    
if run_T1_ge_2q_RO:
    num_steps = 81
    sweep_time =800#ns
    Fs.T1_ge_2q_RO(num_steps,sweep_time,ssm_ge,pi_ge_time,ROIF1, ROIF2,q=which_qubit,ifload = 1)
       
if run_T1_M_ge_2q_RO:
    num_steps = 81
    sweep_time =100000#ns
    Fs.T1_M_ge_2q_RO(num_steps,sweep_time,ssm_ge,pi_ge_time,ROIF1, ROIF2,q=which_qubit,ifload = 1)
 
if run_T1_modified:
    num_steps = 51
    sweep_time =40000 #ns
    Fs.T1_ge_modified(num_steps,sweep_time)

if run_ramsey:
    piecho = 0
    osc_num = 4
    num_steps = 51
    sweep_time =20000 #ns
    if which_qubit == 0:
        Fs.ramsey(num_steps,sweep_time,piecho,osc_num,ssm_ge,pi_ge_time,ROIF2,q=which_qubit)
        #Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
    if which_qubit == 1:
        #Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF1,which_qubit)
        Fs.ramsey(num_steps,sweep_time,piecho,osc_num,ssm_ge,pi_ge_time,ROIF1,q=which_qubit)

if run_ramsey_raman:
    piecho = 0
    osc_num = 4
    num_steps = 51
    sweep_time =5000 #ns
    
    detun = -0.03 
    ssm_coax = ROIF2 + detun
    amp_coax = 1
    ssm_Q2 = ssm_ge + detun #- 0.0055 # +0.001 for amp_Q2=0.5 # excitation -detun; decay + detun
    amp_Q2 = 0.5 #0.5 #0.54 
        
    if which_qubit == 0:
        Fs.ramsey_raman(num_steps,sweep_time,piecho,osc_num,ssm_ge,pi_ge_time,ssm_coax,amp_coax,ssm_Q2,amp_Q2,ROIF2,q=which_qubit)
        #Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)
    if which_qubit == 1:
        #Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF1,which_qubit)
        Fs.ramsey_raman(num_steps,sweep_time,piecho,osc_num,ssm_ge,pi_ge_time,ssm_coax,amp_coax,ssm_Q2,amp_Q2,ROIF1,q=which_qubit)
    
    
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


##    
#else:
#    #Do actual heterodyne
#    wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
#    n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
#       
##    print("Qubit 1: P_Q1")
##    P_Q1 = prob_vs_pats_1[0];plt.plot(P_Q1);plt.title('Q1 thresholded');plt.show()
##    I_Q1 = rec_avg_vs_pats_1[0];
##    Q_Q1 = rec_avg_vs_pats_1[1];
##    p1_readout= np.matmul(scale_matrix1,prob_vs_pats_1) # scaled
##    plt.plot(I_Q1);plt.title('I Q1');plt.show()
##    plt.plot(Q_Q1);plt.title('Q Q1');plt.show()
##
##    plt.plot(p1_readout[0],label='|g>');plt.plot(p1_readout[1],label='|e>')
##    plt.legend();plt.title('Q1 data scaled with matrix');plt.show()
#
#
###    
#    print("Qubit 2:")
#    P_Q2 = prob_vs_pats_2[0];plt.plot(P_Q2);plt.title('Q2 thresholded');plt.show()
#    I_Q2 = rec_avg_vs_pats_2[0];
#    Q_Q2 = rec_avg_vs_pats_2[1];
#    p2_readout= np.matmul(scale_matrix2,prob_vs_pats_2)    
#    plt.plot(I_Q2);plt.title('I Q2');plt.show()
#    plt.plot(Q_Q2);plt.title('Q Q2');plt.show()       
#
#    plt.plot(p2_readout[0],label='|g>');plt.plot(p2_readout[1],label='|e>')
#    plt.legend();plt.title('Q2 data scaled with matrix');plt.show()
#

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
#if which_qubit == 0:
#    I = I_Q2; Q = Q_Q2
#    y = P_Q2
#if which_qubit == 1:
#    I = I_Q1; Q = Q_Q1
#    y = P_Q1
    

if run_rabi:
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,P_Q2,guess_vals=[5, 0.5, -200, 1.1, -500]) # [6,0.05,0.07,90,0.5]
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
    Qrange = abs(np.max(Q_Q2)-np.min(Q_Q2))
    Irange = abs(np.max(I_Q2)-np.min(I_Q2))
    if Qrange>Irange:
        freq_index = np.where(Q_Q2 == np.max(Q_Q2))
        print("Q")
        plt.plot(freq,Q_Q2)
    if Irange>Qrange:
        freq_index = np.where(I_Q2 == np.max(I_Q2))
        print("I")
        plt.plot(freq,I_Q2)
#    freq_index = np.where(y == np.max(y))
    ssm_ge = freq[freq_index]
    print(ssm_ge, freq_index)
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
        T1_ge_fit_vals,_,_,_ = analysis.fit_exp_decay(times,Q_Q2,guess_vals=[-400,0.001,-1700])
    if Irange>Qrange:
        T1_ge_fit_vals,_,_,_ = analysis.fit_exp_decay(times,I_Q2,guess_vals=[400,1,-1700])
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
    
    
    
if run_T1ef:
    times = np.linspace(0,sweep_time/1000,num_steps)
    analysis.fit_exp_decay(times,y,guess_vals=[-4,0.08,0.5])
    
if run_T1_modified:
    times = np.linspace(0,sweep_time/1000,num_steps)
    analysis.fit_exp_decay(times,I,guess_vals=[-4,0.08,0.5])

if run_ramsey or run_ramsey_ef or run_ramsey_raman:
    times = np.linspace(0,sweep_time/1000,num_steps)
    Qrange = abs(np.max(Q_Q2)-np.min(Q_Q2))
    Irange = abs(np.max(I_Q2)-np.min(I_Q2))
    if Qrange>Irange:
        T2_fit_vals,_,_,_ = analysis.fit_sine_decay(times,Q_Q2,guess_vals=[0.8,0.4,0.08,-260,-1600])
    if Irange>Qrange:
        T2_fit_vals,_,_,_ = analysis.fit_sine_decay(times,I_Q2,guess_vals=[0.8,0.4,0.08,-260,-1600])
    T2_fit_vals,_,_,_ = analysis.fit_sine_decay(times,P_Q2,guess_vals=[0.4,1,180,90,0.5])
    T2 = 1/T2_fit_vals[1]
    print("T2* = {} \u03BCs".format(T2))
#    analysis.fit_exp_decay(times,Q,guess_vals=[5,1,120])
    
#### run keithley sweep ####
if run_sweep:
#    center_cur = -0.1
    #num_steps=101#may need to remove eventually
    out_Q1, out_I1,  out_Q2, out_I2 = chevron.sweep_keithley(11.1,11.35,21,num_steps,reps,ro_dur,ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2,qubit_1_thr,qubit_2_thr)
    save_basename = '\specq1_amp_'+str(wx_amps[1])+'coupler_p5_2'+ str(wx_offs[1])+'.txt'
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
    start_freq = qubit_bnc
    stop_freq = qubit_bnc + 0.06
    out_Q, out_I = chevron.sweep_bnc_freq(start_freq, stop_freq,sweep_steps,num_steps,reps,ro_dur, IQ_angle)
    #np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q1_I',out_I)
    #np.savetxt(save_dir+'\sweep_qubit_BNC_freq_7dB_'+str(start_freq)+'to'+str(stop_freq)+'_q1_Q',out_Q)
    
if run_bnc_sweep_2q:
    sweep_steps = 51
    start_freq = 4.32- 0.002 #qubit_bnc - 0.005
    stop_freq = 4.32+ 0.002 #qubit_bnc + 0.005
    out_Q1, out_I1,out_Q2, out_I2 = chevron.sweep_bnc_freq_2q(start_freq, stop_freq,sweep_steps,num_steps,reps,ro_dur, ROIF1,ROIF2, deg_1 = 0, deg_2 = 0,verbose=True)
    np.savetxt(save_dir+'\coupler_p5_8usRabi_chev'+str(start_freq)+'to'+str(stop_freq)+'_q1_I',out_I1)
    np.savetxt(save_dir+'\coupler_p5_8usRabi_chev'+str(start_freq)+'to'+str(stop_freq)+'_q1_Q',out_Q1)
    np.savetxt(save_dir+'\coupler_p5_8usRabi_chev'+str(start_freq)+'to'+str(stop_freq)+'_q2_I',out_I2)
    np.savetxt(save_dir+'\coupler_p5_8usRabi_chev'+str(start_freq)+'to'+str(stop_freq)+'_q2_Q',out_Q2)
    
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
    
    
if run_sweep_keithley:

    start_current = -0.5
    stop_current = 0.5
    num_points = 21
    sweep_vals = np.linspace(start_current, stop_current, num_points)
    
    num_steps = 201
    f1 = -0.5
    f2 = -0.05 
    freq = np.linspace(f1,f2,num_steps)
    Fs.spectroscopy_ge(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.05, q = which_qubit)

#    ro_freq = np.linspace(6.889,6.88888,num_points)
    pop_e = np.zeros((num_points,num_steps))
    ssm_ge_tot = np.zeros(num_points)
    for i,ss in enumerate(sweep_vals):
        print(i)
        keithley2401.set_current(ss, step_size_mA=0.001)
#        input_resonator_freq = 6.88888 #7.0137 #ro_freq[i]
#        bnc.set_bnc_output(input_resonator_freq, bnc_addr=target_bnc_black_address) #in GHz
#        keithley2401.set_voltage(ss, step_size_mV=10)
        
        daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_pi,counts_pi,rec_readout,rec_avg_all,rec_all = daq_programs.run_daq(num_steps,reps,ro_dur,IQangle=100)
        
#        daq_params, t_histo, p,a,b,n_readout = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=[-120, -94],num_patterns=num_steps, 
#                                                                                               num_records_per_pattern=reps,authr=0,fg=2,ro_dur=ro_dur,IQangle=100)
        #daq_params, t_histo, p,a,b,n_readout = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold=m,num_patterns=num_steps, 
        #                                                                                       num_records_per_pattern=reps,authr=0,fg=2,ro_dur=ro_dur,IQangle=IQ_angle)
        
#        m = daq_params.threshold
#        y=p[1]+p[2]
#        pop_e[i] = y
        
        plt.plot(Q)
        plt.show()
        plt.plot(I)
        plt.show()

        ssm_ge_tot[i] = freq[np.argmax(Q)]
        print(ssm_ge_tot[i])
        

    plt.plot(sweep_vals,ssm_ge)
    plt.show()
    
    
if run_spec_ROcavity:
    sweep_steps = 21
    start_freq = 6.77 - 0.0005
    stop_freq = 6.77 + 0.0005
    
#    RO_LO = 6.77#6.85 #6.6


    freq = np.linspace(start_freq,stop_freq,sweep_steps)
# qubit at ground state
    Fs.spectroscopy_ROcavity(num_steps=3,ssm_ge=ssm_ge,pi_ge_time=pi_ge_time,spec_amp=.05*4,ROIF = ROIF2, estate = 0, q = 0)
#    out_IQ_g, out_phase_g = copy_chevron.sweep_bnc_freq_cavity_spec(start_freq, stop_freq, ROIF2, sweep_steps,IQ_angle,3,reps,ro_dur, wx_amps)
    out_IQ_g, out_phase_g = copy_chevron.sweep_bnc_freq_cavity_spec(start_freq,stop_freq,sweep_steps, ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2, 3,reps,ro_dur,qubit_1_thr,qubit_2_thr,wx_amps)


# qubit at excited state
    Fs.spectroscopy_ROcavity(num_steps=3,ssm_ge=ssm_ge,pi_ge_time=pi_ge_time,spec_amp=.05*4, ROIF = ROIF2, estate = 1, q = 0)
#    out_IQ_e, out_phase_e = copy_chevron.sweep_bnc_freq_cavity_spec(start_freq, stop_freq, ROIF2, sweep_steps,IQ_angle,3,reps,ro_dur, wx_amps)
    out_IQ_e, out_phase_e = copy_chevron.sweep_bnc_freq_cavity_spec(start_freq,stop_freq,sweep_steps, ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2, 3,reps,ro_dur,qubit_1_thr,qubit_2_thr,wx_amps)
        
    
    freq_g = freq[np.argmin(out_IQ_g)]
    freq_e = freq[np.argmin(out_IQ_e)]
    print(freq_g)
    print(freq_e)

    plt.plot(freq,out_IQ_g,'b-',label='g')
    plt.plot(freq,out_IQ_e,'r-',label='e')
    plt.xlabel('freq (GHz)')
    plt.ylabel('I^2+Q^2')
    plt.legend()
    plt.show()
    
    plt.plot(freq,out_phase_g,'b-',label='g')
    plt.plot(freq,out_phase_e,'r-',label='e')
    plt.xlabel('freq (GHz)')
    plt.ylabel('arctan(Q/I)')
    plt.legend()
    plt.show()
    
#    np.savetxt(save_dir+'ROcavity_spec_with_qubitGroundState_'+'probe_amp_0p05'+'_freq'+str(start_freq)+','+str(stop_freq)+'_out_IQ_g',out_IQ_g)
#    np.savetxt(save_dir+'ROcavity_spec_with_qubitGroundState_'+'probe_amp_0p05'+'_freq'+str(start_freq)+','+str(stop_freq)+'_out_phase_g',out_phase_g)
#    np.savetxt(save_dir+'ROcavity_spec_with_qubitExcitedState_'+'probe_amp_0p05'+'_freq'+str(start_freq)+','+str(stop_freq)+'_out_IQ_e',out_IQ_e)
#    np.savetxt(save_dir+'ROcavity_spec_with_qubitExcitedState_'+'probe_amp_0p05'+'_freq'+str(start_freq)+','+str(stop_freq)+'_out_phase_e',out_phase_e)


    
if run_spec_ROcavity_powerSweep:
    sweep_steps = 1
    start_freq = 6.8875 #7.015 #7.0161
    stop_freq = 6.89 #7.017 #7.01635
    freq = np.linspace(start_freq,stop_freq,sweep_steps)
    
    amp_start = 0.25 #0.001
    amp_stop = 0.25
    sweep_steps_2 = 1 #100    
    sweep_vals = np.linspace(amp_start,amp_stop,sweep_steps_2)
    
    out_IQ_tot = np.zeros((sweep_steps_2, sweep_steps))
    out_phase_tot = np.zeros((sweep_steps_2, sweep_steps))

    for i,val in enumerate(sweep_vals):
    
        print('current sweep val = ', val)
        print('current index = ', i)
        
        # qubit at ground state
        Fs.spectroscopy_ROcavity(num_steps=3,ssm_ge=ssm_ge,pi_ge_time=pi_ge_time,spec_amp=val,estate = 0, q = 0)
        out_IQ, out_phase = chevron.sweep_bnc_freq_cavity_spec(start_freq, stop_freq, sweep_steps,3,reps,ro_dur)

        # qubit at excited state
#        Fs.spectroscopy_ROcavity(num_steps=3,ssm_ge=ssm_ge,pi_ge_time=pi_ge_time,spec_amp=val,estate = 1, q = 1)
#        out_IQ, out_phase = chevron.sweep_bnc_freq_cavity_spec(start_freq, stop_freq, sweep_steps,3,reps,ro_dur)

        
        plt.plot(freq,out_phase,'b-')
        plt.xlabel('freq (GHz)')
        plt.ylabel('phase')
        plt.legend()
        plt.show()
        
        out_IQ_tot[i] = out_IQ
        out_phase_tot[i] = out_phase
    
        #   plt.figure(figsize=(4, 3))
    plt.imshow(out_IQ_tot, extent=[start_freq, stop_freq, amp_stop, amp_start],aspect='auto')
    plt.xlabel('freq (GHz)')
    plt.ylabel('amp (a.u.)')
    plt.show()
    
    out_phase_tot_corr = np.zeros((sweep_steps_2, sweep_steps))
    out_IQ_tot_corr = np.zeros((sweep_steps_2, sweep_steps))
    for i in range(sweep_steps_2):
        out_phase_tot_corr[i] = out_phase_tot[i] - out_phase_tot[i,-1]
        out_IQ_tot_corr[i] = out_IQ_tot[i] - out_IQ_tot[i,-1]
        #   plt.figure(figsize=(4, 3))
    plt.imshow(out_phase_tot_corr, extent=[start_freq, stop_freq, amp_stop, amp_start],aspect='auto')
    plt.xlabel('freq (GHz)')
    plt.ylabel('amp (a.u.)')
    plt.show()

    plt.imshow(out_IQ_tot_corr, extent=[start_freq, stop_freq, amp_stop, amp_start],aspect='auto')
    plt.xlabel('freq (GHz)')
    plt.ylabel('amp (a.u.)')
    plt.show()

#    plt.plot(freq,out_phase_tot_corr[20],'b-',label='amp=0.02')
#    plt.plot(freq,out_phase_tot_corr[40],'r-',label='amp=0.04')
#    plt.plot(freq,out_phase_tot_corr[60],'k-',label='amp=0.06')
#    plt.plot(freq,out_phase_tot_corr[80],'g-',label='amp=0.08')
#    plt.legend()
#    plt.xlabel('freq')
#    plt.ylabel('phase')
#    plt.show()

#    np.savetxt(save_dir+'ROcavity_spec_with_qubitGroundState_sweepProbePower_'+'probe_amp'+str(amp_start)+','+str(amp_stop)+'_freq'+str(start_freq)+','+str(stop_freq)+'_out_IQ_tot',out_IQ_tot)
#    np.savetxt(save_dir+'ROcavity_spec_with_qubitGroundState_sweepProbePower_'+'probe_amp'+str(amp_start)+','+str(amp_stop)+'_freq'+str(start_freq)+','+str(stop_freq)+'_out_phase_tot',out_phase_tot)
#
#    out_IQ_tot_e = np.zeros((sweep_steps_2, sweep_steps))
#    IQ_difference = np.zeros((sweep_steps_2, sweep_steps))
#    out_IQ_tot_e = np.loadtxt(save_dir+"ROcavity_spec_with_qubitExcitedState_sweepProbePower_probe_amp0.001,0.1_freq7.015,7.017_out_phase_tot")
#    IQ_difference = out_IQ_tot_e - out_phase_tot
#    
#    plt.imshow(IQ_difference,cmap='bwr', extent=[start_freq, stop_freq, amp_stop, amp_start],aspect='auto')
#    plt.xlabel('freq (GHz)')
#    plt.ylabel('amp (a.u.)')
#    plt.colorbar()
#    plt.show()
    
    
if run_raman_cavityPowerSweep:

    detun = -0.03 #-0.005 # -0.01 # -0.015
    ssm_cavity = ROIF2 + detun
    
    sweep_steps = 51 #81
    start_freq = 4.32  #4.2305 #4.229 #7.0161
    stop_freq = 4.32 + 0.05 #4.231 #7.01635
    freq = np.linspace(start_freq,stop_freq,sweep_steps)
    
    amp_start = 0.2#0.17 #0.0
    amp_stop = 0.9 #0.3
    sweep_steps_2 = 36
    sweep_vals = np.linspace(amp_start,amp_stop,sweep_steps_2)
    
    out_Q = np.zeros((sweep_steps_2, sweep_steps))
    out_I = np.zeros((sweep_steps_2, sweep_steps))
    pop_e_tot = np.zeros((sweep_steps_2, sweep_steps))
    sigma_z_tot = np.zeros((sweep_steps_2, sweep_steps))
    for i,val in enumerate(sweep_vals):
    
        print('current sweep val = ', val)
        print('current index = ', i)
        
        Fs.raman_ROcavity(num_steps=3,ssm_ge=ssm_ge,ge_amp = val, amp_cavity=1, ssm_cavity=ssm_cavity, rabi_time=10000, ROIF=ROIF2, q = 0)
        pop_e = copy_chevron.sweep_bnc_freq_raman(start_freq, stop_freq, sweep_steps,ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2,3,reps,ro_dur,qubit_1_thr,qubit_2_thr,wx_amps)
        
        
        pop_e_tot[i] = 1 - pop_e
        sigma_z_tot[i] = 1 - 2*pop_e

#        out_Q[i] = Q
#        out_I[i] = I
#        
        plt.plot(freq,pop_e,'b-')
        plt.xlabel('freq (GHz)')
        plt.ylabel('pop_e')
        plt.legend()
        plt.show()
        
    
        #   plt.figure(figsize=(4, 3))
    plt.imshow(pop_e_tot,  extent=[start_freq, stop_freq, amp_stop, amp_start],aspect='auto')
    plt.xlabel('qubit drive freq (GHz)')
    plt.ylabel('qubit drive amp (a.u.)')
    plt.colorbar()
    plt.show()
    
    plt.imshow(sigma_z_tot, cmap='bwr', extent=[start_freq, stop_freq, amp_stop, amp_start],aspect='auto')
    plt.xlabel('qubit drive freq (GHz)')
    plt.ylabel('qubit drive amp (a.u.)')
    plt.colorbar()
    plt.show()
    
#    np.savetxt(save_dir+'raman_spec_sweepQubitProbePower_'+'probe_amp'+str(amp_start)+','+str(amp_stop)+'_qubti_drive_freq'+str(start_freq)+','+str(stop_freq)+'_pop_e_tot'+'_032924',pop_e_tot)
#    np.savetxt(save_dir+'raman_spec_sweepCavityProbePower_'+'probe_amp'+str(amp_start)+','+str(amp_stop)+'_qubit_drive_freq'+str(start_freq)+','+str(stop_freq)+'_out_Q_tot',out_Q)

if run_raman_cavityDetuningSweep:

#    detun = -0.05-0.05*2
#    ssm_cavity = detun
    amp_cavity = 0.3
    
    detun_start = -0.02 #0.0
    detun_stop = -0.06 #0.3
    sweep_steps_2 = 9
    sweep_vals = np.linspace(detun_start,detun_stop,sweep_steps_2)

    # BNC sweep steps for qubit drive
    freq_span = 0.01
    sweep_steps = 21 #81
    freq = np.linspace(-freq_span/2,freq_span/2,sweep_steps)

  
    reps = 500
    
    out_Q = np.zeros((sweep_steps_2, sweep_steps))
    out_I = np.zeros((sweep_steps_2, sweep_steps))
    pop_e_tot = np.zeros((sweep_steps_2, sweep_steps))
    sigma_z_tot = np.zeros((sweep_steps_2, sweep_steps))
    
    for i,val in enumerate(sweep_vals):
    
        print('current sweep val = ', val)
        print('current index = ', i)

        ssm_cavity = val
        
        start_freq = 4.5 - val - freq_span/2
        stop_freq = 4.5 - val + freq_span/2

        
        Fs.raman_ROcavity(num_steps=3,ssm_ge=ssm_ge,ge_amp = 1, amp_cavity=amp_cavity, ssm_cavity=val, rabi_time=20000, q = 0)
        pop_e = chevron.sweep_bnc_freq_raman(start_freq, stop_freq, sweep_steps,3,reps,ro_dur)
        
        pop_e_tot[i] = pop_e
        sigma_z_tot[i] = 1 - 2*pop_e

#        out_Q[i] = Q
#        out_I[i] = I
#        
        plt.plot(freq,pop_e,'b-')
        plt.xlabel('freq (GHz)')
        plt.ylabel('pop_e')
        plt.legend()
        plt.show()
        
    
        #   plt.figure(figsize=(4, 3))
    plt.imshow(pop_e_tot, cmap='bwr',  extent=[-freq_span/2, freq_span/2, detun_stop, detun_start],aspect='auto')
    plt.xlabel('qubit drive freq detun (GHz)')
    plt.ylabel('cavity detun (GHz)')
    plt.colorbar()
    plt.show()
    
    plt.imshow(sigma_z_tot, cmap='bwr', extent=[-freq_span/2, freq_span/2, detun_stop, detun_start],aspect='auto')
    plt.xlabel('qubit drive freq detun (GHz)')
    plt.ylabel('cavity detun (GHz)')
    plt.colorbar()
    plt.show()
    
#    np.savetxt(save_dir+'raman_spec_sweepCavityProbeDetun_'+'detun_'+str(detun_start)+','+str(detun_stop)+'_qubti_drive_freq_span'+str(freq_span)+'_pop_e_tot',pop_e_tot)
#    np.savetxt(save_dir+'raman_spec_sweepCavityProbePower_'+'probe_amp'+str(amp_start)+','+str(amp_stop)+'_qubit_drive_freq'+str(start_freq)+','+str(stop_freq)+'_out_Q_tot',out_Q)
    
#    freq_raman_bnc = np.zeros(sweep_steps_2)
#    freq_add_detun = np.zeros(sweep_steps_2)
#    for i, val in enumerate(sweep_vals):
#        
#
#        freq_add_detun[i] = freq[np.argmax(pop_e_tot[i])]
#        freq_raman_bnc[i] = 4.5 - val + freq_add_detun[i]
#        
#    plt.plot(sweep_vals, freq_add_detun,'bs')
#    plt.xlabel('RO cavity detun (GHz)')
#    plt.ylabel('added detun (GHz)')
#    plt.ylim([-0.002, 0.001])
#    plt.show()
#    
#    plt.plot(sweep_vals, freq_raman_bnc,'b-s')
#    plt.xlabel('RO cavity detun (GHz)')
#    plt.ylabel('qubit BNC freq (GHz)')
#    plt.show()
    
if run_tomo_calib:
    num_steps = 51
    sweep_time =200 #ns
    sweep_steps = 19
    sweep_vals = np.linspace(90,99,sweep_steps)
    amplitude = []
    out_tomo = np.zeros( (sweep_steps, num_steps_i))
    for i,val in enumerate(sweep_vals):
        print(val)
#        Fs.full_tomo(off=0,num_steps=num_steps_i,sweep_time=sweep_time,ssm_ge=ssm_ge,pi_ge=pi_ge_time,tomo=1,amp_tomo=val,ph=0,q=which_qubit) #y
        Fs.full_tomo(off=0,num_steps=num_steps,sweep_time=sweep_time,ssm_ge=ssm_ge,pi_ge=pi_ge_time,tomo=0.98744,amp_tomo=1.1231,ph=val,q=which_qubit) #x
        Fs.loading(num_steps)
        wx_programs.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
        daq_params, t_histo, p,a,b,n_readout = daq_programs.run_daq_cluster_threshold(model,Ax=Ax,By=By,num_patterns=num_steps,
                                                                                         num_records_per_pattern=reps,authr=0,fg=3,ro_dur=ro_dur,IQangle=IQ_angle)
        p_readout = p
        p_readout= np.matmul(scale_matrix,p)
        plt.plot(p_readout[0],label='|g>');plt.plot(p_readout[1],label='|e>');plt.plot(p_readout[2],label='|f>')
        plt.legend();plt.title('data scaled with matrix');plt.show()
        y=p_readout[0]
        out_tomo[i] = y
        times = np.linspace(0,sweep_time/1000,num_steps_i)
        popt,peer, y_vals, _ = analysis.fit_sine_decay(times,y,guess_vals=[10,0.5,0.14,0.5,0.5])
        amplitude.append(popt[2])
        np.savetxt(save_dir+'sweep_phase90,99_19steps_200nsrabi_51time_1000avg_scaledmatrix_xtomo',out_tomo)  
        np.savetxt(save_dir+'sweep_phase90,99_19steps_200nsrabi_51time_1000avg_scaledmatrix_amp',amplitude)
        
    plt.imshow(out_tomo)
    plt.show()
    plt.plot(sweep_vals,np.abs(amplitude))   
    plt.show()
    
if run_qubit_spectroscopy_vs_cavity_drive:
    
    detun = -0.015 #-0.015
    ssb_cavity = ROIF2 + detun #-0.249  #+ 0.0065
    amp_cavity = 0.6
    
    reps = 500

    num_steps = 41
    f1 = -0.1538 - 0.01 
    f2 = -0.1538 + 0.01 #-0.05
    freq = np.linspace(f1,f2,num_steps)

    # sweep off-resonance cavity drive 
    amp_cavity_start = 0
    amp_cavity_stop = 1
    sweep_steps = 11    
    sweep_vals = np.linspace(amp_cavity_start,amp_cavity_stop,sweep_steps)
#
    T1_time = []
    out_Q = np.zeros((sweep_steps, num_steps))
    
    spec_freq = []; spec_bandwidth = []
    
    for i,val in enumerate(sweep_vals):
    
        print('current sweep val = ', val)
        print('current index = ', i)
        Fs.spectroscopy_ge_with_cavityDrive(num_steps,ssm_start=f1,ssm_stop=f2,spec_amp=0.01, amp_cavity=val, ssm_cavity=ssb_cavity, ROIF=ROIF2, q = which_qubit)
        
        wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave

#        # No thresholding
#        daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_pi,counts_pi,rec_readout,rec_avg_all,rec_all = daq_programs.run_daq(num_steps,reps,ro_dur,IQangle=100)
#        
#        plt.plot(Q);plt.show();plt.plot(I);plt.show()
#        out_Q[i] = I
#        ind = np.argmax(out_Q[i])
#        spec_freq.append(freq[ind])


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
        
        out_Q[i] = P_Q2
        ind = np.argmax(out_Q[i])
        spec_freq.append(freq[ind])

    
    power = np.square(sweep_vals)
    
    plt.plot(power,spec_freq,'b-s')
    plt.xlabel('power(a.u.)')
    plt.ylabel('qubit ssm (GHz)')
    plt.show()
    
    #   plt.figure(figsize=(4, 3))
    plt.imshow(out_Q, extent=[f1, f2, amp_cavity_stop**2, amp_cavity_start**2],aspect='auto')
    plt.xlabel('ssm_ge (GHz)')
    plt.ylabel('power (a.u.)')
    plt.show()

#    np.savetxt(save_dir+'spectroscopy_under_cavityDrive_'+'detun_'+str(detun)+','+str(amp_cavity_start)+','+str(amp_cavity_stop)+'ssm_range'+str(f1)+'to'+'f2',out_Q)
    
    
if run_ramsey_coax_drive:
    ssb_coax = 0.02 #-0.25
    ssb_start = ssb_coax - 0.05
    ssb_stop = ssb_coax + 0.05
    sweep_steps = 11 # 41
    sweep_vals = np.linspace(ssb_start,ssb_stop,sweep_steps)
    num_steps = 81
    t2_time = 5000
    amp_coax = 1
    T2_time = []
    freq_fit = []
    out_I = np.zeros( (sweep_steps, num_steps))
    for i,val in enumerate(sweep_vals):
        print('current sweep val = ', val)
        print('current index = ', i)
        Fs.ramsey_coax_drive(num_steps = num_steps,t1_time = t2_time,pi_echo_coef=0,osc=4,ssm_ge=ssm_ge,ssm_coax=val,amp_coax=amp_coax,pi_ge=pi_ge_time,q=which_qubit)
        
        wx_programs.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
        daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_pi,counts_pi,rec_readout,rec_avg_all,rec_all = daq_programs.run_daq(num_steps,reps,ro_dur,IQangle=IQ_angle)
        plt.plot(Q);plt.show();plt.plot(I);plt.show()
        out_I[i] = I
        
        times = np.linspace(0,t2_time/1000,num_steps)
        T2_fit_vals,_,_,_ = analysis.fit_sine_decay(times,Q,guess_vals=[0.8,0.6,2,1,50])
        T2 = 1/T2_fit_vals[1]
        freq = T2_fit_vals[0]
        T2_time.append(T2)
        freq_fit.append(freq)
        print("T2* = {} \u03BCs".format(T2))
#        daq_params, t_histo, p,a,b,n_readout = daq_programs.run_daq_cluster_threshold(model,Ax=Ax,By=By,num_patterns=num_steps,
#                                                                                         num_records_per_pattern=reps,authr=0,fg=3,ro_dur=ro_dur,IQangle=IQ_angle)
#        p_readout = p
#        p_readout= np.matmul(scale_matrix,p)
#        plt.plot(p_readout[0],label='|g>');plt.plot(p_readout[1],label='|e>');plt.plot(p_readout[2],label='|f>')
#        plt.legend();plt.title('data scaled with matrix');plt.show()
#        y=p_readout[1]
#        plt.plot(y)

#        np.savetxt(save_dir+'ramsey_coax_drive_'+str(ssb_start)+','+str(ssb_stop)+'_coax_amp'+str(amp_coax)+'t2sweeptime'+str(t2_time)+'_fullsequence',out_I)
#        np.savetxt(save_dir+'ramsey_coax_drive_'+str(ssb_start)+','+str(ssb_stop)+'_coax_amp'+str(amp_coax)+'t2sweeptime'+str(t2_time)+'_T2_times',T2_time)
#        np.savetxt(save_dir+'ramsey_coax_drive_'+str(ssb_start)+','+str(ssb_stop)+'_coax_amp'+str(amp_coax)+'t2sweeptime'+str(t2_time)+'_freq_fit',freq_fit)
        
    plt.plot(sweep_vals,np.abs(freq_fit),label='freq');plt.plot(sweep_vals,T2_time,label='T2');
    plt.xlabel('ssb freq (GHz)')
#    plt.ylim([0,5])
    plt.legend();plt.show()
    
#    plt.figure(figsize=(4, 3))
    plt.imshow(out_I, extent=[0, t2_time/1000, ssb_stop, ssb_start],aspect='auto')
    plt.xlabel('time (us)')
    plt.ylabel('ssb freq (GHz)')
    plt.show()


if run_T1_ge_raman_cooling:
    
    detun = -0.03 #-0.01 #-0.005 # -0.015 
    ssb_coax = ROIF2 + detun #-0.249  #+ 0.0065
    amp_coax = 1 
    ssb_Q2 = ssm_ge + detun  # excitation -detun; decay + detun
    amp_Q2 =  0.5 #0.5 #0.34 #0.17 # 0.54 

    num_steps = 21
    t1_time = 20000
    sweep_steps = 1 # 41

    # sweep detuning
    detun_add_start = 0.005 #-0.00475 #0.001 #-0.01# 0.026 #-0.02 #-0.005 #-0.003 # -0.1 #-0.125 #-0.2 # 0.15 #ssb_coax - 0.03#0.094 #ssb_coax - 0.03 
    detun_add_stop = 0.005#0.030 #0.094 #ssb_coax - 0.02 #0.0066
    sweep_vals = np.linspace(detun_add_start,detun_add_stop,sweep_steps)
#
    reps = 1000

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
    
    for i,val in enumerate(sweep_vals):
    
        print('current sweep val = ', val)
        print('current index = ', i)
        Fs.T1_ge_tunable_dissipation_gaussian_pulse(num_steps= num_steps, sweeptime = t1_time,ssm_ge=ssm_ge,pi_ge_time=pi_ge_time,ssm_coax=ssb_coax, amp_coax=amp_coax, ssm_Q2=ssb_Q2 + val, amp_Q2=amp_Q2, ROIF=ROIF2, q=which_qubit, ifload = 1)
        
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
        
        p2_readout= np.matmul(scale_matrix2,prob_vs_pats_2)    
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
    
#    np.savetxt(save_dir+'t1_decay_heating_20us_detun-30MHz_ampCavity1.0_ampQubit0.5_qubitDetunAdd-4.75MHz_pop_e_022624',pop_e)
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


