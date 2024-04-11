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

#sys.path.append(r"C:\Users\crow104\Documents\Python Scripts\sequence_generator")
from generator_nonHtrial import *
#from generator import *
import wx_programs
import os
#import expt_parameters
#import seq_experiments
#import seq_programs
import daq_programs
import analysis
import math


import pyvisa
import daq_programs as dp
from tkinter import Tk
import tkinter as tki
from tkinter.filedialog import askopenfilename, asksaveasfilename, asksaveasfile
import scipy.io as io
#import data_analysis as da
import time
from IPython import get_ipython
#import data_analysis as da
from scipy import optimize as opt
from matplotlib.patches import Rectangle
#import x_readout_calibration as xr
import warnings
#import vaunix_control as vc
from scipy.optimize import fsolve


pi = np.pi
import matplotlib.pyplot as plt



def daqnon(num_patterns=51,pi_f=1,n_g=3,au=0,m_r=1,n_p=2,num_avgs=1,prev_threshold=[0,0],num_records_per_patterns=1000,times=1000,fun="dec"):
    popt1=[]
    peer1=[]
    y_vals1=[]
    popt2=[]
    peer2=[]
    y_vals2=[] 
    popt1p=[]
    peer1p=[]
    y_vals1p=[] 
    popt2p=[]
    peer2p=[]
    y_vals2p=[] 
    loading(num_steps = pi_f*num_patterns)
    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[-.052,+.053, 0.057, -0.061])
    for k in range(num_avgs):
        
        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
            num_patterns=pi_f*num_patterns, num_records_per_pattern=num_records_per_patterns,authr=au,fg=n_g)
        

        if k is 0:
            p_readout = p
#                   
        else:
            p_readout += p
#
    p_post = analysis.p_readout_postselected(p_readout)
#    p_post=analysis.p_readout_postselected_pief(p_readout)
    x=times
    y=p_readout[m_r]
    yp=p_post[n_p]     
    plt.plot(x,y)
    plt.plot(x,yp)
    if fun=="sin":
        popt1,peer1, y_vals1, _ = analysis.fit_sine_decay(x,y,guess_vals=None) 
    if fun=="dec":    
        popt2,peer2, y_vals2, _ = analysis.fit_exp_decay(x,y,guess_vals=None)
    if n_g == 3:
        if fun=="sin":
            popt1p,peer1p, y_vals1p, _ = analysis.fit_sine_decay(x,yp,guess_vals=None) 

    return x,y,yp,popt1,popt2,popt1p,popt2p,prev_threshold



def es_transport(ssm_ge = 0.3865,ssm_ef = 0.0912): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    rabi_time = 4000
    pi_ge=34
    pi_ef=28 
    pi_hf=26

#    ssm_ge = 0.3885
#    ssm_ef = 0.09175
    ssm_hf = 0.205

#    p_pi_ge = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#
#    p_pi_ef = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
#    p_pi_ef.phase = 180
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
#
    p_pi_ge = Pulse(start=16995-pi_ge/2, duration=-pi_ge/2, amplitude=.5, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    
    pt_pi_ef_r = Pulse(start=16995-pi_ge/2, duration=0, amplitude=.5, ssm_freq=ssm_ge, phase=180,phase_ini=0, t_loop=rabi_time, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
    pt_pi_ef_r.phase = 270
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    pt_pi_ef = Pulse(start=16995-pi_ge/2, duration=0, amplitude=0.5, ssm_freq=ssm_ge, phase=270,phase_ini=1*np.pi/2, t_loop=rabi_time, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef)
    pt_pi_ef.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef)
    
    p_pi_ge_r = Pulse(start=16995, duration=-pi_ge/2, amplitude=.5, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge_r)
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,5000:7000], aspect='auto', extent=[5000,7000,200,0])
#        plt.plot(channel1_ch[50,:],'b--o')
 
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\test"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom



def es_drive(pi_ge=34): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    rabi_time = 2000
#    pi_ge=34
    pi_ef=28 
    pi_hf=26

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.205

#    p_pi_ge = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#
#    p_pi_ef = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
#    p_pi_ef.phase = 180
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)

    p_pi_ge = Pulse(start=16995-pi_ge/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase =90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#
    p_pi_ge = Pulse(start=16995, duration=-pi_ge/2, amplitude=.5, ssm_freq=ssm_ge, phase=270)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    
    pt_pi_ef_r = Pulse(start=16995, duration=0, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
    pt_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)


    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,5000:7000], aspect='auto', extent=[5000,7000,200,0])
#        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\test"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def es_drive_jump(amp_ef=0.3): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 101
#    detunlin = - 0.1
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
# note: amp_ef 0.3 for eigenstate with less dissipation, phase 270, 360
    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    rabi_time = 200
    pi_ge=34
    pi_ef=28 
    pi_hf=26
    
    ssm_ef_detun = -0.01
#    amp_ef=0.2
    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.205

#    p_pi_ge = Pulse(start=16995-pi_ef, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
#    p_pi_ge.phase =90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    ###
#    p_pi_ge = Pulse(start=16995-32-pi_ef, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#    p_pi_ge.phase =90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#
#    p_pi_ef = Pulse(start=16995-pi_ef, duration=-32, amplitude=amp_ef, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
#    p_pi_ef.phase = 180
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
#    
#    pt_pi_ef_r = Pulse(start=16995-pi_ef, duration=0, amplitude=.5, ssm_freq=ssm_ef+ssm_ef_detun, phase=0, detunlinear = ssm_ef_detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
#    pt_pi_ef_r.phase = 90
#    the_seq.add_sweep(channel=4,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
##
#    p_pi_ef = Pulse(start=16995, duration=-pi_ef, amplitude=0.25, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='phase_linear_detun', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
#    p_pi_ef.phase =90
#    the_seq.add_sweep(channel=2,  sweep_name='phase_linear_detun', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)  
# 
#    p_pi_ef = Pulse(start=16995, duration=-pi_ef, amplitude=0.25, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
#    p_pi_ef.phase =90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)     
    
    ###
    #### ###
    p_pi_ge = Pulse(start=16995-32-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase =90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_ef/2, duration=-32, amplitude=amp_ef, ssm_freq=ssm_ef, phase=270, clock_freq=ssm_ef-ssm_ef_detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
    p_pi_ef.phase = 360
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
    
    pt_pi_ef_r = Pulse(start=16995-pi_ef/2, duration=0, amplitude=.5, ssm_freq=ssm_ef+ssm_ef_detun, phase=0, clock_freq=ssm_ef-ssm_ef_detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
    pt_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
    
    p_pi_ef = Pulse(start=16995, duration=-pi_ef/2, amplitude=0.5, ssm_freq=ssm_ef, phase=270, clock_freq=ssm_ef-ssm_ef_detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ef)
    p_pi_ef.phase = 360
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ef)  
    #### ###
    
#    p_pi_ef = Pulse(start=16995, duration=-pi_ef, amplitude=0.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ef) 
#    p_pi_ef.phase =90
#    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ef)
#    
#    pt_pi_ef_r = Pulse(start=16995+pi_ef/2, duration=pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=3, sweep_name='phase_linear_detun', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
#    pt_pi_ef_r.phase = 90
#    the_seq.add_sweep(channel=4,  sweep_name='phase_linear_detun', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)


    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        plt.plot(channel1_ch[10,:],'b--o')
#        plt.plot(channel2_ch[2,:],'y--*')
#        plt.plot(channel3_ch[10,:],'r-s')
#        plt.xlim(6910,7000)
#        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\test"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def es_drive_ef(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 2000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=34
#
##    ssm_ge = 0.3885
##    ssm_ef = 0.0917
#    ssm_hf = 0.243

#    p_pi_ge = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#
#    p_pi_ef = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
#    p_pi_ef.phase = 180
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)

    p_pi_ge = Pulse(start=16995-3*pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase =90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)

    p_pi_ge = Pulse(start=16995-pi_ef/2, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase =90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#
    p_pi_ge = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    
    pt_pi_ef_r = Pulse(start=16995, duration=0, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
    pt_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)


    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,5000:7000], aspect='auto', extent=[5000,7000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom


def phase_meas_esdrive(rabi_time = 100,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
    pi_ge=16
    pi_ef=28 
    pi_hf=26

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.205
#    off_set=0


#
    p_pi_ge = Pulse(start=16995-2*pi_ge-1*pi_ef-rabi_time, duration=-2*pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ef = Pulse(start=16995-2*pi_ge-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    

    p2_pi_ge = Pulse(start=16995-pi_ge-pi_ef/2-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=270)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ge-pi_ef/2, duration=-rabi_time, amplitude=.3, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)

    
    p_pi_ge_r = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge_r)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def phase_meas_esdrive_diff_ang(rabi_time = 100,off_set=0,ang=160): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
    pi_ge=16
    pi_ef=28 
    pi_hf=26

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.205
#    off_set=0
    
    


#
    p_pi_ge = Pulse(start=16995-2*pi_ge-1*pi_ef-rabi_time, duration=-2*pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ef = Pulse(start=16995-2*pi_ge-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    

    p2_pi_ge = Pulse(start=16995-pi_ge-pi_ef/2-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=270+ang)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 0+ang
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ge-pi_ef/2, duration=-rabi_time, amplitude=.3, ssm_freq=ssm_ge, phase=0+ang)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90+ang
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)

    
    p_pi_ge_r = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=90+ang)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 180+ang
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge_r)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def xdrive(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=270,ort=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270


#
    p_pi_ge = Pulse(start=16995-3*pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90+ort
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_ef/2, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90+ort
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
    


    p2_pi_ge = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+esph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph+ort
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995, duration=0, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width',start=0, stop=-rabi_time,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase =180+ort
    the_seq.add_sweep(channel=2,  sweep_name='width',start=0, stop=-rabi_time,initial_pulse=pt_pi_ge_r)



    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def test_phase_meas_esdrive_hf(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270


#
    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-pi_ef-rabi_time, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-rabi_time, duration=-pi_hf/2, amplitude=amp, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    

    p2_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+esph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2, duration=-rabi_time, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)

    
    p_pi_ef_r = Pulse(start=16995-pi_hf/2-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=180+esph+eph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 270+esph+eph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)


    rabi_ef = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=amp, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)
    rabi_ef.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def phase_meas_esdrive_hf(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270


#
    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-pi_ef-rabi_time, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-rabi_time, duration=-pi_hf/2, amplitude=amp, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    

    p2_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+esph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2, duration=-rabi_time, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)

    
    p_pi_ef_r = Pulse(start=16995-pi_hf/2-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=180+esph+eph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 270+esph+eph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)


    rabi_ef = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=amp, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)
    rabi_ef.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def pi2prep(): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    rabi_time = 0
    pi_ge=34
    pi_ef=28 
    pi_hf=26

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.205

    p_pi_ge = Pulse(start=16995-pi_ef, duration=-pi_ge, amplitude=0.52, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ge_r = Pulse(start=16995, duration=-pi_ef, amplitude=0, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='amplitude' , start=0, stop=1 ,initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='amplitude' , start=0, stop=1,initial_pulse=p_pi_ge_r)
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none' ,initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,6900:7000], aspect='auto', extent=[6900,7000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\test"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def loading(num_steps = 51):
      

    file_length = 30000#50000 #30000#18000
#    num_steps = 51*3
    the_seq = Sequence(file_length, num_steps)      
#    write_dir = r"C:\Data\2019\encircling\rabi_ef"
#    write_dir = r"C:\Data\2019\encircling\phase_measurement"
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)

    the_seq.load_sequence('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
#    the_seq.load_sequence('128.252.0.100', base_name='foo', file_path=write_dir, num_offset=0,amp34=.4)    
##END geom
 
def loading_jump():
      

    file_length = 18000
    num_steps = 2*301
    the_seq = Sequence(file_length, num_steps)      
    write_dir = r"C:\Data\2019\encircling\test"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def phase_meas_transport(rabi_time = 100,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
    pi_ge=16
    pi_ef=28 
    pi_hf=26

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.205
#    off_set=0


#
    p_pi_ge = Pulse(start=16995-2*pi_ge-1*pi_ef-rabi_time, duration=-2*pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ef = Pulse(start=16995-2*pi_ge-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    

    p2_pi_ge = Pulse(start=16995-pi_ge-pi_ef/2-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ge-pi_ef/2, duration=-rabi_time, amplitude=1, ssm_freq=ssm_ge, phase=0,phase_ini=0, t_loop=rabi_time, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    pt_pi_ge = Pulse(start=16995-pi_ge-pi_ef/2, duration=-rabi_time, amplitude=1, ssm_freq=ssm_ge, phase=90,phase_ini=1*np.pi/2, t_loop=rabi_time, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge)
    pt_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge)
    
    p_pi_ge_r = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=270)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge_r)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def phase_meas_transport_hf(rabi_time = 100,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
    pi_ge=34
    pi_ef=28 
    pi_hf=40
    sphase = 90

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.243
#    off_set=0


#
    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ge = Pulse(start=16995-pi_hf-2*pi_ef-rabi_time, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-rabi_time, duration=-pi_hf/2, amplitude=.5, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    
    p_pi_ef = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=270+sphase)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 0+sphase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    
    pt_pi_ge_r = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2, duration=-rabi_time, amplitude=.1, ssm_freq=ssm_ef, phase=0,phase_ini=0, t_loop=rabi_time, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=pt_pi_ge_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    pt_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2, duration=-rabi_time, amplitude=.1, ssm_freq=ssm_ef, phase=90,phase_ini=1*np.pi/2, t_loop=rabi_time, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=pt_pi_ge)
    pt_pi_ge.phase = 180
    the_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=pt_pi_ge)
    
    p_pi_ge_r = Pulse(start=16995-pi_hf/2-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90+sphase)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 180+sphase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge_r)

    p_pi_ef_r = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=.5, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)

    p_pi_ge = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6750:7000], aspect='auto', extent=[6750,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def phase_meas_partial_transport_hf(rabi_time = 100,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
    tloop=1500
    pi_ge=34
    pi_ef=28 
    pi_hf=40
    sphase = 0
    par=rabi_time/tloop
    

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.243
#    off_set=0


#
    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ge = Pulse(start=16995-pi_hf-2*pi_ef-rabi_time, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-rabi_time, duration=-pi_hf/2, amplitude=.5, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    
    p_pi_ef = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=270+sphase)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 0+sphase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    
    pt_pi_ge_r = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2, duration=-rabi_time, amplitude=.5, ssm_freq=ssm_ef, phase=0,phase_ini=0, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=pt_pi_ge_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    pt_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2, duration=-rabi_time, amplitude=.5, ssm_freq=ssm_ef, phase=90,phase_ini=1*np.pi/2, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=pt_pi_ge)
    pt_pi_ge.phase = 180
    the_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=pt_pi_ge)
    
    p_pi_ge_r = Pulse(start=16995-pi_hf/2-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90+sphase+360/par)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 180+sphase+360/par
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge_r)

    p_pi_ef_r = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=.5, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)

    p_pi_ge = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6750:7000], aspect='auto', extent=[6750,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def phase_meas_partial_transport(rabi_time = 100,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
    tloop=1500
#    phase_ini=np.pi/2
#    rabi_time = 4000
    par=rabi_time/tloop
    pi_ge=16
    pi_ef=28 
    pi_hf=26

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.2
#    off_set=0


#
    p_pi_ge = Pulse(start=16995-2*pi_ge-1*pi_ef-rabi_time, duration=-2*pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ef = Pulse(start=16995-2*pi_ge-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    

    p2_pi_ge = Pulse(start=16995-pi_ge-pi_ef/2-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ge-pi_ef/2, duration=-rabi_time, amplitude=.5, ssm_freq=ssm_ge, phase=0,phase_ini=0, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    pt_pi_ge = Pulse(start=16995-pi_ge-pi_ef/2, duration=-rabi_time, amplitude=.5, ssm_freq=ssm_ge, phase=90,phase_ini=1*np.pi/2, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge)
    pt_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge)
    
    p_pi_ge_r = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=270 + 360*1/par)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 0 + 360*1/par
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge_r)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def xtransport_ef_tVar(ssm_ge=.3865,ssm_ef=.0912,off_set=0,ph=0,ph2=0,amplz=.5,epamp=.15,rabi_time=500,pi_ge=34,pi_ef=28): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    tloop=2000
#    rabi_time = 2000
    tloop = rabi_time
    par=rabi_time/tloop
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=34
#    par=4/3
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.243
#    off_set=0


    p_pi_ge = Pulse(start=16995-pi_ef-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none' ,initial_pulse=p_pi_ge)

    

    p2_pi_ge = Pulse(start=16995-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 180+ph
    the_seq.add_sweep(channel=2,  sweep_name='none' ,initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2, duration=-rabi_time, amplitude=epamp, ssm_freq=ssm_ef, phase=0,phase_ini=0, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='none' ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='none' ,initial_pulse=pt_pi_ge_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    pt_pi_ge = Pulse(start=16995-pi_ef/2, duration=-rabi_time, amplitude=epamp, ssm_freq=ssm_ef, phase=90,phase_ini=1*np.pi/2, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='none' ,initial_pulse=pt_pi_ge)
    pt_pi_ge.phase = 180
    the_seq.add_sweep(channel=4,  sweep_name='none' ,initial_pulse=pt_pi_ge)
    
    p_pi_ge_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=amplz, ssm_freq=ssm_ef, phase=90+360*1/par+ph2+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 180+360*1/par+ph2+ph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge_r)


    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,5500:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        plt.imshow(channel4_ch[0:200,5500:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def xtransport_ef_tVar_z(ssm_ge=.3865,ssm_ef=.0912,off_set=0,ph=0,ph2=0,amplz=.5,epamp=.15,rabi_time=500,pi_ge=34,pi_ef=28): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    tloop=2000
#    rabi_time = 2000
    tloop = rabi_time
    par=rabi_time/tloop
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=34
#    par=4/3
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.243
#    off_set=0


    p_pi_ge = Pulse(start=16995-pi_ef/2-rabi_time-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none' ,initial_pulse=p_pi_ge)

    p2_pi_ge = Pulse(start=16995-rabi_time-pi_ef/2, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 180+ph
    the_seq.add_sweep(channel=2,  sweep_name='none' ,initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2, duration=-rabi_time, amplitude=epamp, ssm_freq=ssm_ef, phase=0,phase_ini=0, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='none' ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='none' ,initial_pulse=pt_pi_ge_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    pt_pi_ge = Pulse(start=16995-pi_ef/2, duration=-rabi_time, amplitude=epamp, ssm_freq=ssm_ef, phase=90,phase_ini=1*np.pi/2, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='none' ,initial_pulse=pt_pi_ge)
    pt_pi_ge.phase = 180
    the_seq.add_sweep(channel=4,  sweep_name='none' ,initial_pulse=pt_pi_ge)
    

    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.imshow(channel4_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
##        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def xtransport_ef_z(ssm_ge=.3865,ssm_ef=.0912,off_set=0,ph=0,ph2=0,amplz=.5,epamp=.15,rabi_time = 2000): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    tloop=rabi_time
#    rabi_time = 2000
    par=rabi_time/tloop
    pi_ge=34
    pi_ef=28 
    pi_hf=34
#    par=4/3
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.243
    off_set=0


    p_pi_ge = Pulse(start=16995, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)

    

    p2_pi_ge = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 180+ph
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995, duration=0, amplitude=epamp, ssm_freq=ssm_ef, phase=0,phase_ini=0, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    pt_pi_ge = Pulse(start=16995, duration=0, amplitude=epamp, ssm_freq=ssm_ef, phase=90,phase_ini=1*np.pi/2, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge)
    pt_pi_ge.phase = 180
    the_seq.add_sweep(channel=4,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge)
    



    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def xtransport_ef(ssm_ge=.3865,ssm_ef=.0912,off_set=0,ph=0,ph2=0,amplz=.5,epamp=.15,rabi_time = 2000): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    tloop=rabi_time 
#    rabi_time = 2000
    par=rabi_time/tloop
    pi_ge=34
    pi_ef=28 
    pi_hf=34
#    par=4/3
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.243
#    off_set=0


    p_pi_ge = Pulse(start=16995-pi_ef, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    
#    p_pi_ge = Pulse(start=16995-pi_ef, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    
#    p_pi_ef = Pulse(start=16995-2*pi_ge-pi_ef/2, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ef)
#    p_pi_ef.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ef)
    

    p2_pi_ge = Pulse(start=16995-pi_ef/2, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 180+ph
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2, duration=0, amplitude=epamp, ssm_freq=ssm_ef, phase=0,phase_ini=0, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    pt_pi_ge = Pulse(start=16995-pi_ef/2, duration=0, amplitude=epamp, ssm_freq=ssm_ef, phase=90,phase_ini=1*np.pi/2, t_loop=tloop, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge)
    pt_pi_ge.phase = 180
    the_seq.add_sweep(channel=4,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge)
    
    p_pi_ge_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=amplz, ssm_freq=ssm_ef, phase=270+360*1/par+ph2+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 0+360*1/par+ph2+ph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge_r)

#    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',start=0, stop=360,initial_pulse=p_pi_ef_r)
#    p_pi_ef_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',start=0, stop=360,initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom


def encirclingLEP(rabi_time = 2000,pi_ge=32,pi_ef=28,ssm_ge = 0.3865,ssm_ef = 0.912,off_set=0,num_steps = 51,phase_tomo=0,amp_tomo=0.5,pi_time=0,pi2_time=0,min_amp=.01,detun=0.005,amp_enc=.25,ph=0,ph_ini=0,coef_in=1,amp_pief=0.5): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 201
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    tloop=rabi_time
#    rabi_time = 2000
    par=rabi_time/tloop
#    pi_time=0
#    pi_ge=31
#    pi_ef=29 
#    pi_hf=26
##    par=4/
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.205
#    off_set=0


#   encircling a Type-II LEP in the parameter space of (J, \Delta) of g-e submanifold
#   initial state: 1/2(|g>+|e>)+1/sqrt(2)|f>, or 1/2(|g>-|e>)+1/sqrt(2)|f>
#   readout: complex number \rho_fg/\rho_fe
#
    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-pi_ge-pi_ge/2-pi_ef/2, duration= -pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)

    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-pi_ge-pi_ge/2, duration= -pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    
    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-pi_ge, duration= -pi_ge/2, amplitude=.5*coef_in, ssm_freq=ssm_ge, phase=270+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 360+ph
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)

#    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_ge, duration=0, amplitude=amp_enc, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
#    pt_pi_ge_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
#    
#    p2_pi_ge = Pulse(start=16995-pi_time-pi_ef/2, duration=-pi_ge, amplitude=amp_tomo, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_ge, duration=0, amplitude=amp_enc, ssm_freq=ssm_ge, phase=0,phase_ini=ph_ini, t_loop=tloop, ff=-1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=16995-pi_time-pi_ef/2, duration=-pi2_time, amplitude=amp_tomo, ssm_freq=ssm_ge, phase=0,phase_ini=ph_ini, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)


    p_pi_ge = Pulse(start=16995-pi_time, duration= -pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90+phase_tomo)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180+phase_tomo
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ge)
    
    p2_pi_ef = Pulse(start=16995, duration=-pi_time, amplitude=amp_pief, ssm_freq=ssm_ef, phase=0,phase_ini=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
    p2_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)

    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel1_ch[0:200,16800:17000], aspect='auto', extent=[16800,17000,200,0])
#        plt.plot(channel1_ch[50,16600:17000],'b--o')   
#        plt.plot(channel2_ch[50,16000:17000],'r--o')   
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"

    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def tomography_calibration(amp = 0.5, prep_phase=0, mixer_orth=0, ro_phase=0, pi2_ef_on_off=1, pi_ef_on_off=0,off_set=2): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 8000
    num_steps = 1
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
    pi_ge= 34 #35 #34
    pi_ef= 30 #30 #14
    pi2_ef = 15 #14
#    prep_phase = 0
#    ro_phase = 90
    dur_pi_ef = pi_ef*pi_ef_on_off
    dur_pi2_ef = pi2_ef*pi2_ef_on_off
#    rabi_time = 50
    ssm_ge = 0.3813 #0.1742 #0.3870 #0.3885
    ssm_ef = 0.0862 #-0.1213 #0.0912 #0.09137 #0.0917
    
    p_pi_ge = Pulse(start=6995-dur_pi_ef-dur_pi2_ef-pi2_ef-pi_ef*0, duration= -pi_ge, amplitude=amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90+mixer_orth
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=6995-dur_pi_ef-dur_pi2_ef-pi2_ef, duration= -pi_ef*0, amplitude=amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90+mixer_orth
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ef)
    
#    p_pi_ef = Pulse(start=6995-dur_pi_ef-dur_pi2_ef, duration= -pi_ef, amplitude=amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',  initial_pulse=p_pi_ef)
#    p_pi_ef.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ef)
    
#    p_pi_ef = Pulse(start=6995-rabi_time-pi_ef, duration= -pi2_ef, amplitude=0.5, ssm_freq=ssm_ef, phase=180+phases) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='phase', start=0, stop=360, initial_pulse=p_pi_ef)
#    p_pi_ef.phase = 270+phases
#    the_seq.add_sweep(channel=2,  sweep_name='phase', start=0, stop= 360, initial_pulse=p_pi_ef)

    p_pi_ef = Pulse(start=6995-dur_pi_ef-dur_pi2_ef, duration= -pi2_ef, amplitude=amp, ssm_freq=ssm_ef, phase=0+prep_phase) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90+prep_phase
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ef)
    
    rabi_ef = Pulse(start=6995-dur_pi_ef, duration=-dur_pi2_ef, amplitude=amp, ssm_freq=ssm_ef, phase=0+ro_phase) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=rabi_ef)
    rabi_ef.phase = 90+ro_phase+mixer_orth
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=rabi_ef)

    p_pi_ef = Pulse(start=6995, duration= -dur_pi_ef, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90+mixer_orth
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ef)


#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)

    #main readout
    main_pulse = Pulse(start = 7000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6840:7000], aspect='auto', extent=[6840,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True, pr=0)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
## End geom

def rabi_ge(num_steps=51,sweep_time=200,ssm_ge=-0.15,ROIF1=0,ROIF2 = 0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = 0.251 #ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 0.261 #1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
   
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    #HOMO
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp )
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    #HET   
        #Q1 Readout
    #if q == 1:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
    
    

    #Q2 Readout
    #if q == 0:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
   
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom
    
#def rabi_ge_withRamanDrive(num_steps=51,sweep_time=10000,ssm_ge=-0.15,ROIF1=0,ROIF2=0,ssm_coax=-0.25, amp_coax=1, ssm_Q2=0.5, amp_Q2=0.2, q=0,ifsideband=1): #this is pulsed readout to ring up and ring down cavity dfor e state
#    file_length = 18000
##    num_steps = 51
#    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
#    
##    sweep_time = 200#6000 #300 #1000 #3000
#    ## channels   
##    pi_ge = pi_ge_time_setting
#    
#    ge_amp = 0.0251*0 #ge_amp_setting
##    ssm_ge = ssm_ge_setting
#    ro_pulse_dur = 5000
#   
#    readout_amp = 0.035 #0.0261 #1#0.5# 1
#    readout_dur = ro_pulse_dur#8000 #13000 #1000
#    
##    phase_offset = mixer_offset
#    
#    coaxtime = sweep_time #500*4
#    pulse_len_add = 2000
#    
##    if q == 0: #qubit 2
##    rabi_ge = Pulse(start=file_length-readout_dur, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
##    ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
###        rabi_ge.phase = 90+phase_offset
###        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
##   
##    if q == 1: #qubit 1
##        rabi_ge = Pulse(start=file_length-readout_dur, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
##        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
###        rabi_ge.phase = 90+phase_offset
###        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#
#    
#
#    if q == 0: #qubit 2        
#        
#         # off-resonance drive for RO (of Q1)
#        coax_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
##        coax_drive.phase = 90
##        the_seq.add_sweep(channel=2, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
#        
##        
#        # off-resonance drive for Q1
#        q2_off_resonance_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
##        q2_off_resonance_drive.phase = 90+phase_offset
##        the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
#
#        # ,gaussian_bool=True, ff=1
#
##        pi_ge = Pulse(start=file_length-readout_dur-10, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
##        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
##        pi_ge.phase = 90+phase_offset
##        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#    
#    
#        sigma = 200#40# 200+0*300 #100
#        width = np.int(sigma*5/2)
#
#        sweeptime = coaxtime
#        if ifsideband == 1:
#            for k in np.arange(num_steps):
#                if k >= 0:  # Note (WC): The first pulse with nonzero duration should larger than 2*width
#
#                # for off-resonance drives to qubit and RO cavity, both of which have the same duration and starting time
##                  pulseLength = np.int_(2*width + sweeptime * k/(num_steps-1))
#                    pulseLength = np.int_(sweeptime * k/(num_steps-1)) # for qubit drive
#                    pulseLength2 = np.int_(pulse_len_add + sweeptime * k/(num_steps-1)) # for RO cavity drive
#                    pulseStart = file_length - readout_dur - pulseLength #file_length-pi_ge_time-pulseLength - 1000
#                    pulseStart2 = file_length - readout_dur - pulseLength2 #file_length-pi_ge_time-pulseLength - 1000
#
#                    tdata = np.arange(pulseLength)
#                    tdata2 = np.arange(pulseLength2)
#                    pdata1 = np.zeros(pulseLength)
#                    pdata2 = np.zeros(pulseLength2)
#
#                    # gaussian shape                    
#                    pdata1[0:width] = amp_Q2 * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
#                    pdata2[0:width] = amp_coax * np.exp(-(tdata2[0:width]-width)**2/(2*sigma**2))
#                    
#                    # sine shape
##                    pdata1[0:width] = amp_Q2 * np.sin(np.pi/2/width*tdata[0:width]) 
##                    pdata2[0:width] = amp_coax * np.sin(np.pi/2/width*tdata[0:width])
#                    
##                    secondStop = np.int_(width+sweeptime * k/(num_steps-1))
#                    secondStop = np.int_(-width+sweeptime * k/(num_steps-1))
#                    secondStop2 = np.int_(-width+ +pulse_len_add + sweeptime * k/(num_steps-1))
#
#
#                    pdata1[width:secondStop] = amp_Q2
#                    pdata2[width:secondStop2] = amp_coax
#                    
#                    pdata1[secondStop:pulseLength] = (amp_Q2 *
#                          np.exp(-(tdata[secondStop:pulseLength]-secondStop)**2/(2*sigma**2)))
#                    pdata2[secondStop2:pulseLength2] = (amp_coax *
#                          np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
#                          
#                    
##                    the_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
#                    ringupdown_seq.channel_list[3][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
#                    
#                    ringupdown_seq.channel_list[0][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
##                    the_seq.channel_list[1][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
#                
##                    plt.figure()
##                    plt.plot(pdata1)
##                    plt.show()
#
#
#    if q == 0: #qubit 2
#        rabi_ge = Pulse(start=file_length-readout_dur, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
#        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
##        rabi_ge.phase = 90+phase_offset
##        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#   
#    if q == 1: #qubit 1
#        rabi_ge = Pulse(start=file_length-readout_dur, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
#        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
##        rabi_ge.phase = 90+phase_offset
##        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    
#    #HET   
#        #Q1 Readout
#    #if q == 1:
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
#    
#    
#
#    #Q2 Readout
#    #if q == 0:
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
##    main_pulse.phase = 90
##    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#    
##    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
##    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
#   
#    
#    ## markers
#    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
#    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
#    
#    ##create the gate for ch1 an ch2
#
#    ## view output
#    if True:
#        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
#        channel2_ch = ringupdown_seq.channel_list[1][0]
#        channel3_ch = ringupdown_seq.channel_list[2][0]
#        channel4_ch = ringupdown_seq.channel_list[3][0]
#        marker1 = ringupdown_seq.channel_list[0][2]
#        
#        plt.plot(channel4_ch[10,11500:11700])
#        
##        channel = channel1_ch + channel3_ch + marker1
##        plt.figure()
##        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
##        plt.show()
##        
##        plt.figure()
##        plt.imshow(channel[:,6000:8000], aspect='auto')
##        plt.show()
#        
#    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#    ringupdown_seq.load_sequence('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0)
#    
###END geom

def rabi_ge_withRamanDrive(num_steps=51,sweep_time=10000,ssm_ge=-0.15,ROIF1=0.1,ROIF2=0.1,ssm_coax=-0.03, amp_coax=1, ssm_Q2=-0.003, amp_Q2=0.6, q=0,ifsideband=1): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
#    pi_ge = pi_ge_time_setting
    
    ge_amp = 0.0251 #ge_amp_setting
#    ssm_ge = ssm_ge_setting
    ro_pulse_dur = 5000
   
    readout_amp = 0.035 #0.0261 #1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
#    phase_offset = mixer_offset
    
    coaxtime = sweep_time #500*4
    pulse_len_add = 2000
    
#    if q == 0: #qubit 2
#    rabi_ge = Pulse(start=file_length-readout_dur, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
#    ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
##        rabi_ge.phase = 90+phase_offset
##        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#   
#    if q == 1: #qubit 1
#        rabi_ge = Pulse(start=file_length-readout_dur, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
#        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
##        rabi_ge.phase = 90+phase_offset
##        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)

    

    if q == 0: #qubit 2        
        
         # off-resonance drive for RO (of Q1)
        coax_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
#        coax_drive.phase = 90
#        the_seq.add_sweep(channel=2, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
        
#        
        # off-resonance drive for Q1
        q2_off_resonance_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=q2_off_resonance_drive)
#        q2_off_resonance_drive.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)

        # ,gaussian_bool=True, ff=1

#        pi_ge = Pulse(start=file_length-readout_dur-10, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    
        sigma = 200#40# 200+0*300 #100
        width = np.int(sigma*5/2)

        sweeptime = coaxtime
        if ifsideband == 1:
            for k in np.arange(num_steps):
                if k >= 0:  # Note (WC): The first pulse with nonzero duration should larger than 2*width

                # for off-resonance drives to qubit and RO cavity, both of which have the same duration and starting time
#                  pulseLength = np.int_(2*width + sweeptime * k/(num_steps-1))
                    pulseLength = np.int_(sweeptime * k/(num_steps-1)) # for qubit drive
                    pulseLength2 = np.int_(pulse_len_add + sweeptime * k/(num_steps-1)) # for RO cavity drive
                    pulseStart = file_length - readout_dur - pulseLength #file_length-pi_ge_time-pulseLength - 1000
                    pulseStart2 = file_length - readout_dur - pulseLength2 #file_length-pi_ge_time-pulseLength - 1000

                    tdata = np.arange(pulseLength)
                    tdata2 = np.arange(pulseLength2)
                    pdata1 = np.zeros(pulseLength2)
                    pdata2 = np.zeros(pulseLength2)

                    # gaussian shape                    
                    pdata1[0:width] = amp_Q2 * np.exp(-(tdata2[0:width]-width)**2/(2*sigma**2))
                    pdata2[0:width] = amp_coax * np.exp(-(tdata2[0:width]-width)**2/(2*sigma**2))
                    
                    # sine shape
#                    pdata1[0:width] = amp_Q2 * np.sin(np.pi/2/width*tdata[0:width]) 
#                    pdata2[0:width] = amp_coax * np.sin(np.pi/2/width*tdata[0:width])
                    
#                    secondStop = np.int_(width+sweeptime * k/(num_steps-1))
                    secondStop = np.int_(-width+sweeptime * k/(num_steps-1))
                    secondStop2 = np.int_(-width+ +pulse_len_add + sweeptime * k/(num_steps-1))


                    pdata1[width:secondStop2] = amp_Q2
                    pdata2[width:secondStop2] = amp_coax
                    
                    pdata1[secondStop2:pulseLength2] = (amp_Q2 *
                          np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                    pdata2[secondStop2:pulseLength2] = (amp_coax *
                          np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                          
                    
#                    the_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                    ringupdown_seq.channel_list[3][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata1
                    
                    ringupdown_seq.channel_list[0][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
#                    the_seq.channel_list[1][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
                
#                    plt.figure()
#                    plt.plot(pdata1)
#                    plt.show()


    if q == 0: #qubit 2
        
        rabi_ge = Pulse(start=file_length-readout_dur, duration=-22, amplitude=0.251, ssm_freq=ssm_ge, phase=90) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        
        rabi_ge = Pulse(start=file_length-readout_dur, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
   
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    
    #HET   
        #Q1 Readout
    #if q == 1:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
    
    

    #Q2 Readout
    #if q == 0:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
   
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
#        plt.plot(channel4_ch[3,8000:-5000])
#        plt.plot(channel1_ch[3,8000:-5000],'r-')
#        plt.show()
        
#        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,6000:8000], aspect='auto')
#        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0)
    
##END geom

def rabi_tomography(num_steps=51,sweep_time=200,ssm_ge=-0.15,pi_ge=40, dur_pi2_ge=20,ro_phase=0,ROIF1=0,ROIF2 = 0,q=0,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state

    totlength = sweep_time + 8000 #4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
    file_length = 30000
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
#    pi_ge = pi_ge_time_setting
    
    ge_amp = 0.251 #ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 0.035 #0.261#0.5# 1
    readout_dur = 5000 #ro_pulse_dur#8000 #13000 #1000
    mixer_offset = 0
    phase_offset = mixer_offset
    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-dur_pi2_ge-10, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) 
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
   
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-dur_pi2_ge, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    

    pi2_ge = Pulse(start=file_length-readout_dur, duration=-dur_pi2_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0+ro_phase) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi2_ge)
#    rabi_ef.phase = 90+ro_phase+mixer_orth
#    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=pi2_ge)


    #Readout
    #HOMO
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp )
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    #HET   
        #Q1 Readout
    #if q == 1:
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=ROIF1, phase=-file_length*ROIF1*360 ) 
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse) 
    
    

    #Q2 Readout
    #if q == 0:
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=ROIF2, phase=-file_length*ROIF2*360 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp,ssm_freq=-ROIF2, phase=90 )
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
   
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,6000:8000], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#    ringupdown_seq.load_sequence('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0)
    
##END geom


def encirclingLEP1(rabi_time = 2000,pi_ge=32,ssm_ge = 0.3865,ROIF=0.1,off_set=0,num_steps = 51,phase_tomo=0,amp_tomo=0.5,pi_time=0,pi2_time=0,min_amp=.01,detun=0.005,amp_enc=.25,ph=0,ph_ini=0,coef_in=1,tloop=2000): #this is pulsed readout to ring up and ring down cavity dfor e state
#    file_length = 18000

    totlength = rabi_time + 8000 #4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
    file_length = 30000
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    
    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    
    tloop=rabi_time # rabi_time is the total encircling duration, tloop is duration per encirclement. For one encirclement, set tloop = rabi_time
    
#    rabi_time = 2000
    par=rabi_time/tloop
#    pi_time=0H
#    pi_ge=31
#    pi_ef=29 
#    pi_hf=26
##    par=4/
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.205
#    off_set=0
    readout_dur = 5000
    readout_amp = 0.035# 0.261

#   encircling a Type-I LEP in the parameter space of (J, \Delta) of g-e submanifold
#   initial state: +x, -X
#
    p_pi_ge = Pulse(start=file_length- readout_dur-pi2_time, duration= -pi_ge/2, amplitude=0.251, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 360*0 + 180
#    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    
    
    pt_pi_ge_r = Pulse(start=file_length- readout_dur-pi2_time, duration=0, amplitude=amp_enc, ssm_freq=ssm_ge, phase=0,phase_ini=ph_ini, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
#    pt_pi_ge_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=file_length- readout_dur, duration=-pi2_time, amplitude=amp_tomo, ssm_freq=ssm_ge, phase=0+phase_tomo,phase_ini=ph_ini, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)

        
#    p2_pi_ef = Pulse(start=16995, duration=-pi_time, amplitude=amp_pief, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ef)
#    p2_pi_ef.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ef)

    #main readout 
#    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
#    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp, ssm_freq=ROIF, phase=-file_length*ROIF*360) #original readout_amp=1, duration = 1000     -1000
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
        
    ##create the gate for ch1 an ch2
#    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    the_seq.channel_list[0][1] = qubit_gate
#    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel1_ch[0:200,16800:17000], aspect='auto', extent=[16800,17000,200,0])
#        plt.plot(channel1_ch[3,14600:17000],'b--o')   
#        plt.plot(channel2_ch[50,16000:17000],'r--o')   
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
##    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
    
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#    the_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#    the_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
##END geom

def encirclingLEP1_raman_drive(rabi_time = 20000,pi_ge=44,ssm_ge = -0.15,ROIF=0.1,q=0,ssm_coax=-0.03, amp_coax=1, ssm_Q2=-0.03, amp_Q2=0.5, off_set=0,num_steps = 51,phase_tomo=0,amp_tomo=0.25,pi_time=0,pi2_time=10,min_amp=.01,detun=0.005,amp_enc=.2*0,ph=0,ph_ini=0,coef_in=1,tloop=2000,ifsideband=1,ifload=1): #this is pulsed readout to ring up and ring down cavity dfor e state
    
    totlength = rabi_time + 8000 #4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
    file_length = 30000#50000# 30000
    
        
    readout_dur = 5000 #ro_pulse_dur#8000
#    num_steps = 101
    the_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   

# pulse setting    
    ge_amp =  0.15 #0.251 #ge_amp_setting
#    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
#    pi2_ge_time = pi2_ge_time_setting
#    pi_ef_time = pi_ef_time_setting
#    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
#    ssm_ef = ssm_ef_setting
    readout_amp = 0.035 #0.261 

# raman drive setting    
    coaxtime = rabi_time #500*4
    pulse_len_add = 2000

# LEP encircling setting    
    tloop=rabi_time # rabi_time is the total encircling duration, tloop is duration per encirclement. For one encirclement, set tloop = rabi_time
#    rabi_time = 2000
    par=rabi_time/tloop
    
    
#   encircling a Type-I LEP in the parameter space of (J, \Delta) of g-e submanifold
#   initial state: +x, -X
#
    
    
    if q == 1: #qubit 1

        # off-resonance drive for filter (of Q2)
        coax_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=coax_drive)
        coax_drive.phase = 90+phase_offset
        the_seq.add_sweep(channel=2, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=coax_drive)
        
        # off-resonance drive for Q2
        q2_off_resonance_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
        q2_off_resonance_drive.phase = 90+phase_offset
        the_seq.add_sweep(channel=2, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
        
        
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        pi_ge.phase = 90+phase_offset
        the_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 0: #qubit 2        
        
         # off-resonance drive for RO (of Q1)
        coax_drive = Pulse(start=file_length-readout_dur-pi2_time, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=90) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
#        coax_drive.phase = 90
#        the_seq.add_sweep(channel=2, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
        
#        
        # off-resonance drive for Q1
        q2_off_resonance_drive = Pulse(start=file_length-readout_dur-pi2_time, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=90) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=4, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=q2_off_resonance_drive)
#        q2_off_resonance_drive.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)

        # ,gaussian_bool=True, ff=1

#        pi_ge = Pulse(start=file_length-readout_dur-10, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    
        sigma = 200#200#200 #40# 200+0*300 #100
        width = np.int(sigma*5/2)

        sweeptime = coaxtime
        if ifsideband == 1:
            for k in np.arange(num_steps):
                if k <= 2:  # Note (WC): The first pulse with nonzero duration should larger than 2*width
                   sigma = 50#200#200 #40# 200+0*300 #100
                   width = np.int(sigma*5/2)
  
            # for off-resonance drives to qubit and RO cavity, both of which have the same duration and starting time
#                  pulseLength = np.int_(2*width + sweeptime * k/(num_steps-1))
                pulseLength = np.int_(sweeptime * k/(num_steps-1)) # for qubit drive
                pulseLength2 = np.int_(pulse_len_add + sweeptime * k/(num_steps-1)) # for RO cavity drive
                pulseStart = file_length - readout_dur - np.int(pi2_time) - pulseLength #file_length-pi_ge_time-pulseLength - 1000
                pulseStart2 = file_length - readout_dur - np.int(pi2_time) - pulseLength2 #file_length-pi_ge_time-pulseLength - 1000

                tdata = np.arange(pulseLength)
                tdata2 = np.arange(pulseLength2)
                pdata1 = np.zeros(pulseLength2)
                pdata2 = np.zeros(pulseLength2)

                # gaussian shape                    
#                pdata1[0:width] = amp_Q2 * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
#                pdata2[0:width] = amp_coax * np.exp(-(tdata2[0:width]-width)**2/(2*sigma**2))
                
                # sine shape
                pdata1[0:width] = amp_Q2 * np.sin(np.pi/2/width*tdata2[0:width]) 
                pdata2[0:width] = amp_coax * np.sin(np.pi/2/width*tdata2[0:width])
                
#                    secondStop = np.int_(width+sweeptime * k/(num_steps-1))
                secondStop = np.int_(-width+sweeptime * k/(num_steps-1))
                secondStop2 = np.int_(-width+ +pulse_len_add + sweeptime * k/(num_steps-1))


                pdata1[width:secondStop2] = amp_Q2
                pdata2[width:secondStop2] = amp_coax
                
                # gaussian shape
                pdata1[secondStop2:pulseLength2] = (amp_Q2 *
                      np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                pdata2[secondStop2:pulseLength2] = (amp_coax *
                      np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                      
                # sine shape
#                    pdata1[secondStop:pulseLength] = amp_Q2 *
#                          np.cos(-np.pi/2/width*(tdata[secondStop:pulseLength]-secondStop))
#                    pdata2[secondStop2:pulseLength2] = (amp_coax *
#                          np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                
                
                
#                    the_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                the_seq.channel_list[3][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata1
                
                the_seq.channel_list[0][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
#                    the_seq.channel_list[1][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
            
#                    plt.figure()
#                    plt.plot(pdata1)
#                    plt.show()


    p_pi_ge = Pulse(start=file_length- readout_dur-pi2_time, duration= -pi_ge/2, amplitude=0.251, ssm_freq=ssm_ge, phase=-90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 360*0 + 180
#    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    
    
    pt_pi_ge_r = Pulse(start=file_length- readout_dur-pi2_time, duration=0, amplitude=amp_enc, ssm_freq=ssm_ge, phase=0,phase_ini=ph_ini, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
#    pt_pi_ge_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=file_length- readout_dur, duration=-pi2_time, amplitude=amp_tomo, ssm_freq=ssm_ge, phase=0+phase_tomo,phase_ini=ph_ini, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    
#    p2_pi_ef = Pulse(start=16995, duration=-pi_time, amplitude=amp_pief, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ef)
#    p2_pi_ef.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ef)



#    print(np.shape(the_seq.channel_list))    
#    Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp, ssm_freq=ROIF, phase=-file_length*ROIF*360) #original readout_amp=1, duration = 1000     -1000
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]
        
#        plt.plot(channel1_ch[1,file_length-readout_dur-coaxtime-3000:file_length-readout_dur-coaxtime+4500+4000])
#        plt.plot(channel4_ch[5,20000:-5000],'b-')
#        plt.plot(channel1_ch[5,20000:-5000],'r-')
#        plt.plot(channel4_ch[1,file_length-readout_dur-coaxtime-3000:file_length-readout_dur-coaxtime+4500+4000])
#        plt.show()
#        plt.plot(the_seq.channel_list[0][0][4,pulseStart-50:pulseStart+150])
#        plt.plot(channel3_ch[3,:])

        channel = channel4_ch #channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#        the_seq.load_sequence('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def encirclingLEP1_tunable_raman_drive(rabi_time = 20000,pi_ge=44,ssm_ge = -0.15,ROIF=0.1,q=0,ssm_coax=-0.03, amp_coax=1, ssm_Q2=-0.03, amp_Q2=0.5, off_set=0,num_steps = 51,phase_tomo=0,amp_tomo=0.25,pi_time=0,pi2_time=10,min_amp=.01,detun=0.005,amp_enc=.2*0,ph=0,ph_ini=0,coef_in=1,tloop=2000,ifsideband=1,ifload=1): #this is pulsed readout to ring up and ring down cavity dfor e state
    
    totlength = rabi_time + 8000 #4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
    file_length = 50000#50000# 30000
    
        
    readout_dur = 5000 #ro_pulse_dur#8000
#    num_steps = 101
    the_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   

# pulse setting    
    ge_amp =  0.15 #0.251 #ge_amp_setting
#    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
#    pi2_ge_time = pi2_ge_time_setting
#    pi_ef_time = pi_ef_time_setting
#    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
#    ssm_ef = ssm_ef_setting
    readout_amp = 0.035 #0.261 

# raman drive setting    
    coaxtime = rabi_time #500*4
    pulse_len_add = 2000

# LEP encircling setting    
    tloop=rabi_time # rabi_time is the total encircling duration, tloop is duration per encirclement. For one encirclement, set tloop = rabi_time
#    rabi_time = 2000
    par=rabi_time/tloop
    
    
#   encircling a Type-I LEP in the parameter space of (J, \Delta) of g-e submanifold
#   initial state: +x, -X
#
    
    
    if q == 1: #qubit 1

        # off-resonance drive for filter (of Q2)
        coax_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=coax_drive)
        coax_drive.phase = 90+phase_offset
        the_seq.add_sweep(channel=2, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=coax_drive)
        
        # off-resonance drive for Q2
        q2_off_resonance_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
        q2_off_resonance_drive.phase = 90+phase_offset
        the_seq.add_sweep(channel=2, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
        
        
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        pi_ge.phase = 90+phase_offset
        the_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 0: #qubit 2        
        
#         # off-resonance drive for RO (of Q1)
#        coax_drive = Pulse(start=file_length-readout_dur-pi2_time, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=90) #pulse is also a class p is an instance
#        the_seq.add_sweep(channel=1, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
##        coax_drive.phase = 90
##        the_seq.add_sweep(channel=2, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
        
#        
        # off-resonance drive for Q1
        q2_off_resonance_drive = Pulse(start=file_length-readout_dur-pi2_time, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=90) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=4, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=q2_off_resonance_drive)
#        q2_off_resonance_drive.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)

        # ,gaussian_bool=True, ff=1

#        pi_ge = Pulse(start=file_length-readout_dur-10, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
    
    
        sigma = 200#200#200 #40# 200+0*300 #100
        width = np.int(sigma*5/2)

        sweeptime = coaxtime
        if ifsideband == 1:
            for k in np.arange(num_steps):
                if k <= 2:  # Note (WC): The first pulse with nonzero duration should larger than 2*width
                   sigma = 50#200#200 #40# 200+0*300 #100
                   width = np.int(sigma*5/2)
  
            # for off-resonance drives to qubit and RO cavity, both of which have the same duration and starting time
#                  pulseLength = np.int_(2*width + sweeptime * k/(num_steps-1))
                pulseLength = np.int_(sweeptime * k/(num_steps-1)) # for qubit drive
                pulseLength2 = np.int_(pulse_len_add + sweeptime * k/(num_steps-1)) # for RO cavity drive
                pulseStart = file_length - readout_dur - np.int(pi2_time) - pulseLength #file_length-pi_ge_time-pulseLength - 1000
                pulseStart2 = file_length - readout_dur - np.int(pi2_time) - pulseLength2 #file_length-pi_ge_time-pulseLength - 1000

                tdata = np.arange(pulseLength)
                tdata2 = np.arange(pulseLength2)
                pdata1 = np.zeros(pulseLength2)
                pdata2 = np.zeros(pulseLength2)

                # gaussian shape                    
#                pdata1[0:width] = amp_Q2 * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
#                pdata2[0:width] = amp_coax * np.exp(-(tdata2[0:width]-width)**2/(2*sigma**2))
                
                # sine shape
                pdata1[0:width] = amp_Q2 * np.sin(np.pi/2/width*tdata2[0:width]) 
                pdata2[0:width] = amp_coax * np.sin(np.pi/2/width*tdata2[0:width])
                
#                    secondStop = np.int_(width+sweeptime * k/(num_steps-1))
                secondStop = np.int_(-width+sweeptime * k/(num_steps-1))
                secondStop2 = np.int_(-width+ +pulse_len_add + sweeptime * k/(num_steps-1))


                pdata1[width:secondStop2] = amp_Q2
                pdata2[width:secondStop2] = amp_coax
                
                # gaussian shape
                pdata1[secondStop2:pulseLength2] = (amp_Q2 *
                      np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                pdata2[secondStop2:pulseLength2] = (amp_coax *
                      np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                      
                # sine shape
#                    pdata1[secondStop:pulseLength] = amp_Q2 *
#                          np.cos(-np.pi/2/width*(tdata[secondStop:pulseLength]-secondStop))
#                    pdata2[secondStop2:pulseLength2] = (amp_coax *
#                          np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                
                
                
#                    the_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                the_seq.channel_list[3][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata1
                
                the_seq.channel_list[0][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
#                    the_seq.channel_list[1][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
            
#                    plt.figure()
#                    plt.plot(pdata1)
#                    plt.show()

#        coax_drive = Pulse(start=file_length-readout_dur-pi2_time, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=90) #pulse is also a class p is an instance
#        the_seq.add_sweep(channel=1, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)


# detuned cavity drive (tunable amplitude-- use encircling pulse )
    cavity_drive = Pulse(start=file_length- readout_dur-pi2_time, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0,phase_ini=180, t_loop=tloop, ff=1,detun_NH=0,jmin=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=cavity_drive)
#    pt_pi_ge_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)

    p_pi_ge = Pulse(start=file_length- readout_dur-pi2_time, duration= -pi_ge/2, amplitude=0.251, ssm_freq=ssm_ge, phase=-90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 360*0 + 180
#    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    
    
    pt_pi_ge_r = Pulse(start=file_length- readout_dur-pi2_time, duration=0, amplitude=amp_enc, ssm_freq=ssm_ge, phase=0,phase_ini=ph_ini, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
#    pt_pi_ge_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=file_length- readout_dur, duration=-pi2_time, amplitude=amp_tomo, ssm_freq=ssm_ge, phase=0+phase_tomo,phase_ini=ph_ini, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    
#    p2_pi_ef = Pulse(start=16995, duration=-pi_time, amplitude=amp_pief, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ef)
#    p2_pi_ef.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ef)




#    print(np.shape(the_seq.channel_list))    
#    Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp, ssm_freq=ROIF, phase=-file_length*ROIF*360) #original readout_amp=1, duration = 1000     -1000
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]
        
#        plt.plot(channel1_ch[1,file_length-readout_dur-coaxtime-3000:file_length-readout_dur-coaxtime+4500+4000])
#        plt.plot(channel4_ch[5,20000:-5000],'b-')
        plt.plot(channel1_ch[50,78000:-5000],'r-')
#        plt.plot(channel4_ch[1,file_length-readout_dur-coaxtime-3000:file_length-readout_dur-coaxtime+4500+4000])
        plt.show()
#        plt.plot(the_seq.channel_list[0][0][4,pulseStart-50:pulseStart+150])
#        plt.plot(channel3_ch[3,:])

        channel = channel4_ch #channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#        the_seq.load_sequence('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def T1_ge_tunable_dissipation_gaussian_pulse_vary_amp(num_steps= 51, sweeptime = 50000,ssm_ge=-.15,pi_ge_time=100,ssm_coax=-0.25, amp_coax=1, ssm_Q2=0.1, amp_Q2=0.2, ROIF=0.1, q=0, ifsideband=1, ifload = 0): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 8000 #4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
    
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
        
    readout_dur = 5000 #ro_pulse_dur#8000
#    num_steps = 101
    the_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp =  0.251 #ge_amp_setting
    
#    ef_amp = 1
#    pi_ge_time = pi_ge_time_setting
#    pi2_ge_time = pi2_ge_time_setting
#    pi_ef_time = pi_ef_time_setting
#    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
#    ssm_ef = ssm_ef_setting
    readout_amp = 0.035 #0.261 
#    oscNum = 6
#    phase_offset = mixer_offset
    
    coaxtime = sweeptime
    pulse_len_add = 2000
    

# LEP encircling setting    
    tloop=np.int(sweeptime/4) # rabi_time is the total encircling duration, tloop is duration per encirclement. For one encirclement, set tloop = rabi_time


    
    
    if q == 1: #qubit 1

        # off-resonance drive for filter (of Q2)
        coax_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=coax_drive)
        coax_drive.phase = 90+phase_offset
        the_seq.add_sweep(channel=2, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=coax_drive)
        
        # off-resonance drive for Q2
        q2_off_resonance_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
        q2_off_resonance_drive.phase = 90+phase_offset
        the_seq.add_sweep(channel=2, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)
        
        
        pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        pi_ge.phase = 90+phase_offset
        the_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
    if q == 0: #qubit 2

#        # off-resonance drive for filter (of Q2)
#        coax_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0,gaussian_bool=True, ff=1) #pulse is also a class p is an instance
#        the_seq.add_sweep(channel=3, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=coax_drive)
#        coax_drive.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='width', start=-500, stop=-coaxtime,initial_pulse=coax_drive)
        
        
#         # off-resonance drive for RO (of Q1)
#        coax_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0) #pulse is also a class p is an instance
#        the_seq.add_sweep(channel=1, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
##        coax_drive.phase = 90
##        the_seq.add_sweep(channel=2, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=coax_drive)
        
#        
        # off-resonance drive for Q1
        q2_off_resonance_drive = Pulse(start=file_length-readout_dur, duration=0, amplitude=amp_Q2, ssm_freq=ssm_Q2, phase=0) #pulse is also a class p is an instance
        the_seq.add_sweep(channel=4, sweep_name='width', start=-pulse_len_add, stop=-coaxtime-pulse_len_add,initial_pulse=q2_off_resonance_drive)
#        q2_off_resonance_drive.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-coaxtime,initial_pulse=q2_off_resonance_drive)

        # ,gaussian_bool=True, ff=1


        sigma = 200#200#200 #40# 200+0*300 #100
        width = np.int(sigma*5/2)

#        sweeptime = coaxtime
        if ifsideband == 1:
            for k in np.arange(num_steps):
                if k <= 2:  # Note (WC): The first pulse with nonzero duration should larger than 2*width
                   sigma = 50#200#200 #40# 200+0*300 #100
                   width = np.int(sigma*5/2)
  
            # for off-resonance drives to qubit and RO cavity, both of which have the same duration and starting time
#                  pulseLength = np.int_(2*width + sweeptime * k/(num_steps-1))
                pulseLength = np.int_(sweeptime * k/(num_steps-1)) # for qubit drive
                pulseLength2 = np.int_(pulse_len_add + sweeptime * k/(num_steps-1)) # for RO cavity drive
                pulseStart = file_length - readout_dur - pulseLength #file_length-pi_ge_time-pulseLength - 1000
                pulseStart2 = file_length - readout_dur - pulseLength2 #file_length-pi_ge_time-pulseLength - 1000

                tdata = np.arange(pulseLength)
                tdata2 = np.arange(pulseLength2)
                pdata1 = np.zeros(pulseLength2)
                pdata2 = np.zeros(pulseLength2)

                # gaussian shape                    
#                pdata1[0:width] = amp_Q2 * np.exp(-(tdata[0:width]-width)**2/(2*sigma**2))
#                pdata2[0:width] = amp_coax * np.exp(-(tdata2[0:width]-width)**2/(2*sigma**2))
                
                # sine shape
                pdata1[0:width] = amp_Q2 * np.sin(np.pi/2/width*tdata2[0:width]) 
                pdata2[0:width] = amp_coax * np.sin(np.pi/2/width*tdata2[0:width])
                
#                    secondStop = np.int_(width+sweeptime * k/(num_steps-1))
                secondStop = np.int_(-width+sweeptime * k/(num_steps-1))
                secondStop2 = np.int_(-width+ +pulse_len_add + sweeptime * k/(num_steps-1))


                pdata1[width:secondStop2] = amp_Q2
                pdata2[width:secondStop2] = amp_coax
                
                # gaussian shape
                pdata1[secondStop2:pulseLength2] = (amp_Q2 *
                      np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                pdata2[secondStop2:pulseLength2] = (amp_coax *
                      np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                      
                # sine shape
#                    pdata1[secondStop:pulseLength] = amp_Q2 *
#                          np.cos(-np.pi/2/width*(tdata[secondStop:pulseLength]-secondStop))
#                    pdata2[secondStop2:pulseLength2] = (amp_coax *
#                          np.exp(-(tdata2[secondStop2:pulseLength2]-secondStop2)**2/(2*sigma**2)))
                
                
                
#                    the_seq.channel_list[2][0][k,pulseStart:pulseStart+pulseLength] *= pdata1
                the_seq.channel_list[3][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata1
                
                the_seq.channel_list[0][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
#                    the_seq.channel_list[1][0][k,pulseStart2:pulseStart2+pulseLength2] *= pdata2
            
#                    plt.figure()
#                    plt.plot(pdata1)
#                    plt.show()

# detuned cavity drive (tunable amplitude-- use encircling pulse )
    cavity_drive = Pulse(start=file_length- readout_dur, duration=0, amplitude=amp_coax, ssm_freq=ssm_coax, phase=0,phase_ini=180, t_loop=tloop, ff=1,detun_NH=0,jmin=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= - sweeptime ,initial_pulse=cavity_drive)
#    pt_pi_ge_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)

    pi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
#        pi_ge.phase = 90+phase_offset
#        the_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)

        
#    print(np.shape(the_seq.channel_list))    
#    Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp, ssm_freq=ROIF, phase=-file_length*ROIF*360) #original readout_amp=1, duration = 1000     -1000
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        marker1 = the_seq.channel_list[0][2]
        
#        plt.plot(channel1_ch[1,file_length-readout_dur-coaxtime-3000:file_length-readout_dur-coaxtime+4500+4000])
#        plt.plot(channel1_ch[50,0:-5000],'r-')
#        plt.plot(channel4_ch[50,0:-5000],'b-')
#        plt.plot(channel4_ch[1,file_length-readout_dur-coaxtime-3000:file_length-readout_dur-coaxtime+4500+4000])
#        plt.show()
#        plt.plot(the_seq.channel_list[0][0][4,pulseStart-50:pulseStart+150])
#        plt.plot(channel3_ch[3,:])

        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
        
    
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        the_seq.load_sequence('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0)
#END T1_GE_tunable_dissipation
        
        
def floqEncirclingLEP1(rabi_time = 2000,pi_ge=32,pi_ef=28,ssm_ge = 0.3865,ssm_ef = 0.912,off_set=0,num_steps = 51,phase_tomo=0,amp_tomo=0.5,pi_time=0,pi2_time=0,min_amp=.01,detun=0.005,amp_enc=.25,ph=0,ph_ini=0,coef_in=1,amp_pief=0.5,tloop=2000): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 201
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    
#    tloop=rabi_time # rabi_time is the total encircling duration, tloop is duration per encirclement. For one encirclement, set tloop = rabi_time
    
#    rabi_time = 2000
    par=rabi_time/tloop
#    pi_time=0
#    pi_ge=31
#    pi_ef=29 
#    pi_hf=26
##    par=4/
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.205
#    off_set=0


#   encircling a Type-I LEP in the parameter space of (J, \Delta) of g-e submanifold
#   initial state: +x, -X
#
    p_pi_ge = Pulse(start=16995-pi_time-pi2_time, duration= -pi_ge, amplitude=0.51-0.01, ssm_freq=ssm_ge, phase=0 + 0*90+0*270)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 360*0 + 180*0 + 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)

#    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-pi_ge-pi_ge/2, duration= -pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 180
#    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
#    
#    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-pi_ge, duration= -pi_ge/2, amplitude=.5*coef_in, ssm_freq=ssm_ge, phase=270+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 360+ph
#    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)

#    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_ge, duration=0, amplitude=amp_enc, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
#    pt_pi_ge_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
#    
#    p2_pi_ge = Pulse(start=16995-pi_time-pi_ef/2, duration=-pi_ge, amplitude=amp_tomo, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)
    
    
    pt_pi_ge_r = Pulse(start=16995-pi2_time-pi_time, duration=0, amplitude=amp_enc, ssm_freq=ssm_ge, phase=0,phase_ini=ph_ini, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi2_time, amplitude=amp_tomo, ssm_freq=ssm_ge, phase=0+phase_tomo,phase_ini=ph_ini, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+phase_tomo
    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)

    
#    pt_pi_ge_r = Pulse(start=16995-pi_time-pi2_time, duration=0, amplitude=amp_enc, ssm_freq=ssm_ge, phase=0,phase_ini=ph_ini, t_loop=tloop, ff=-1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
#    pt_pi_ge_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
#    
#    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi2_time, amplitude=amp_tomo, ssm_freq=ssm_ge, phase=0+phase_tomo,phase_ini=ph_ini, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)

#    p_pi_ge = Pulse(start=16995-pi_time, duration= -pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=90+phase_tomo)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 180+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ge)
    
    p2_pi_ef = Pulse(start=16995, duration=-pi_time, amplitude=amp_pief, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ef)
    p2_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ef)

    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel1_ch[0:200,16800:17000], aspect='auto', extent=[16800,17000,200,0])
#        plt.plot(channel1_ch[50,14600:17000],'b--o')   
#        plt.plot(channel2_ch[50,16000:17000],'r--o')   
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"

    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def encirclingLEP1_tVar(rabi_time = 2000,pi_ge=32,pi_ef=28,ssm_ge = 0.3865,ssm_ef = 0.912,off_set=0,phase_tomo=0,amp_tomo=0.5,pi_time=0,pi2_time=0,min_amp=.01,detun=0.005,amp_enc=.25,ph=0,ph_ini=0,coef_in=1,amp_pief=0.5): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    tloop=rabi_time
#    rabi_time = 2000
    par=rabi_time/tloop
#    pi_time=0
#    pi_ge=31
#    pi_ef=29 
#    pi_hf=26
##    par=4/
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.205
#    off_set=0


#   encircling a Type-I LEP in the parameter space of (J, \Delta) of g-e submanifold
#   initial state: +x, -X
#
    p_pi_ge = Pulse(start=16995-pi_time-pi2_time-rabi_time, duration= -pi_ge/2, amplitude=0.51, ssm_freq=ssm_ge, phase=270*0 + 90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ge)
    p_pi_ge.phase = 360*0 + 180
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi2_time-pi_time, duration=-rabi_time, amplitude=amp_enc, ssm_freq=ssm_ge, phase=0,phase_ini=ph_ini, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=16995-pi_time, duration= -pi2_time, amplitude=amp_tomo, ssm_freq=ssm_ge, phase=0+phase_tomo)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+phase_tomo
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)
        
    p2_pi_ef = Pulse(start=16995, duration=-pi_time, amplitude=amp_pief, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ef)
    p2_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ef)

    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel1_ch[0:200,16800:17000], aspect='auto', extent=[16800,17000,200,0])
#        plt.plot(channel1_ch[50,16600:17000],'b--o')   
#        plt.plot(channel2_ch[50,16000:17000],'r--o')   
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"

    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def encirclingNHHEP(rabi_time = 2000,pi_ge=32,pi_ef=28,ssm_ge = 0.3865,ssm_ef = 0.912,off_set=0,num_steps = 51,phase_tomo=0,amp_tomo=0.5,pi_time=0,pi2_time=0, min_amp=.01,detun=0.005,amp_enc=.25,ph=0,ph_ini=0,coef_in=1,amp_pief=0.5): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 201
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    tloop=rabi_time
#    rabi_time = 2000
    par=rabi_time/tloop
#    pi_time=0
#    pi_ge=32
#    pi_ef=28 
#    pi_hf=26
##    par=4/
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.205
#    off_set=0

    
#
    p_pi_ge = Pulse(start=16995-pi2_time-pi_ef/2-pi_time-pi_ef, duration= -pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)

    p_pi_ge = Pulse(start=16995-pi2_time-pi_time-pi_ef/2, duration= -pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    
    p_pi_ge = Pulse(start=16995-pi2_time-pi_time, duration= -pi_ef/2, amplitude=.5*coef_in, ssm_freq=ssm_ef, phase=270+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 360+ph
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)


#    
    pt_pi_ge_r = Pulse(start=16995-pi2_time-pi_time, duration=0, amplitude=amp_enc, ssm_freq=ssm_ef, phase=0,phase_ini=ph_ini, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi2_time, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=ph_ini, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+phase_tomo
    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)

    
    p2_pi_ge = Pulse(start=16995, duration=-pi_time, amplitude=amp_pief, ssm_freq=ssm_ef, phase=0,phase_ini=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)

    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel1_ch[0:200,16900:17000], aspect='auto', extent=[16800,17000,200,0])
#        plt.plot(channel1_ch[50,16600:17000],'b--o')   
#        plt.plot(channel2_ch[50,16000:17000],'r--o')   
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"

    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom


def encirclingEP(rabi_time = 2000,pi_ge=32,pi_ef=28,ssm_ge = 0.3865,ssm_ef = 0.912,off_set=0,num_steps = 51,phase_tomo=0,amp_tomo=0.5,pi_time=0,min_amp=.01,detun=0.005,amp_enc=.25,ph=0,ph_ini=0,coef_in=1,amp_pief=0.5): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 201
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    tloop=rabi_time
#    rabi_time = 2000
    par=rabi_time/tloop
#    pi_time=0
#    pi_ge=32
#    pi_ef=28 
#    pi_hf=26
##    par=4/
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.205
#    off_set=0

    
#
    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_ef/2-pi_time-pi_ef, duration= -pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)

    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-pi_ef/2, duration= -pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    
    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time, duration= -pi_ef/2, amplitude=.5*coef_in, ssm_freq=ssm_ef, phase=0+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90+ph
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0 , stop= -rabi_time ,initial_pulse=p_pi_ge)


#    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time, duration=0, amplitude=amp_enc, ssm_freq=ssm_ef, phase=90,phase_ini=ph_ini, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time ,initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=ph_ini, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+phase_tomo
    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)

    
    p2_pi_ge = Pulse(start=16995, duration=-pi_time, amplitude=amp_pief, ssm_freq=ssm_ef, phase=0,phase_ini=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)

    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel1_ch[0:200,16900:17000], aspect='auto', extent=[16800,17000,200,0])
#        plt.plot(channel1_ch[50,16600:17000],'b--o')   
#        plt.plot(channel2_ch[50,16000:17000],'r--o')   
        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"

    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def phase_meas_encircling(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=0,detun=0.005,min_amp=0.01,amp_enc=0.25,amp_h=.51,coef_in=.97,exph=6,plt_s=0,in_a=0,in_b=0,orth_phase=0,extime=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
    tloop=rabi_time
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270
    pi_time=0
    


#
    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-rabi_time-extime, duration=-pi_ge, amplitude=.518, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-pi_ef-rabi_time-extime, duration=-pi_ef, amplitude=.514, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-rabi_time-extime, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    

    p2_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2-rabi_time-extime, duration=-pi_ef/2, amplitude=.514*coef_in, ssm_freq=ssm_ef, phase=0+esph+exph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph+exph+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_hf/2-pi_ef-extime, duration=-rabi_time, amplitude=amp_enc, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)
    
#    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase' , stop= -rabi_time, initial_pulse=p2_pi_ge)

    
    p_pi_ef_r = Pulse(start=16995-pi_hf/2-pi_ef-extime, duration=-pi_ef/2, amplitude=.514*coef_in, ssm_freq=ssm_ef, phase=esph+eph+exph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90+esph+eph+exph+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)


    rabi_ef = Pulse(start=16995-pi_ef-extime, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)
    rabi_ef.phase = 0+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=.514, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
    if plt_s is 0:
        pass
    else:
        plt.plot(channel1_ch[plt_s,in_a:in_b],'b--o')
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)



def st_phase_meas_encircling(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=0,detun=0.005,min_amp=0.01,amp_enc=0.25,amp_h=.51,coef_in=.97,exph=6,plt_s=0,in_a=0,in_b=0,orth_phase=0,extime=0,ampf=.514,ampe=.518): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
    tloop=rabi_time
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270
    pi_time=0
    


#
    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-rabi_time-extime, duration=-pi_ge, amplitude=ampe, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-pi_ef-rabi_time-extime, duration=-pi_ef, amplitude=ampf, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-rabi_time-extime, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    

    p2_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2-rabi_time-extime, duration=-pi_ef/2, amplitude=ampf*coef_in, ssm_freq=ssm_ef, phase=0+esph+exph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph+exph+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_hf/2-pi_ef-extime, duration=-rabi_time, amplitude=amp_enc, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)
    
#    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase' , stop= -rabi_time, initial_pulse=p2_pi_ge)

    
    p_pi_ef_r = Pulse(start=16995-pi_hf/2-pi_ef-extime, duration=-pi_ef/2, amplitude=ampf*coef_in, ssm_freq=ssm_ef, phase=esph+eph+exph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90+esph+eph+exph+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)


    rabi_ef = Pulse(start=16995-pi_ef-extime, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)
    rabi_ef.phase = 0+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=ampf, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
    if plt_s is 0:
        pass
    else:
        plt.plot(channel1_ch[plt_s,in_a:in_b],'b--o')
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)







def check_phase_enc(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 1000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=0,exph=0,detun=0.005,min_amp=0.01,amp_enc=0.25,amp_h=.5,coef_in=1): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
    tloop=rabi_time
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270
    pi_time=0


#
    p_pi_ge = Pulse(start=16995-3*pi_ef, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_ef-pi_ef, duration=-pi_ef, amplitude=.5*coef_in, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
#    
#    p_pi_hf = Pulse(start=16995-pi_ef-pi_ef, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_hf)
#    p_pi_hf.phase = 0
#    the_seq.add_sweep(channel=2,  sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p_pi_hf)
#    

    p2_pi_ge = Pulse(start=16995-pi_ef-pi_ef/2, duration=-pi_ef/2, amplitude=.5*coef_in, ssm_freq=ssm_ef, phase=0+esph+exph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph+exph
    the_seq.add_sweep(channel=2,  sweep_name='start',start=0, stop=-rabi_time,initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_ef, duration=0, amplitude=amp_enc, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width',start=0, stop=-rabi_time,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='width',start=0, stop=-rabi_time,initial_pulse=pt_pi_ge_r)
    
#    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase' , stop= -rabi_time, initial_pulse=p2_pi_ge)
    p2_pi_ge = Pulse(start=16995-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+esph+eph,phase_ini=0, t_loop=tloop,ff=None,detunlinear=detun,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph+eph
    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase'  ,start=0, stop= -rabi_time, initial_pulse=p2_pi_ge)
    
#    p_pi_ef_r = Pulse(start=16995-pi_hf/2-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=180+esph+eph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
#    p_pi_ef_r.phase = 270+esph+eph
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)

#
#    rabi_ef = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',start=0, stop=360,initial_pulse=rabi_ef)
#    rabi_ef.phase = 0
#    the_seq.add_sweep(channel=2,  sweep_name='none',start=0, stop=360,initial_pulse=rabi_ef)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.plot(channel1_ch[50,13000:17000],'b--o')
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)

def xmp_projection_phase(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=0,detun=0.005,min_amp=0.01,amp_enc=0.25,amp_h=.51,coef_in=.97,exph=6,plt_s=0,in_a=0,in_b=0,orth_phase=0,extime=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
    tloop=rabi_time
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270
    pi_time=0
    


#
    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-pi_ef-pi_hf-rabi_time-extime, duration=-pi_ge, amplitude=.484, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-pi_ef-pi_ef-pi_hf-rabi_time-extime, duration=-pi_ef, amplitude=.49, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-pi_ef-pi_hf-rabi_time-extime, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    

    p2_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-pi_ef/2-pi_hf-rabi_time-extime, duration=-pi_ef/2, amplitude=.49*coef_in, ssm_freq=ssm_ef, phase=0+esph+exph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph+exph+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_ef-pi_time-pi_hf/2-pi_ef-pi_hf-extime, duration=-rabi_time, amplitude=amp_enc, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)
    
#    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase' , stop= -rabi_time, initial_pulse=p2_pi_ge)

    
    p_pi_ef_r = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-pi_hf-extime, duration=-pi_ef/2, amplitude=.49*coef_in, ssm_freq=ssm_ef, phase=0+esph+eph+exph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90+esph+eph+exph+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)


    p_pi_ef_r = Pulse(start=16995-extime-pi_ef-pi_ef, duration=-pi_hf, amplitude=amp_h, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 0+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)

    p_pi_ef_r = Pulse(start=16995-extime-pi_ef, duration=-pi_ef, amplitude=.49, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=.49, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90+orth_phase
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
    if plt_s is 0:
        pass
    else:
        plt.plot(channel1_ch[plt_s,in_a:in_b],'b--o')
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)













#def phase_meas_encircling_evolution(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,t loop=2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=0,detun=0.005,min_amp=0.01,amp_enc=0.25,amp_h=.5): #this is pulsed readout to ring up and ring down cavity dfor e state
#    file_length = 18000
##    num_steps = 37
#    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
#
#    ## channels   
##    tloop=rabi_time
##    phase_ini=np.pi/2
##    rabi_time = 4000
##    pi_ge=34
##    pi_ef=28 
##    pi_hf=40
##
##    ssm_ge = 0.3885
##    ssm_ef = 0.0917
##    ssm_hf = 0.243
###    off_set=0
##    esph=270
#    pi_time=0
#
#
##
#    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
#
#    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-pi_ef-rabi_time, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
#    p_pi_ef.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
#    
#    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-rabi_time, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
#    p_pi_hf.phase = 0
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
#    
#
#    p2_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+esph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+esph
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
#    
#    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_hf/2-pi_ef, duration=-rabi_time, amplitude=amp_enc, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
#    pt_pi_ge_r.phase = 180
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)
#    
##    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
##    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', stop= -rabi_time, initial_pulse=p2_pi_ge)
##    p2_pi_ge.phase = 90+phase_tomo
##    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase' , stop= -rabi_time, initial_pulse=p2_pi_ge)
#
#    
#    p_pi_ef_r = Pulse(start=16995-pi_hf/2-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=180+esph+eph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
#    p_pi_ef_r.phase = 270+esph+eph
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
#
#
#    rabi_ef = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)
#    rabi_ef.phase = 0
#    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)
#
#    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
#    p_pi_ef_r.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
#    
#    #main readout 
#    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
#    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
#    
#    
#    ## markers
#    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
#    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
#    
#    ##create the gate for ch1 an ch2
#    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
#    
#    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    the_seq.channel_list[0][1] = qubit_gate
#    the_seq.channel_list[1][1] = qubit_gate
#    ## view output
#    if True:
#        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
#        channel2_ch = the_seq.channel_list[1][0]
#        channel3_ch = the_seq.channel_list[2][0]
#        channel4_ch = the_seq.channel_list[3][0]
##        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
##        plt.show()
##        
#    ## write output
##    write_dir = r"C:\Data\2019\encircling\python_loading"
#
#    write_dir = r"C:\Data\2019\encircling\phase_measurement"
##    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
def phase_meas_encircling_jmin(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,
                               ssm_hf = .2424,pi_hf=20,eph=0,esph=0,detun=0.005,amp_enc=0.25,amp_h=.5,min_amp=0,max_amp=1): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
    tloop=rabi_time
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270
    pi_time=0


#
    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-pi_ef-rabi_time, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-rabi_time, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    

    p2_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+esph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_hf/2-pi_ef, duration=-rabi_time, amplitude=amp_enc, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='Jmin_var',start=min_amp, stop=max_amp,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='Jmin_var',start=min_amp, stop=max_amp,initial_pulse=pt_pi_ge_r)
    
#    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase' , stop= -rabi_time, initial_pulse=p2_pi_ge)

    
    p_pi_ef_r = Pulse(start=16995-pi_hf/2-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=180+esph+eph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 270+esph+eph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)


    rabi_ef = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=amp_h, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabi_ef)
    rabi_ef.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=rabi_ef)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)    
def general_phase_meas(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=0,detun=0.005,min_amp=0.01): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
    tloop=rabi_time
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270
    pi_time=0


#
    p_pi_ge = Pulse(start=16995-pi_hf-3*pi_ef-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-pi_ef-rabi_time, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef-rabi_time, duration=-pi_hf/2, amplitude=1, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    

    p2_pi_ge = Pulse(start=16995-pi_hf/2-pi_ef-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+esph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)
    
#    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_hf/2-pi_ef, duration=-rabi_time, amplitude=.25, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pt_pi_ge_r)
#    pt_pi_ge_r.phase = 180
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_hf/2-pi_ef, duration=-rabi_time, amplitude=.25, ssm_freq=ssm_ef, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=pt_pi_ge_r)   
#    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase' , stop= -rabi_time, initial_pulse=p2_pi_ge)

    
    p_pi_ef_r = Pulse(start=16995-pi_hf/2-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=180+esph+eph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 270+esph+eph
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)


    rabi_ef = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=1, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)
    rabi_ef.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
def test_phase_meas_encircling(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=20,eph=0,esph=0,detun=0.005,min_amp=0.01): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
    tloop=rabi_time
#    phase_ini=np.pi/2
#    rabi_time = 4000
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
##    off_set=0
#    esph=270
    pi_time=0


#
    p_pi_ge = Pulse(start=16995-2*pi_ef-pi_ef, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1,  sweep_name='start',start=0 , stop= -rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,   sweep_name='start',start=0 , stop= -rabi_time,initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_ef-pi_ef, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1,  sweep_name='start',start=0 , stop= -rabi_time,initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,   sweep_name='start',start=0 , stop= -rabi_time,initial_pulse=p_pi_ef)
    
#    p_pi_hf = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=1, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1,  sweep_name='start',start=0 , stop= -rabi_time,initial_pulse=p_pi_hf)
#    p_pi_hf.phase = 0
#    the_seq.add_sweep(channel=2,   sweep_name='start',start=0 , stop= -rabi_time,initial_pulse=p_pi_hf)
#    

    p2_pi_ge = Pulse(start=16995-pi_ef/2-pi_ef, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+esph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start',start=0 , stop= -rabi_time,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+esph
    the_seq.add_sweep(channel=2,   sweep_name='start',start=0 , stop= -rabi_time,initial_pulse=p2_pi_ge)
    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time-pi_ef, duration=0, amplitude=.25, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width',start=0 , stop= -rabi_time,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='width',start=0 , stop= -rabi_time,initial_pulse=pt_pi_ge_r)
    
#    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', stop= -rabi_time, initial_pulse=p2_pi_ge)
#    p2_pi_ge.phase = 90+phase_tomo
#    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase' , stop= -rabi_time, initial_pulse=p2_pi_ge)

    
    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef/2-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=180+esph+eph,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1))#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='detuned_phase', start=0, stop= -rabi_time,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 270+esph+eph
    the_seq.add_sweep(channel=2,  sweep_name='detuned_phase', start=0, stop= -rabi_time,initial_pulse=p_pi_ef_r)

#    rabi_ef = Pulse(start=16995, duration=-pi_hf/2, amplitude=1, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)
#    rabi_ef.phase = 0
#    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=rabi_ef)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6500:7000], aspect='auto', extent=[6500,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
    
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def encirclingEPJvarySweep(rabi_time = 1000,pi_ge=32,pi_ef=28,ssm_ge = 0.3865,ssm_ef = 0.0912,off_set=0,num_steps = 201,phase_tomo=0,amp_tomo=0,pi_time=0,min_amp=.01,max_amp=0.2,detun=-0.005,amp_enc=.25,ph=0,pl=1,amp_tomo_pi2=.5): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 201
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    tloop=rabi_time
#    rabi_time = 2000
    par=rabi_time/tloop
#    pi_time=0
#    pi_ge=32
#    pi_ef=28 
#    pi_hf=26
##    par=4/3
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.205
#    off_set=0


#
    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_ef/2-pi_time-pi_ef-rabi_time, duration= -pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-pi_ef/2-rabi_time, duration= -pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none' ,initial_pulse=p_pi_ge)
    
    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-rabi_time, duration= -pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90+ph
    the_seq.add_sweep(channel=2,  sweep_name='none' ,initial_pulse=p_pi_ge)


#    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time, duration=-rabi_time, amplitude=amp_enc, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='Jmin_var' ,start=min_amp , stop= max_amp ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='Jmin_var' ,start=min_amp , stop= max_amp ,initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo_pi2, ssm_freq=ssm_ef, phase=0+phase_tomo,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+phase_tomo
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)
    
    p2_pi_ge = Pulse(start=16995, duration=-pi_time, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0,phase_ini=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)

    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel1_ch[0:200,14900:17000], aspect='auto', extent=[16800,17000,200,0])
#        plt.plot(channel3_ch[50,14950:15050],'b--o')   
#        plt.plot(channel2_ch[num_steps-pl,14800:17000],'r--o')   
#        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"

    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)    
def encirclingEPJvary(rabi_time = 1000,pi_ge=32,pi_ef=28,ssm_ge = 0.3865,ssm_ef = 0.0912,off_set=0,num_steps = 201,phase_tomo=0,amp_tomo=0,pi_time=0,min_amp=.01,
                      detun=-0.005,amp_en=.25,ph=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 201
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    tloop=rabi_time
#    rabi_time = 2000
    par=rabi_time/tloop
#    pi_time=0
#    pi_ge=32
#    pi_ef=28 
#    pi_hf=26
##    par=4/3
#    ssm_ge = 0.3865
#    ssm_ef = 0.0912
#    ssm_hf = 0.205
#    off_set=0


#
    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_ef/2-pi_time-pi_ef-rabi_time, duration= -pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-pi_ef/2-rabi_time, duration= -pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none' ,initial_pulse=p_pi_ge)
    
    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_time-rabi_time, duration= -pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0+ph)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90+ph
    the_seq.add_sweep(channel=2,  sweep_name='none' ,initial_pulse=p_pi_ge)


#    
    pt_pi_ge_r = Pulse(start=16995-pi_ef/2-pi_time, duration=-rabi_time, amplitude=amp_en, ssm_freq=ssm_ef, phase=90,phase_ini=0, t_loop=tloop, ff=1,detun_NH=detun,jmin=min_amp)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' ,initial_pulse=pt_pi_ge_r)
    pt_pi_ge_r.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none' ,initial_pulse=pt_pi_ge_r)
    
    p2_pi_ge = Pulse(start=16995-pi_time, duration=-pi_ef/2, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=0+phase_tomo+ph,phase_ini=0, t_loop=tloop,ff=None,detunlinear=0,detun_NH_phase=detun)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+phase_tomo+ph
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)
    
    p2_pi_ge = Pulse(start=16995, duration=-pi_time, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=90,phase_ini=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p2_pi_ge)

    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel1_ch[0:200,14900:17000], aspect='auto', extent=[16800,17000,200,0])
#        plt.plot(channel3_ch[50,14950:15050],'b--o')   
#        plt.plot(channel2_ch[50,14800:17000],'r--o')   
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
#  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def phase_meas_no_transport(rabi_time = 100,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
    pi_ge=16
    pi_ef=28 
    pi_hf=26

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.205
#    off_set=0

#    p_pi_ge = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#

#
    p_pi_ge = Pulse(start=16995-2*pi_ef-rabi_time, duration=-pi_ge, amplitude=1, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ef = Pulse(start=16995-pi_ef/2-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    


    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_ntransport"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def phase_meas_no_transport_hf(rabi_time = 100,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
    pi_ge=34
    pi_ef=28 
    pi_hf=40

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.243
#    off_set=0

#    p_pi_ge = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#

#
    p_pi_ge = Pulse(start=16995-pi_ef-pi_hf-pi_ef-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-rabi_time, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    p_pi_hf = Pulse(start=16995-pi_hf/2-pi_ef-rabi_time, duration=-pi_hf/2, amplitude=.5, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_hf)
    p_pi_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_hf)
    


    p_pi_ef_r = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=.5, ssm_freq=ssm_hf, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)

    p_pi_ef2 = Pulse(start=16995, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef2)
    p_pi_ef2.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef2)    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_ntransport"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom    
    
def phase_meas_xf(rabi_time = 100,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 37
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
#    rabi_time = 4000
    pi_ge=16
    pi_ef=28 
    pi_hf=26

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.205
#    off_set=0

#    p_pi_ge = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#

#
    p_pi_ge = Pulse(start=16995-pi_ef-2*pi_ge-rabi_time, duration=-pi_ge, amplitude=1, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ef = Pulse(start=16995-pi_ef/2-2*pi_ge-rabi_time, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)
    
    p_pi_ge = Pulse(start=16995-pi_ef/2-pi_ge-rabi_time, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ge = Pulse(start=16995-pi_ef/2, duration=-pi_ge, amplitude=.5, ssm_freq=ssm_ge, phase=270)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)
    
    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=.5, ssm_freq=ssm_ef, phase=0)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='phase',start=0, stop=360,initial_pulse=p_pi_ef_r)
    
    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_ntransport"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def es_evolution(amp=.5): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    rabi_time = 400
    pi_ge=34
    pi_ef=28 
    pi_hf=26

    ssm_ge = 0.3885
    ssm_ef = 0.0917
    ssm_hf = 0.205
    

#
    p_pi_ge = Pulse(start=16995-pi_ge/2, duration=-pi_ge/2, amplitude=amp, ssm_freq=ssm_ge, phase=90)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    
    pt_pi_ef_r = Pulse(start=16995-pi_ge/2, duration=0, amplitude=.5, ssm_freq=ssm_ge, phase=180)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)
    pt_pi_ef_r.phase = 270
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=pt_pi_ef_r)


    #main readout 
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,5000:7000], aspect='auto', extent=[5000,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\test"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
def rabi_ge_with_2pi_phase_shift(): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 101
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels   
#    t_loop=2000
#    phase_ini=np.pi/2
    rabi_time = 1000
    ssm_ge = 0.3885
    p_pi_ge_r = Pulse(start=16995, duration=0, amplitude=1, ssm_freq=ssm_ge, phase=0,phase_ini=0, t_loop=rabi_time, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=p_pi_ge_r)
    p_pi_ge_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=p_pi_ge_r)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    p_pi_ge = Pulse(start=16995, duration=0, amplitude=1, ssm_freq=ssm_ge, phase=90,phase_ini=np.pi/2, t_loop=rabi_time, ff=1)#, phase_ini=np.pi/2, t_loop=400, ff=1) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1)
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\test"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def rabi_ef(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,amp_g=.505,amp_tomo=0,p_tomo=0,ph_tomo=0,sig=-1,ef_dur=0,coef_f=0,phas=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
#    pi_ge=34
#    pi_ef=28
#    rabi_time = 2000
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
    p_pi_ge = Pulse(start=16995-p_tomo-ef_dur/2, duration= -pi_ge, amplitude=amp_g, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
 
#    p_pi_ge = Pulse(start=16995-p_tomo-ef_dur/2, duration= -ef_dur, amplitude=.5, ssm_freq=ssm_ef, phase=90) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
#    p_pi_ge.phase = 180
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)

    p_pi_ge = Pulse(start=16995-p_tomo, duration= -ef_dur/2, amplitude=.5*coef_f, ssm_freq=ssm_ef, phase=90+phas+180) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 180+phas+180
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)

#    p_pi_ef = Pulse(start=16995, duration= -pi_ef, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
#    p_pi_ef.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)

    rabi_ef = Pulse(start=16995-p_tomo, duration=0, amplitude=amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_ef)
    rabi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_ef)




    g_ge = Pulse(start=16997, duration=-p_tomo, amplitude=amp_tomo, ssm_freq=ssm_ef, phase=ph_tomo+0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
    g_ge.phase = 90+ph_tomo
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)

    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,16840:17000], aspect='auto', extent=[6840,7000,200,0])
#        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def unbroken_es(phases=0,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
    pi_ge=34
    pi_ef=14
    rabi_time = 50
    ssm_ge = 0.3885
    ssm_ef = 0.0917
    p_pi_ge = Pulse(start=16995-pi_ef-rabi_time, duration= -pi_ge, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-rabi_time, duration= -pi_ef, amplitude=0.5, ssm_freq=ssm_ef, phase=0+phases) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90+phases
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ef)

    rabi_ef = Pulse(start=16995, duration=rabi_time, amplitude=0, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='amplitude', start=0, stop=.3,initial_pulse=rabi_ef)
    rabi_ef.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='amplitude', start=0, stop=.3,initial_pulse=rabi_ef)



#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)

    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6840:7000], aspect='auto', extent=[6840,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\rabi_ef"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def rabi_hf(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,rabi_time = 2000,pi_ge=34,pi_ef=28,ssm_hf = .2424,amp_hf=0,hf_dur=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 101
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
#    pi_ge=34
#    pi_ef=28
#    rabi_time = 100
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = .2424
    p_pi_ge = Pulse(start=16995-2*pi_ef-hf_dur, duration=-pi_ge, amplitude=.495, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_ef-hf_dur, duration=-pi_ef, amplitude=.487, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)

    p_pi_ef = Pulse(start=16995-pi_ef, duration=-hf_dur, amplitude=amp_hf, ssm_freq=ssm_hf, phase=90+90) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
    p_pi_ef.phase = 0+90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=p_pi_ef)
    #p.phase = 90 #make the pulse phase 90 degrees to get the single sideband modulation
    #rabi_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-200,initial_pulse=p)
    rabi_hf = Pulse(start=16995-pi_ef, duration=0, amplitude=amp, ssm_freq=ssm_hf, phase=90) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_hf)
    rabi_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_hf)

    pi_ef_in = Pulse(start=16995, duration=-pi_ef, amplitude=0.487, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=pi_ef_in)
    pi_ef_in.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=pi_ef_in)

    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
#        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def rabi_det(): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 101
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
    off_set=0
    ## channels  
    pi_ge=34
    rabi_time = 100
    ssm_ge = 0.01
    rabi_ge = Pulse(start=16995, duration=-50, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=rabi_ge)
    
    rabi_ge = Pulse(start=16995, duration=0, amplitude=0.5, ssm_freq=ssm_ge+.012, phase=0, detunlinear=0.012) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_ge)


    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.plot(channel2_ch[20,6800:7000])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
    write_dir = r"C:\Data\2019\encircling\rabi"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
#  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    

def rabi(ssm_ge = .3925, num_steps = 51,amp = 0.5, rabi_time = 100,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 201
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
#    off_set=0
    ## channels  
#    pi_ge=34
#    rabi_time = 100
#    ssm_ge = .387

    rabi_ge = Pulse(start=16995, duration=0, amplitude=amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance

#    rabi_ge = Pulse(start=16995, duration=0, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_ge)


    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,16000:17000], aspect='auto', extent=[16000,17000,200,0])
#        plt.plot(channel1_ch[50,16800:17000],'b--o')
#        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True, pr=0)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
#    wx_programs.wx_set_and_amplitude_and_offset()
##END geom

    
def two_level_ep(ssm_ge = .3885, num_steps = 51,amp = 0.5): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 201
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
    off_set=0
    ## channels  
    pi_ge=35
    rabi_time = 2000
#    ssm_ge = 0.3885+0.001
    
    rabi_ge2 = Pulse(start=16995, duration=-pi_ge, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-rabi_time,initial_pulse=rabi_ge2)
    rabi_ge2.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-rabi_time,initial_pulse=rabi_ge2)
    
    rabi_ge = Pulse(start=16995, duration=0, amplitude=amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='width', start=0, stop=-rabi_time,initial_pulse=rabi_ge)

    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
#        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
    write_dir = r"C:\Data\2019\encircling\rabi"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
    
def calibrating_offset(): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
    off_set=2
    ## channels  
    pi_ge=34
    pi_ef=28
    rabi_time=6999
    ssm_ef=0.0935
    ssm_ge = 0.1
    rabi_ge = Pulse(start=6999, duration=-rabi_time, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=rabi_ge)

#    g_ge = Pulse(start=6997, duration=-(pi_ge+2), amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)
#    
#    rabi_ge = Pulse(start=16995-pi_ef, duration=-pi_ge, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabi_ge)
#    rabi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=rabi_ge)
#    
#    rabi_ef = Pulse(start=16995, duration=-pi_ge, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabi_ef)
#    rabi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=rabi_ef)
    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
    write_dir = r"C:\Data\2019\encircling\rabi"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
    
def no_pi_pi_pipi(off_set=2,pi_ge=35,pi_ef=40,ssm_ge = 0.3855,ssm_ef=0.0912,p1=35,p2=40,p3=0,ssm_hf=.222): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
#    off_set=2
    ## channels  
#    pi_ge=34
#    pi_ef=28
#    ssm_ef=0.092
#    ssm_ge = 0.3885
#    rabi_ge = Pulse(start=16995, duration=-pi_ge, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabi_ge)
#    rabi_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=rabi_ge)
###
    
    rabi_ge = Pulse(start=16995-p2-p3, duration=-p1, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabi_ge)
    rabi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=rabi_ge)
    
    rabi_ef = Pulse(start=16995-p3, duration=-p2, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabi_ef)
    rabi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=rabi_ef)
    
#    rabi_ef = Pulse(start=16995, duration=-p3, amplitude=0.5, ssm_freq=ssm_hf, phase=90) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=rabi_ef)
#    rabi_ge.phase = 0
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=rabi_ef)

    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def readout_only(): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 1
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
    off_set=2
    ## channels  
    pi_ge=34
    pi_ef=28
    ssm_ef=.092
    ssm_ge = 0.3885

    #main readout
    main_pulse = Pulse(start = 0,duration = 8000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
    write_dir = r"C:\Data\2019\encircling\rabi"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def t1(ssm_ge = 0.3906,off_set=0,num_steps=51,t1_time = 1000,t_min=0,pi_ge=34,amp_g=.5): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
#    off_set=0
#    pi_ge=34 #36
#    t1_time = 1000
#    ssm_ge = 0.3885 #0.3925 #0.386

    t1_ge = Pulse(start=16995-t_min, duration= -pi_ge, amplitude=amp_g, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t1_ge)
    t1_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t1_ge)

#    g_ge = Pulse(start=6997, duration=50, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)

    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
#    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
#    
#    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    the_seq.channel_list[0][1] = qubit_gate
#    the_seq.channel_list[1][1] = qubit_gate
#    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:100,6000:7000], aspect='auto', extent=[6000,7000,100,0])
#        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"

    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
#    wx_programs.wx_set_and_amplitude_and_offset()
##END geom
def ramsey(ssm_ge = 0.3885,off_set=0,num_steps = 51,t1_time = 1000,pi_ge=34,amp=0.5,pi_echo=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
#    pi_ge=17
#    t1_time = 1000
#    ssm_ge = 0.3885
    t2_ge = Pulse(start=16995-pi_ge/2-pi_echo, duration=-pi_ge/2, amplitude=amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
    t2_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)

#    t2_ge = Pulse(start=16995-pi_ge/2, duration=-pi_echo, amplitude=amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
#    t2_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
    
    t2_ge = Pulse(start=16995, duration=-pi_ge/2, amplitude=amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=t2_ge)
    t2_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=t2_ge)
# some small pulse to open the gate
#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)
    #main readout
    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
#        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
#    wx_programs.wx_set_and_amplitude_and_offset()

##END geom    
def t1_ef(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,t1_time = 6000,t_min=0,pi_ge=34,pi_ef=28,amp_g=.505): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
#    pi_ge=34
#    pi_ef=28   
#    t1_time = 6000
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917

    p_pi_ge = Pulse(start=16995-pi_ef-t_min, duration=-pi_ge, amplitude=amp_g, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ge)

    t1_ef = Pulse(start=16995-t_min, duration=-pi_ef, amplitude=.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t1_ef)
    t1_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t1_ef)

#    g_ge = Pulse(start=6997, duration=100, amplitude=5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=g_ge)
#    g_ge.phase = 90
#    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=g_ge)

    #main readout

    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)

#    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
#    the_seq.add_sweep(channel=2, marker=1, sweep_name='none',initial_pulse=main_pulse)    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate
    the_seq.channel_list[1][1] = qubit_gate
    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:25,6000:7000], aspect='auto', extent=[6000,7000,25,0])
#        plt.plot(channel2_ch[0:0,6000:7000])
#        plt.show()
#        
    ## write output
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
#    write_dir = r"C:\Data\2019\encircling\
## 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
#  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def ramsey_ef_phase(): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
    pi_ge=34
    pi_ef=28    
    t1_time = 1000
    ssm_ge = 0.3885
    ssm_ef = 0.0917
    p_pi_ge = Pulse(start=16995-pi_ef, duration=-pi_ge, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ge)

    t2_ef = Pulse(start=16995-pi_ef/2, duration=-pi_ef/2, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
    t2_ef.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef/2, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=360*5 ,initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=360*5 ,initial_pulse=p_pi_ef_r)


    #main readout

    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
#    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
#    
#    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
#        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\ramsey_ef"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
  
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def ramsey_ef(ssm_ge = 0.3865, ssm_ef = 0.0917,off_set=0,t1_time = 1000,num_steps = 51,pi_ge=34,pi_ef=28,pi_echo=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
#    pi_ge=34
#    pi_ef=28    
#    t1_time = 1000
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
    p_pi_ge = Pulse(start=16995-pi_ef-pi_echo, duration=-pi_ge, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ge)

    t2_ef = Pulse(start=16995-pi_ef/2-pi_echo, duration=-pi_ef/2, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
    t2_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)

    t2_ef = Pulse(start=16995-pi_ef/2, duration=-pi_echo, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
    t2_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
    
    p_pi_ge = Pulse(start=16995, duration=-pi_ef/2, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p_pi_ge)


    #main readout

    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
#    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
#    
#    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
#        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def t1_hf(ssm_ge = 0.3885,ssm_ef = 0.0917,ssm_hf = 0.243,pi_ge=32,pi_ef=28,pi_hf=36,amp=.5,t1_time = 4000,num_steps = 101,off_set=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 101
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=34
#    t1_time = 4000
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf = 0.243
    
    p_pi_ge = Pulse(start=16995-pi_ef-pi_hf-pi_ef, duration=pi_ge, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef, duration=pi_ef, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=p_pi_ef)

    t1_hf = Pulse(start=16995-pi_ef, duration=pi_hf, amplitude=amp, ssm_freq=ssm_hf, phase=90) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t1_hf)
    t1_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t1_hf)

    p_pi_ef_r = Pulse(start=16995, duration=pi_ef, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ef_r)


    #main readout

    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,16000:17000], aspect='auto', extent=[16000,17000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True)
  
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def ramsey_hf(ssm_ge = 0.3865,ssm_ef = 0.0917,off_set=0,num_steps = 51,amp=.5,t2_time= 1000,pi_ge=34,pi_ef=28,ssm_hf = .2424,pi_hf=40,xp=0,p2ef=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
#    num_steps = 51
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
#    pi_ge=34
#    pi_ef=28 
#    pi_hf=40
#    t2_time = 1000
#    ssm_ge = 0.3885
#    ssm_ef = 0.0917
#    ssm_hf=.2424
    
    p_pi_ge = Pulse(start=16995-2*pi_ef-pi_hf-p2ef-xp, duration=-pi_ge, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t2_time,initial_pulse=p_pi_ge)
    p_pi_ge.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t2_time,initial_pulse=p_pi_ge)

    p_pi_ef = Pulse(start=16995-pi_hf-pi_ef-p2ef-xp, duration=-pi_ef, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t2_time,initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t2_time,initial_pulse=p_pi_ef)
    
    
 
    t2_hf = Pulse(start=16995-pi_ef-pi_hf/2-p2ef-xp, duration=-pi_hf/2, amplitude=amp, ssm_freq=ssm_hf, phase=90) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t2_time,initial_pulse=t2_hf)
    t2_hf.phase = 0
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t2_time,initial_pulse=t2_hf)

    p_pi_ef = Pulse(start=16995-pi_hf/2-pi_ef-p2ef/2-xp, duration=-p2ef/2, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t2_time,initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t2_time,initial_pulse=p_pi_ef)

    p_pi_ef = Pulse(start=16995-pi_hf/2-pi_ef-xp, duration=-p2ef/2, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ef)
    p_pi_ef.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ef)

    t2_hf = Pulse(start=16995-pi_ef-pi_hf/2, duration=-xp, amplitude=amp, ssm_freq=ssm_hf, phase=90+90) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t2_time/2,initial_pulse=t2_hf)
    t2_hf.phase = 0+90
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t2_time/2,initial_pulse=t2_hf)
#    
    pi2_hf = Pulse(start=16995-pi_ef, duration=-pi_hf/2, amplitude=amp, ssm_freq=ssm_hf, phase=90+90) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none' , initial_pulse=pi2_hf)
    pi2_hf.phase = 0+90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=pi2_hf)
    

    p_pi_ef_r = Pulse(start=16995, duration=-pi_ef, amplitude=0.5, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=p_pi_ef_r)
    p_pi_ef_r.phase = 90
    the_seq.add_sweep(channel=2,  sweep_name='none', initial_pulse=p_pi_ef_r)


    #main readout

    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
#        plt.plot(channel2_ch[4,16500:17000])
#        plt.imshow(channel2_ch[0:200,15000:17000], aspect='auto', extent=[15000,17000,200,0])
#        plt.show()

#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off_set, write_binary=True,pr=0)
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom
def ch12_leakage(pi_ge=34,pi_ef=28,pi_hf=50,ssm_ge = 0.1734,ssm_ef = 0.1105,ssm_hf = 0.22265,x=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

#    ## channels  
#    pi_ge=34
#    pi_ef=28
#    pi_hf=50
    t2_time = 17000
#    ssm_ge = 0.3885
#    ssm_ef = 0.1105
#    ssm_hf = 0.22265
#    
    p2_pi_ge = Pulse(start=17000, duration=-t2_time, amplitude=1, ssm_freq=ssm_ge, phase=0+x) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90-x
    the_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=p2_pi_ge)


    #main readout

    main_pulse = Pulse(start = 17000,duration =1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel1_ch[0:200,6000:7000], aspect='auto', extent=[1000,17000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\phase_measurement"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
#    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

def ch34_leakage(): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 18000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class

    ## channels  
    pi_ge=34
    pi_ef=28 
    pi_hf=50
    t2_time = 7000
    ssm_ge = 0.3885
    ssm_ef = 0.1105
    ssm_hf = 0.22265
    
    p2_pi_ge = Pulse(start=7000, duration=-t2_time, amplitude=0.5, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t2_time,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90
    the_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t2_time,initial_pulse=p2_pi_ge)


    #main readout

    main_pulse = Pulse(start = 17000,duration = 1000, amplitude= 1 )
    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-2000, duration=1000, amplitude=1)
    the_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2
    the_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
    channel1_channel = the_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    channel2_channel = the_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
    qubit_gate = create_gate(both_ch1_ch2)
    the_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = the_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = the_seq.channel_list[1][0]
        channel3_ch = the_seq.channel_list[2][0]
        channel4_ch = the_seq.channel_list[3][0]
        plt.imshow(channel2_ch[0:200,6000:7000], aspect='auto', extent=[6000,7000,200,0])
        plt.show()
        
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\Data\2019\encircling\ramsey_hf"
# 
    the_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    the_seq.load_sequence('128.252.134.15', base_name='foo', file_path=write_dir, num_offset=0)
##END geom

if __name__ == '__main__':
    pass
#    phase_meas_encircling(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                    ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=180
#                                                    ,detun=-0.005,min_amp=0,amp_enc=.25,amp_h=.5)
#    check_phase_enc(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                    ,off_set=0,num_steps = 51,amp=.5,rabi_time = 1000
#                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=180,esph=180
#                                                    ,detun=-0.005,min_amp=.5,amp_enc=1,amp_h=.5)
#    check_phase_enc()    
#    readout_only()
#    ch12_leakage()
#    calibrating_offset()
#    no_pi_pi_pipi()
#    rabi_ge_with_2pi_phase_shift()
#    rabi()
#    ramsey()
#    t1_ef()
#    rabi_ef()
#    ramsey_ef()
#    t1_hf()
#    rabi_hf()
#    ramsey_hf()
#    rabi_det()
#    es_transport()
#    es_evolution()
#    phase_meas_transport()
#    xtransport()
#    phase_meas_xf()
#    t1()
#    pi2prep()
#    loading()
#    phase_meas_no_transport()
#    ramsey_ef_phase()
#    es_drive()
#    phase_meas_esdrive()
#    phase_meas_partial_transport()
#    phase_meas_partial_transport_hf()
#    xtransport_ef()
#    es_drive_ef()
#    phase_meas_no_transport_hf()
#    phase_meas_transport_hf()
#    phase_meas_esdrive_diff_ang()
#    unbroken_es()
#    two_level_ep()
#    es_drive_jump()
#    loading_jump()
#    wx_programs.wx_set_and_amplitude_and_offset()
#    xtransport_ef_tVar()
#    encirclingEP()
#    phase_meas_encircling()
#    phase_meas_esdrive_hf()
#    xdrive()
#    encirclingEPJvary()
#    phase_meas_encircling_jmin()
#    encirclingLEP()
#    encirclingLEP1()
#    encirclingLEP1_tVar()
#    floqEncirclingLEP1()
#    rabi_ge_withRamanDrive()
#    encirclingLEP1_raman_drive()
#    encirclingLEP1_tunable_raman_drive()
#    T1_ge_tunable_dissipation_gaussian_pulse_vary_amp()