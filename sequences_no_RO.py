# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 15:00:27 2023

@author: crow104
"""
import numpy as np
from generator import *
import matplotlib.pyplot as plt

def T1_ge(num_steps= 101, sweeptime = 15000,q): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
        
#    num_steps = 101
    seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    pi_ge = Pulse(start=file_length-q.readout_dur, duration=-q.pi_ge_time, amplitude=q.ge_amp, ssm_freq=q.ssm_ge, phase=0) #pulse is also a class p is an instance
    seq.add_sweep(channel=q.ch, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)

    ## view output
    if True:
        channel1_ch = seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = seq.channel_list[1][0]
        channel3_ch = seq.channel_list[2][0]
        channel4_ch = seq.channel_list[3][0]
        marker1 = seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
        plt.figure()
        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,:], aspect='auto')
#        plt.colorbar()
        plt.show()
    return seq
#END T1_GE

def T1_ef(num_steps= 101, sweeptime = 15000,ssm_ge=-.15,ssm_ef=-0.30,pi_ge_time=20,pi_ef_time=20,q=0, ifload = 1): #this is pulsed readout to ring up and ring down cavity dfor e state
    totlength = sweeptime + 4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
    
#    if file_length > 80000:
#        raise ValueError('File length too long. Make it less than 80000')
#        file_length = 80000
        
    readout_dur = ro_pulse_dur#8000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 100
    ## channels   
    
    ge_amp = ge_amp_setting
    ef_amp = ge_amp_setting
#    pi_ge_time = pi_ge_time_setting
#    pi2_ge_time = pi2_ge_time_setting
#    pi_ef_time = pi_ef_time_setting
#    pi2_ef_time = pi2_ef_time_setting
##    ssm_ge = ssm_ge_setting\
#    ssm_ef = ssm_ef_setting
    readout_amp = 1 
#    oscNum = 6
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    
    if q == 0: #0 for qubit 2
        pi_ge = Pulse(start=file_length-readout_dur-0*pi_ge_time-1*pi_ef_time-50, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        pi_ge.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
        pi_ef = Pulse(start=file_length-readout_dur-0*pi_ge_time-0*pi_ef_time-50, duration=-pi_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ef)
        pi_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ef)
        
    if q == 1: #0 for qubit 1
        pi_ge = Pulse(start=file_length-readout_dur-0*pi_ge_time-1*pi_ef_time-50, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        pi_ge.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ge)
        
        pi_ef = Pulse(start=file_length-readout_dur-0*pi_ge_time-0*pi_ef_time-50, duration=-pi_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ef)
        pi_ef.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweeptime, initial_pulse=pi_ef)
    
#    pi_ef = Pulse(start=file_length-readout_dur-pi_ge_time-50, duration=-pi_ef_time, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ef)
#    pi_ef.phase = 90+phase_offset
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ef)
    
#    pi_ge = Pulse(start=file_length-readout_dur-50, duration=-pi_ge_time, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', initial_pulse=pi_ge)
#    pi_ge.phase = 90+phase_offset
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge)
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1, marker = 2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
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
        plt.imshow(channel[:,:], aspect='auto')
    
#        plt.colorbar()
#        plt.show()
#        plt.plot(channel1_ch[0,31500:32000],'b--o')
#        plt.show()

#        
    if ifload:
        write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
        ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
        ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
        
#end T1_ef
        
        
def ch12_leakage(ssb = -0.205, offset=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 8000
    num_steps = 3
    the_seq = Sequence(file_length, num_steps) #this creates something called the_seq that is an instance of a sequence class
    ## channels  
#    pi_ge=34
#    pi_ef=28
#    pi_hf=50
    pulse_time = 7000
#    ssm_ge = 0.3885
#    ssm_ef = 0.1105
#    ssm_hf = 0.22265
    
    p2_pi_ge = Pulse(start=7000, duration=-pulse_time, amplitude=1, ssm_freq=ssb, phase=0) #pulse is also a class p is an instance
    the_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-pulse_time,initial_pulse=p2_pi_ge)
    p2_pi_ge.phase = 90+offset
    the_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-pulse_time,initial_pulse=p2_pi_ge)


    #main readout

#    main_pulse = Pulse(start = 7000,duration = 1000, amplitude= 1 )
#    the_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)
#    
    
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
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    the_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    the_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])    


def pi_nopi(coef=0,coefpief=0,off=0,ro_dur=8000,ro_amp=1,pi_ge=24,pi_ef=20,ssm_ge = -0.04575,ssm_ef=-0.04575-0.15): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
    num_steps = 3
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
#    pi_ge = pi_ge_time_setting    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = ro_amp#0.5# 1
    readout_dur = ro_pulse_dur#ro_dur #8000#13000 #1000
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
 
    
    pi_2x = Pulse(start=file_length-readout_dur-coefpief*pi_ef, duration=-(pi_ge)/2*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=90) #pulse is also a class p is an instance
    ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=pi_2x)

#    rabi_ef = Pulse(start=file_length-readout_dur, duration=-pi_ef*coefpief, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#    rabi_ef.phase = 90+phase_offset_ef
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
    
#    rabi_ge = Pulse(start=file_length-readout_dur, duration=-pi_ge*coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    rabi_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=4, sweep_name='none', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
    
    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1 ) #original readout_amp=1, duration = 1000     -1000
    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    ##create the gate for ch1 an ch2

    ## view output
    if False:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
        
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length-3000-300:file_length-3000+50], aspect='auto')
#        plt.show()
#        
#        plt.figure()
#        plt.imshow(channel[:,:], aspect='auto')
#        plt.show()
#        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=off, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])

def rabi_ge(num_steps=51,sweep_time=200,ssm_ge=-0.15,ROIF2 = 0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
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
    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp,ssm_freq=ROIF2, phase=0 ) 
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
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom
    
def rabi_ge_2qubit(num_steps=51,sweep_time=200,ssm_ge1=-0.15,ssm_ge2=-0.15, ROIF1 =0.1, ROIF2 = 0.1): #for herterodyne
    file_length = 16000
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    pi_ge = pi_ge_time_setting
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    
#    qubit 1 rabi
    rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge1, phase=0)
    ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
#    
    #qubit 2 rabi
    rabi_ge = Pulse(start=file_length-readout_dur-100, duration=0, amplitude=ge_amp, ssm_freq=ssm_ge2, phase=0)
    ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)

#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp )
#    ringupdown_seq.add_sweep(channel=1,marker = 2, sweep_name='none',initial_pulse=main_pulse)
    #Q1 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF1, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    #main_pulse.phase = 90
    #ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF2, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    # main_pulse.phase = 90
    #ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
 
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
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom    
    
def rabi_ef(num_steps=51,sweep_time=200,pi_ge = 20,ssm_ge = -0.150,ssm_ef=-0.300,ef_amp=1,ROIF= 0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 30000
#    num_steps = 101
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#20000# #300 #1000 #3000
    ## channels   
    
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
#    pi_ge = pi_ge_time_setting
#    ef_amp = ge_amp_setting
#    ssm_ef = ssm_ef_setting
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    buffer = 50
    if q == 0: #qubit 2
        #first pi_ge pulse
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        #drive rabi e-f
        rabi_ef = Pulse(start=file_length-readout_dur-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#        rabi_ef.phase = 90+phase_offset_ef
#        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        
        #second pi_ge-pulse    
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
        
    if q == 1: #qubit 1
        #first pi_ge pulse
        pi_ge_pulse = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        #drive rabi e-f
        rabi_ef = Pulse(start=file_length-readout_dur-pi_ge-buffer, duration=0, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
#        rabi_ef.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ef)
        
        #second pi_ge-pulse    
        pi_ge_pulse = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=pi_ge_pulse)
#        pi_ge_pulse.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
  
#    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp ) #original readout_amp=1, duration = 1000     -1000
#    ringupdown_seq.add_sweep(channel=1,marker = 1, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    

    
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
        plt.imshow(channel[:,:], aspect='auto')
        plt.show()
        
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    
##END geom
    
def spectroscopy_ge(num_steps=101,ssm_start=-.05,ssm_stop=-.15,spec_amp=.5,ROIF =0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    file_length = 16000 #ns
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
#    sweep_time = 200#6000 #300 #1000 #3000
    ## channels   
    
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#2000 #13000 #1000
    
    phase_offset = mixer_offset
    
    if q == 1: #qubit 1
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-7000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
##    
    if q == 0: #qubit 2
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-7000, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#        rabi_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
        
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000     -1000
#    ringupdown_seq.add_sweep(channel=1,marker=2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
#    
    #HET
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
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
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
    return ringupdown_seq
##END geom
    
    
def spectroscopy_ef(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_start=-.15,ssm_stop=-.25,spec_amp=.5,ROIF=0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
    sweep_time = 7000#200#6000 #300 #1000 #3000
    totlength = sweep_time + 4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    ## channels       
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    
    #first pi_ge pulse
    if q == 1: #qubit 1
        pi_ge_pulse = Pulse(start=file_length-readout_dur-sweep_time-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', initial_pulse=pi_ge_pulse)
    #    pi_ge_pulse.phase = 90+phase_offset
    #    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
        
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-sweep_time, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
    #    rabi_ge.phase = 90+phase_offset_ef
    #    ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
#    
    
    if q == 0: #qubit 2
        pi_ge_pulse = Pulse(start=file_length-readout_dur-sweep_time-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
    #    pi_ge_pulse.phase = 90+phase_offset
    #    ringupdown_seq.add_sweep(channel=2, sweep_name='none', initial_pulse=pi_ge_pulse)
        
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=-sweep_time, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
    #    rabi_ge.phase = 90+phase_offset_ef
    #    ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=rabi_ge)
    
    
    
    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
   #HET
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    main_pulse.phase = 90
    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
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
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])


def parametric_coupling(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_start=-.15,ssm_stop=-.25,spec_amp=.5,ROIF1=0,ROIF2=0,q=0):
    
    sweep_time = 1000#200#6000 #300 #1000 #3000
    totlength = sweep_time + 4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
#    num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    ## channels       
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    
    if q==0:
        pi_ge_pulse = Pulse(start=file_length-readout_dur-sweep_time-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=4, sweep_name='none', initial_pulse=pi_ge_pulse)
                
        parametric_12 = Pulse(start=file_length-readout_dur-10, duration=-sweep_time, amplitude=spec_amp, ssm_freq=0, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=2, sweep_name='ssm_freq', start=ssm_start, stop=ssm_stop,initial_pulse=parametric_12)
    
    
     #HET
    #main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
    #ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    #Q1 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF1, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
  
   
    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF2, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)

#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
       ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
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
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
##END geom
    
  
def parametric_coupling_time_domain(num_steps=101,ssm_ge = -0.2,pi_ge =20,ssm_para=0,spec_amp=.5,ROIF1=0,ROIF2=0,q=0,sweep_time=0):
    
    #sweep_time = 1000#200#6000 #300 #1000 #3000
    totlength = sweep_time + 4000
    file_length = 10000 * (np.int(np.ceil(totlength/10000))+1)
    #num_steps = 51
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called rabi_seq that is an instance of a sequence class
    
    ## channels       
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
   
    readout_amp = 1#0.5# 1
    readout_dur = ro_pulse_dur#8000 #13000 #1000
    
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef
    
    if q==0:
        pi_ge_pulse = Pulse(start=file_length-readout_dur-10, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0)
        ringupdown_seq.add_sweep(channel=4, sweep_name='start',start=0, stop=-sweep_time, initial_pulse=pi_ge_pulse)
        
        rabi_ge = Pulse(start=file_length-readout_dur-10, duration=0, amplitude=spec_amp, ssm_freq=ssm_para, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=2, sweep_name='width', start=0, stop=-sweep_time,initial_pulse=rabi_ge)
        
    #HET
    #main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
    #ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1.3*readout_amp,ssm_freq=ROIF1, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
  
   
    #Q2 Readout
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= 1*readout_amp,ssm_freq=ROIF2, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)

#    main_pulse.phase = 90
#    ringupdown_seq.add_sweep(channel=2, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
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
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])
##END geom
    
def ramsey(num_steps = 101,t1_time = 100000,pi_echo_coef=0,osc=0,ssm_ge=-0.15,pi_ge=20,ROIF=0,q=0): #this is pulsed readout to ring up and ring down cavity dfor e state
#    file_length = 10000*(np.int(t1_time/10000)+1) #64000
    totlength = t1_time + 4000
    file_length = 10000 * (np.int(totlength/10000)+2)
    #file_length = 100000#100000# 
#    num_steps = 51
    ge_amp = ge_amp_setting
#    ssm_ge = ssm_ge_setting
#    pi_ge = pi_ge_time_setting
    phase_offset = mixer_offset

#    pi_echo=0
    readout_amp = 1#0.5 # 1
    readout_dur = ro_pulse_dur#8000 #1000
    buffer =0#550#2050
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called ringupdown_seq that is an instance of a sequence class
    ## channels  
#    pi_ge=17
#    t1_time = 1000
#    ssm_ge = 0.3885
    if q == 0: #qubit 2
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2-pi_echo_coef*pi_ge, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
    
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge*pi_echo_coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=4, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
#            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=2,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=4, sweep_name='none',initial_pulse=t2_ge)
#            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=t2_ge)
            
            
            
    if q == 1: #qubit 1
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2-pi_echo_coef*pi_ge, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ge)
    
        t2_ge = Pulse(start=file_length-readout_dur-buffer-pi_ge/2, duration=-pi_ge*pi_echo_coef, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
#        t2_ge.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ge)
        if osc !=0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
#            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=4,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ge)
        elif osc == 0:
            t2_ge = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge/2, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=t2_ge)
#            t2_ge.phase = 90+phase_offset
#            ringupdown_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=t2_ge)
# some small pulse to open the gate
#    g_ge = Pulse(start=6997, duration=100, amplitude=0.5E-20, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=g_ge)
#    g_ge.phase = 90
#    ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=g_ge)
   
#    #trigger pulse to open switch gate
#    gate_trigger = Pulse(start=file_length- readout_dur, duration=readout_dur, amplitude=1)#-500
#    ringupdown_seq.add_sweep(channel=1, marker=2, sweep_name='none', initial_pulse=gate_trigger )
    
    #Readout
#    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude= readout_amp ) #original readout_amp=1, duration = 1000 #-500     -1000
#    ringupdown_seq.add_sweep(channel=1, marker =2, sweep_name='none',initial_pulse=main_pulse)# , marker=2
    main_pulse = Pulse(start = file_length- readout_dur,duration= readout_dur, amplitude=1.3*readout_amp,ssm_freq=ROIF, phase=0 ) 
    ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=main_pulse)
    
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2
#    ringupdown_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = ringupdown_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = ringupdown_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    ringupdown_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length:file_length], aspect='auto')
#        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,file_length-readout_dur-1000-4000:file_length-readout_dur], aspect='auto')
        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])

def ramsey_ef(num_steps = 101,t1_time = 100000,pi_echo_coef=0,osc=0,ssm_ge=-0.15,ssm_ef=-0.3,pi_ge=20,pi_ef = 30,ef_amp=1,q=0):
#    file_length = 10000*(np.int(t1_time/10000)+1) #64000
    totlength = t1_time + 4000
    file_length = 10000 * (np.int(totlength/10000)+2)
    #file_length = 100000#100000# 
#    num_steps = 51
    ge_amp = ge_amp_setting
    phase_offset = mixer_offset
    phase_offset_ef = mixer_offset_ef

    readout_dur = ro_pulse_dur#8000 #1000
    buffer =0#550#2050
    ringupdown_seq = Sequence(file_length, num_steps) #this creates something called ringupdown_seq that is an instance of a sequence class

    if q == 0: #qubit 2
        #first pi_ge pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ef/2-pi_echo_coef*pi_ef-pi_ef/2, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
        #first pi_ef/2 pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ef/2-pi_echo_coef*pi_ef, duration=-pi_ef/2, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
        #echo pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ef/2, duration=-pi_ef*pi_echo_coef, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=1, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=2,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
        
        #second pi_ef/2 pulse
        if osc !=0:
            t2_ef = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ef/2, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=1, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ef)
            t2_ef.phase = 90+phase_offset_ef
            ringupdown_seq.add_sweep(channel=2,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ef)
        elif osc == 0:
            t2_ef = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ef/2, amplitude=ef_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=1, sweep_name='none',initial_pulse=t2_ef)
            t2_ef.phase = 90+phase_offset_ef
            ringupdown_seq.add_sweep(channel=2,  sweep_name='none',initial_pulse=t2_ef)
            
        #second pi_ge pulse
#        t2_ef = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
#        ringupdown_seq.add_sweep(channel=1, sweep_name='none', start=0, stop=-t1_time,initial_pulse=t2_ef)
#        t2_ef.phase = 90+phase_offset
#        ringupdown_seq.add_sweep(channel=2,  sweep_name='none', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
    if q == 1: #qubit 1
        #first pi_ge pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge-pi_ef/2-pi_echo_coef*pi_ef-pi_ef/2, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
        #first pi_ef/2 pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge-pi_ef/2-pi_echo_coef*pi_ef, duration=-pi_ef/2, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
        #echo pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge-pi_ef/2, duration=-pi_ef*pi_echo_coef, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset_ef
        ringupdown_seq.add_sweep(channel=4,  sweep_name='start', start=0, stop=-t1_time/2,initial_pulse=t2_ef)
        
        #second pi_ef/2 pulse
        if osc !=0:
            t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge, duration=-pi_ef/2, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ef)
            t2_ef.phase = 90+phase_offset_ef
            ringupdown_seq.add_sweep(channel=4,  sweep_name='phase', start=0, stop=osc*360,initial_pulse=t2_ef)
        elif osc == 0:
            t2_ef = Pulse(start=file_length-readout_dur-buffer-pi_ge, duration=-pi_ef/2, amplitude=ge_amp, ssm_freq=ssm_ef, phase=0) #pulse is also a class p is an instance
            ringupdown_seq.add_sweep(channel=3, sweep_name='none',initial_pulse=t2_ef)
            t2_ef.phase = 90+phase_offset_ef
            ringupdown_seq.add_sweep(channel=4,  sweep_name='none',initial_pulse=t2_ef)
            
        #second pi_ge pulse
        t2_ef = Pulse(start=file_length-readout_dur-buffer, duration=-pi_ge, amplitude=ge_amp, ssm_freq=ssm_ge, phase=0) #pulse is also a class p is an instance
        ringupdown_seq.add_sweep(channel=3, sweep_name='none', start=0, stop=-t1_time,initial_pulse=t2_ef)
        t2_ef.phase = 90+phase_offset
        ringupdown_seq.add_sweep(channel=4,  sweep_name='none', start=0, stop=-t1_time,initial_pulse=t2_ef)
        
    ## markers
    alazar_trigger = Pulse(start=file_length-readout_dur-1000, duration=1000, amplitude=1)
    ringupdown_seq.add_sweep(channel=3, marker=1, sweep_name='none', initial_pulse=alazar_trigger )
    
    
    ##create the gate for ch1 an ch2
#    ringupdown_seq.add_gate(source_1=1, source_2=2, destination_tuple=(1,1))
    
#    channel1_channel = ringupdown_seq.channel_list[0][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    channel2_channel = ringupdown_seq.channel_list[1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
#    both_ch1_ch2 = channel1_channel**2 + channel2_channel**2
#    qubit_gate = create_gate(both_ch1_ch2)
#    ringupdown_seq.channel_list[0][1] = qubit_gate

    ## view output
    if True:
        channel1_ch = ringupdown_seq.channel_list[0][0] #[channel name -1][0:channel, 1:marker 1, 2:marker 2]
        channel2_ch = ringupdown_seq.channel_list[1][0]
        channel3_ch = ringupdown_seq.channel_list[2][0]
        channel4_ch = ringupdown_seq.channel_list[3][0]
        marker1 = ringupdown_seq.channel_list[0][2]
#        plt.imshow(channel2_ch[0:200,6800:7000], aspect='auto', extent=[6800,7000,200,0])
#        plt.show()
        channel = channel1_ch + channel3_ch + marker1
#        plt.figure()
#        plt.imshow(channel[:,file_length:file_length], aspect='auto')
#        plt.show()
        
        plt.figure()
        plt.imshow(channel[:,file_length-readout_dur-1000-4000:file_length-readout_dur], aspect='auto')
        plt.show()
    ## write output
#    write_dir = r"C:\Data\2019\encircling\python_loading"
#    ringupdown_seq.write_sequence(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0)
    write_dir = r"C:\arbsequences\strong_dispersive_withPython\test_pulse_ringupdown_bin"
    ringupdown_seq.write_sequence_to_disk(base_name='foo', file_path=write_dir, use_range_01=False,num_offset=0, write_binary=True)
    ringupdown_seq.load_sequence_from_disk('128.252.134.31', base_name='foo', file_path=write_dir, num_offset=0, ch_amp=[1,1,1,1])