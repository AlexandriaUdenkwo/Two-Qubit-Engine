# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:14:04 2020

@author: crow104
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import datetime
    
import expt_parameters
import seq_experiments
import seq_programs
import daq_programs
import analysis
import math
import black_nonHermitian_ma

from Nop_class import Nop

import wx_programs
import keithley2401


if __name__ == '__main__':
    expt = expt_parameters.expt_parameters()
    seq = Nop("seq")
    msmt = Nop("measurement")
    sweep = Nop("sweep")
    save_dir = r"Z:\Maryam\Ering data\cal/"
    ssb_ge=.387#+.002-.0025#-.212
    ssb_ef=.0914#-.212
    ssb_hf= .243#+.212 #.2418 .4535#
#    ssb_ge=.3903
#    ssb_ef=0.0945
    p_ge=36 #34 #35
    p_ef=28
    p_hf=38

    seq.comment = "measurement"
    seq.num_patterns = 101
#    0.3
    seq.sweep_time = 2000
    seq.num_records_per_pattern = 500
    seq.times = np.linspace(0.0, seq.sweep_time, seq.num_patterns)*1e-3
    seq.num_avgs=5
    
    
#    pi_tomo=0
#    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = rtime,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
#                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,amp_tomo=0.5
#                                                 ,pi_time=pi_tomo,min_amp=-1,max_amp=0,detun=0.005,amp_enc=.25,ph=0,pl=1)
    
#    black_nonHermitian_ma.ch12_leakage(pi_ge=p_ge,pi_ef=p_ef,ssm_ge = ssb_ef,ssm_ef = ssb_ef)    
#    black_nonHermitian_ma.ramsey_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,t1_time = seq.sweep_time,num_steps =seq.num_patterns,pi_ge=p_ge,pi_ef=p_ef)
#    black_nonHermitian_ma.t1_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,t1_time = seq.sweep_time,t_min=0,pi_ge=p_ge,pi_ef=p_ef,amp_g=.505)   
#    black_nonHermitian_ma.ramsey(ssm_ge = ssb_ge,off_set=0,num_steps = seq.num_patterns,t1_time = seq.sweep_time,pi_ge=p_ge,amp=.5)     
#    black_nonHermitian_ma.t1(ssm_ge = ssb_ge,off_set=0,num_steps=seq.num_patterns,t1_time = seq.sweep_time,t_min=0,pi_ge=p_ge,amp_g=.5)
#    black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,
#                                  rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.5,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
#    black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=seq.num_patterns,num_steps = seq.num_patterns,amp=.5,
#                                  rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.5,amp_tomo=0.5,p_tomo=p_ef,ph_tomo=0,sig=1)
#    black_nonHermitian_ma.rabi(ssm_ge = ssb_ge, num_steps = seq.num_patterns,amp = 0.5,rabi_time = seq.sweep_time)  
#    black_nonHermitian_ma.phase_meas_esdrive_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=180)
#    black_nonHermitian_ma.xdrive(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=188,ort=-5)
#    black_nonHermitian_ma.ramsey_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,t2_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf)
#    black_nonHermitian_ma.t1_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,ssm_hf = ssb_hf,pi_ge=p_ge,pi_ef=p_ef,pi_hf=p_hf,amp=.5,t1_time = seq.sweep_time,num_steps = seq.num_patterns,off_set=0)
#    black_nonHermitian_ma.rabi_hf(ssm_ge = ssb_ge,ssm_ef = ssb_ef,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf)
#    black_nonHermitian_ma.phase_meas_encircling(ssm_ge = ssb_ge,ssm_ef = ssb_ef
##                                                ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
###                                                ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=0,esph=0,detun=0.005,jmin=0.01)
    pi_tomo=0
    black_nonHermitian_ma.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=0,
                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=0.003,amp_enc=.2,ph=10)
####    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
####                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=180,off_set=0,amp_tomo=0.5
####                                                 ,pi_time=pi_tomo,min_amp=.01,max_amp=.5,detun=0.003,amp_enc=.17,ph=15)
#####
    pi_tomo=p_ef
    black_nonHermitian_ma.encirclingEP(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=seq.num_patterns,
                                       amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,detun=0.003,amp_enc=.2,ph=10)
###    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
###                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=180,off_set=seq.num_patterns,amp_tomo=0.5,pi_time=pi_tomo,min_amp=.01,max_amp=.5,detun=-0.005)
###    black_nonHermitian_ma.encirclingEPJvarySweep(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
###                                                 ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=180,off_set=seq.num_patterns,amp_tomo=0.5
###                                                 ,pi_time=pi_tomo,min_amp=.01,max_amp=.5,detun=-0.003,amp_enc=.17,ph=15)
###    black_nonHermitian_ma.encirclingEPJvary(rabi_time = seq.sweep_time,pi_ge=p_ge,pi_ef=p_ef,ssm_ge =ssb_ge
###                                       ,ssm_ef = ssb_ef,num_steps = seq.num_patterns,phase_tomo=0,off_set=201,amp_tomo=0.5,pi_time=pi_tomo)
###
###    black_nonHermitian_ma.test_phase_meas_encircling(ssm_ge = ssb_ge,ssm_ef = ssb_ef
###                                                    ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
###                                                    ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=180,esph=0,detun=-0.005,min_amp=.01)   
#    black_nonHermitian_ma.phase_meas_encircling_jmin(ssm_ge = ssb_ge,ssm_ef = ssb_ef
#                                                     ,off_set=0,num_steps = seq.num_patterns,amp=.5,rabi_time = seq.sweep_time
#                                                     ,pi_ge=p_ge,pi_ef=p_ef,ssm_hf = ssb_hf,pi_hf=p_hf,eph=180,esph=0
#                                                     ,detun=-0.005,amp_enc=.25,amp_h=.5,min_amp=0,max_amp=1)
##    black_nonHermitian_ma.no_pi_pi_pipi(off_set=2*seq.num_patterns,pi_ge=p_ge,pi_ef=p_ef,ssm_ge = ssb_ge,ssm_ef= ssb_ef,p1=0,p2=0)
##    black_nonHermitian_ma.no_pi_pi_pipi(off_set=2*seq.num_patterns+1,pi_ge=p_ge,pi_ef=p_ef,ssm_ge = ssb_ge,ssm_ef= ssb_ef,p1=p_ge,p2=0)
##    black_nonHermitian_ma.no_pi_pi_pipi(off_set=2*seq.num_patterns+2,pi_ge=p_ge,pi_ef=p_ef,ssm_ge = ssb_ge,ssm_ef= ssb_ef,p1=p_ge,p2=p_ef)
    black_nonHermitian_ma.loading(num_steps = 2*seq.num_patterns)
    wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5],offset=[-.052,+.053, 0.057, -0.061])
    prev_threshold = [152.385192421315, 157.385192421315];
#############    
##########
###########    daq_params, t_histo, p_readout,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
###########            num_patterns=seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0)
##########
    for k in range(seq.num_avgs):
        
        daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
            num_patterns=2*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=3)
        

        if k is 0:
            p_readout = p
#                   
        else:
            p_readout += p

#    p_post = analysis.p_readout_postselected(p_readout)
    p_post=analysis.p_readout_postselected_pief(p_readout)
    x=seq.times
#    y=p_post[1]
    y=p_post[2]      
    plt.plot(x,y)
#    popt,peer, y_vals, _ = analysis.fit_sine_decay(x,y,guess_vals=None) 
#    popt,peer, y_vals, _ = analysis.fit_exp_decay(x,y,guess_vals=None)
# 
#    popt, y_vals = analysis.fit_sine(x,y,guess_vals=None)
#    ind=np.argmax(y_vals)
#    print(x[ind])
###                 
##        p_readout = p_readout/seq.num_avgs
##        p_post = analysis.p_readout_postselected(p_readout)
##        plt.plot(p_post[1])
##        plt.show()
###


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
#    sweep.vals=np.linspace(.380,.39, 11)
#    freq=[]
#    fp=[]
#    for ink,ss in enumerate(sweep.vals):
#            print(ink)
##            keithley2401.set_current(ss, step_size_mA=0.001)
##            black_nonHermitian_ma.rabi_ef(ssm_ge = ssb_ge,ssm_ef = ss,off_set=0,num_steps = seq.num_patterns,amp=.5,
##                            rabi_time = seq.sweep_time,pi_ge=p_ge,amp_g=.5,amp_tomo=0.5,p_tomo=0,ph_tomo=0,sig=1)
#            black_nonHermitian_ma.ramsey(ssm_ge = ss,off_set=0,num_steps = seq.num_patterns,t1_time = seq.sweep_time,pi_ge=p_ge,amp=.5) 
#            black_nonHermitian_ma.loading(num_steps = 1*seq.num_patterns)
#            wx_programs.wx_set_and_amplitude_and_offset(amp=[1.5, 1.5, 1.5, 1.5])
#            prev_threshold = [152.385192421315, 156.385192421315];     
#            for k in range(seq.num_avgs):
#                
#                daq_params, t_histo, p,a,b = daq_programs.run_daq_auto_threshold_modify_ec(prev_threshold,
#                    num_patterns=1*seq.num_patterns, num_records_per_pattern=seq.num_records_per_pattern,authr=0,fg=2)
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
##            p_post = analysis.p_readout_postselected(p_readout)
#
#            x=seq.times
#            y=p_readout[1]
#            fp.append(y)   
#            plt.plot(x,y)
#            popt,peer, y_vals, _ = analysis.fit_sine_decay(x,y,guess_vals=None) 
#            freq.append(popt[0])

#             


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