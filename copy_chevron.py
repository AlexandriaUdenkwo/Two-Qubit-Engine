# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 13:15:36 2021

@author: Crow108
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import FPJPA_sequence as Fs
import wx_programs as wx
import keithley2401
import analysis


instrument_path = r"C:\Users\Crow108\Documents\Python\instr\analyzer"
if instrument_path not in sys.path: sys.path.append(instrument_path )
import bnc

analyzer_path = r"C:\Users\Crow108\Documents\Python\instr\python_interface\python_without_WX2184C"
if analyzer_path not in sys.path: sys.path.append(analyzer_path )
import daq_programs_homo

target_bnc_black_address = 'GPIB0::19::INSTR' # readout
target_bnc_qubit_address = 'USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR'
                            
def sweep_bnc_freq(start_freq, stop_freq, num_points,steps_in_seq_,num_averages_,ro_dur_, IQ_angle):
    freq_list = np.linspace(start_freq, stop_freq, num_points)
    steps_in_seq = steps_in_seq_
    num_averages = num_averages_ # 10000
    ro_dur = ro_dur_ #RO duration
    
    out_I, out_Q = np.zeros( (2, num_points, steps_in_seq))
    for i,freq in enumerate(freq_list):
#        bnc.set_bnc_output(freq, bnc_addr=target_bnc_address_12)
        bnc.set_bnc_output(freq,power_dBm=7, bnc_addr=target_bnc_qubit_address)
        
        daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_pi,counts_pi,rec_readout,rec_avg_all,rec_all = daq_programs_homo.run_daq(steps_in_seq,num_averages,ro_dur,IQangle=IQ_angle)
        plt.plot(Q);plt.show();plt.plot(I);plt.show()

                
        out_Q[i] = Q
        out_I[i] = I
        plt.plot(Q)
        plt.show()
        plt.plot(I)
        plt.show()
        print(i)
        #time.sleep(4)
        
    # make plots
#    save_dir = r"Z:\candle qubit/"
#    np.savetxt(save_dir+'4_18_3pm_' + 'chevron_Q',out_Q)
    plt.imshow(out_Q, extent=[0,steps_in_seq,stop_freq,start_freq],aspect='auto' )
    plt.show()
    plt.imshow(out_I, extent=[0,steps_in_seq,stop_freq,start_freq],aspect='auto' )
    return out_Q, out_I

def sweep_bnc_freq_2q(start_freq, stop_freq, num_points,steps_in_seq_,num_averages_,ro_dur_, ROIF1=0.1,ROIF2=0.1, deg_1 = 0, deg_2 = 0,verbose=True):
    freq_list = np.linspace(start_freq, stop_freq, num_points)
    num_steps = steps_in_seq_
    reps = num_averages_ # 10000
    ro_dur = ro_dur_ #RO duration
    
    out_I1, out_Q1,out_I2, out_Q2 = np.zeros( (4, num_points, num_steps))

    for i,freq in enumerate(freq_list):
#        bnc.set_bnc_output(freq, bnc_addr=target_bnc_address_12)
        bnc.set_bnc_output(freq,power_dBm=7, bnc_addr=target_bnc_qubit_address)
        
        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = 0, deg_2 = 0,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur, verbose=True)
        
        I_Q1 = rec_avg_vs_pats_1[0]
        Q_Q1 = rec_avg_vs_pats_1[1]
        
        
        I_Q2 = rec_avg_vs_pats_2[0]
        Q_Q2 = rec_avg_vs_pats_2[1]

                
        out_Q1[i] = Q_Q1
        out_I1[i] = I_Q1
        plt.plot(Q_Q1)
        plt.show()
        plt.plot(I_Q1)
        plt.show()
        
        out_Q2[i] = Q_Q2
        out_I2[i] = I_Q2
        plt.plot(Q_Q2)
        plt.show()
        plt.plot(I_Q2)
        plt.show()
        
        print(i)
        #time.sleep(4)
        
    # make plots
#    save_dir = r"Z:\candle qubit/"
#    np.savetxt(save_dir+'4_18_3pm_' + 'chevron_Q',out_Q)
    plt.imshow(out_Q1, extent=[0,num_steps,stop_freq,start_freq],aspect='auto' )
    plt.show()
    plt.imshow(out_I1, extent=[0,num_steps,stop_freq,start_freq],aspect='auto' )
    plt.imshow(out_Q2, extent=[0,num_steps,stop_freq,start_freq],aspect='auto' )
    plt.imshow(out_I2, extent=[0,num_steps,stop_freq,start_freq],aspect='auto' )
    plt.show()
    return out_Q1, out_I1,out_Q2, out_I2

def sweep_bnc_power(start_power, stop_power, num_points,steps_in_seq_,num_averages_,ro_dur_):
    power_list = np.linspace(start_power, stop_power, num_points)
#    ssm_ge_list = np.linspace(-0.06, -0.09, num_points) #-0.04905, -0.05905
    steps_in_seq = steps_in_seq_
    num_averages = num_averages_ # 10000
    ro_dur = ro_dur_ #RO duration
    freq = 4.2
    out_I, out_Q = np.zeros( (2, num_points, steps_in_seq))
#    out_y = np.zeros( (num_points, steps_in_seq))
    for i,power in enumerate(power_list):
        target_bnc_address =  'USB0::0x03EB::0xAFFF::141-216340000-0292::INSTR'
        bnc.set_bnc_output(freq,power, bnc_addr=target_bnc_address)
        
        daq_params, rec_readout_vs_pats, p_vs_pats,I,Q,a,b = daq_programs_homo.run_daq(steps_in_seq, num_averages, ro_dur)
                
        out_Q[i] = Q
        out_I[i] = I
        plt.plot(Q)
        plt.show()
        plt.plot(I)
        plt.show()
        
#        ssm_ge = ssm_ge_list[i]
#        Fs.ramsey(steps_in_seq,10000,0,0,ssm_ge,24)
#        wx.wx_set_and_amplitude_and_offset(amp=[1, 1, .89,1])
#        m = [-157,-137]
#        daq_params, t_histo, p,a,b = daq_programs_homo.run_daq_auto_threshold_modify_ec(prev_threshold=m,num_patterns=steps_in_seq, 
#                                                                                       num_records_per_pattern=num_averages,authr=1,fg=2)
#        p_readout = p
#        y=p_readout[1]
#        
#        out_y[i] = y
#        plt.plot(y)
#        plt.show()
        print(i)
        #time.sleep(4)
        
    # make plots
    plt.imshow(out_Q, extent=[0,steps_in_seq,stop_power,start_power],aspect='auto' )
    plt.show()
    plt.imshow(out_I, extent=[0,steps_in_seq,stop_power,start_power],aspect='auto' )
#    plt.imshow(out_y, extent=[0,steps_in_seq,stop_power,start_power],aspect='auto' )
#    return out_y
    return out_Q, out_I
    #return out_Q


def sweep_bnc_freq_and_power(start_power, stop_power, num_points,start_freq, stop_freq, num_pointsfreq,num_averages_,ro_dur_):
    power_list = np.linspace(start_power, stop_power, num_points)
    freq_list = np.linspace(start_freq, stop_freq, num_pointsfreq)
#    ssm_ge_list = np.linspace(-0.06, -0.09, num_points) #-0.04905, -0.05905
    num_averages = num_averages_ # 10000
    ro_dur = ro_dur_ #RO duration
    out_bins_pi, out_counts_pi,out_bins_nopi, out_counts_nopi = [],[],[],[]
#    out_y = np.zeros( (num_points, steps_in_seq))
    
    for i,freq in enumerate(freq_list):
        for j,power in enumerate(power_list):
            target_bnc_address =  'USB0::0x03EB::0xAFFF::141-216340000-0292::INSTR'
            bnc.set_bnc_output(freq,power, bnc_addr=target_bnc_address)
            
            
            Fs.pi_nopi(coef=1,ro_dur=0,ro_amp = 0,pi_ge=27,ssm_ge = -0.201);
            wx.wx_set_and_amplitude_and_offset(amp=[1,1,1,1])
            daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_pi,counts_pi = daq_programs_homo.run_daq(3,num_averages,ro_dur);
            Fs.pi_nopi(coef=0,ro_dur=0, ro_amp = 0,pi_ge=27,ssm_ge = -0.201);
            wx.wx_set_and_amplitude_and_offset(amp=[1,1,1,1])
            daq_params,rec_readout_vs_pats, p_vs_pats,I,Q,bins_nopi,counts_nopi = daq_programs_homo.run_daq(3,num_averages,ro_dur);
            
            plt.plot(bins_pi,counts_pi,label='pi')
            plt.plot(bins_nopi,counts_nopi,label='nopi')
            plt.legend()
            plt.show()
            
            out_bins_pi.append(bins_pi)
            out_counts_pi.append(counts_pi)
            out_bins_nopi.append(bins_nopi) 
            out_counts_nopi.append(counts_nopi)
            
    return out_bins_pi, out_counts_pi,out_bins_nopi, out_counts_nopi
        

##END sweep_bnc_freq
def sweep_bnc_black_freq(start_freq, stop_freq, num_points,steps_in_seq_,num_averages_,ro_dur_):
    freq_list = np.linspace(start_freq, stop_freq, num_points)
    steps_in_seq = steps_in_seq_
    num_averages = num_averages_ # 10000
    ro_dur = ro_dur_ #RO duration
    
    out_I, out_Q = np.zeros( (2, num_points, steps_in_seq))
    for i,freq in enumerate(freq_list):
        bnc.set_bnc_output(freq, bnc_addr=target_bnc_address)
        daq_params, rec_readout_vs_pats, p_vs_pats,I,Q,a,b = daq_programs_homo.run_daq(steps_in_seq, num_averages, ro_dur)
                
        out_Q[i] = Q
        out_I[i] = I
        plt.plot(Q)
        plt.show()
        plt.plot(I)
        plt.show()
        print(i)
        #time.sleep(4)
        
    # make plots
#    save_dir = r"Z:\candle qubit/"
#    np.savetxt(save_dir+'4_18_3pm_' + 'chevron_Q',out_Q)
    plt.imshow(out_Q, extent=[0,steps_in_seq,stop_freq,start_freq],aspect='auto' )
    plt.show()
    plt.imshow(out_I, extent=[0,steps_in_seq,stop_freq,start_freq],aspect='auto' )
    return out_Q, out_I
        
    #return out_Q
##END sweep_bnc_freq    

#def sweep_bnc_freq_cavity_spec(start_freq, stop_freq,ROIF2, num_points, IQ_angle,steps_in_seq_,num_averages_,ro_dur_,wx_amps):
def sweep_bnc_freq_cavity_spec(start_freq,stop_freq,num_points, ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2, steps_in_seq_,num_averages_,ro_dur_,qubit_1_thr,qubit_2_thr,wx_amps):
    
    freq_list = np.linspace(start_freq, stop_freq, num_points)
    steps_in_seq = steps_in_seq_
    num_averages = num_averages_ # 10000
    ro_dur = ro_dur_ #RO duration
    
    out_IQ = np.zeros(num_points)
    out_phase = np.zeros(num_points)
    for i,freq in enumerate(freq_list):
        print(i)
        #        bnc.set_bnc_output(freq, bnc_addr=target_bnc_address_12)
        bnc.set_bnc_output(freq, bnc_addr=target_bnc_black_address)
#        bnc.set_bnc_output(freq, bnc_addr=target_bnc_address_6)

        wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=[0,0.5*2,0,0])
#        daq_params, rec_readout_vs_pats, p_vs_pats,I,Q,a,b,c,d,e = daq_programs.run_daq(steps_in_seq, num_averages, ro_dur,IQangle=0)
        
#        rec_avg_all, rec_readout, rec_avg_vs_pats, rec_all_het, bins, counts = daq_programs_homo.run_daq_het(ssm_if=ROIF2, num_patterns=steps_in_seq, num_records_per_pattern=num_averages,ro_dur=ro_dur, verbose=True)
#        I = rec_avg_vs_pats[0];#plt.plot(I);plt.show()
#        Q = rec_avg_vs_pats[1];#plt.plot(Q);plt.show()


        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=steps_in_seq, num_records_per_pattern=num_averages,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
           
                                                                                                                                                                                                    
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
#        plt.plot(I_Q2);plt.title('I Q2');plt.show()
#        plt.plot(Q_Q2);plt.title('Q Q2');plt.show()   


        out_IQ[i] = Q_Q2[0]**2 + I_Q2[0]**2        
        out_phase[i] = np.arctan(Q_Q2[0]/I_Q2[0])    
         
#        out_IQ[i] = Q[0]**2 + I[0]**2        
#        out_phase[i] = np.arctan(Q[0]/I[0])    
        
    plt.plot(freq_list,out_IQ)
    plt.show()        

    return out_IQ, out_phase

def sweep_bnc_freq_raman(start_freq, stop_freq, num_points,ROIF1,ROIF2,IQ_angle_q1,IQ_angle_q2,steps_in_seq_,num_averages_,ro_dur_,qubit_1_thr,qubit_2_thr,wx_amps):
            
    freq_list = np.linspace(start_freq, stop_freq, num_points)
    steps_in_seq = steps_in_seq_
    num_averages = num_averages_ # 10000
    ro_dur = ro_dur_ #RO duration
    
    out_Q = np.zeros(num_points)
    out_I = np.zeros(num_points)
    pop_e = np.zeros(num_points)
    for i,freq in enumerate(freq_list):
        print(i)
        #        bnc.set_bnc_output(freq, bnc_addr=target_bnc_address_12)
        bnc.set_bnc_output(freq, bnc_addr=target_bnc_qubit_address)
#        bnc.set_bnc_output(freq, bnc_addr=target_bnc_address_6)

        wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=[0,0.5*2,0,0])

        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=steps_in_seq, num_records_per_pattern=num_averages,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
           
                                                                                                                                                                                                    
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
#        plt.plot(I_Q2);plt.title('I Q2');plt.show()
#        plt.plot(Q_Q2);plt.title('Q Q2');plt.show()   
        pop_e[i] = P_Q2[0]
    return pop_e #out_Q, out_I
    
def sweep_rabi_amp(start_amp, stop_amp, num_points,ro_dur,reps,num_steps,IQ_angle,wxamps,ssm_ge,rabi_time):
    sweep_list = np.linspace(start_amp, stop_amp, num_points)
    out_I, out_Q = np.zeros( (2, num_points, num_steps))
    times = np.linspace(0,rabi_time/1000,num_steps)
    for i,val in enumerate(sweep_list):
#        se.rabi_ge(num_steps,rabi_time,ssm_ge,val,0)
        print('current step = ', i)
        print('current sweep val = ', val)
        wxamps[2] = val
        wxamps[3] = val
        print('wx amps=',wxamps)
        wx.wx_set_and_amplitude_and_offset(amp=wxamps)
        
        daq_params, rec_readout_vs_pats, p_vs_pats,I,Q,bins_pi,counts_pi,rec_readout,rec_avg_all,rec_all = daq_programs.run_daq(num_steps,reps,ro_dur,IQangle=IQ_angle)
        plt.plot(times,Q);plt.ylabel('Q');plt.xlabel('Time (\u03BCs)');plt.show();
        plt.plot(times,I);plt.ylabel('I');plt.xlabel('Time (\u03BCs)');plt.show()
        
        out_Q[i] = Q
        out_I[i] = I
        
        print(i)
    
        #time.sleep(4)
        
    # make plots
#    save_dir = r"Z:\candle qubit/"
#    np.savetxt(save_dir+'4_18_3pm_' + 'chevron_Q',out_Q)
    plt.imshow(out_Q, extent=[0,num_steps,stop_amp,start_amp],aspect='auto' )
    plt.show()
    plt.imshow(out_I, extent=[0,num_steps,stop_amp,start_amp],aspect='auto' )
    return out_Q, out_I
##END sweep_rabi_amp
    

def sweep_ssb_freq(start_freq, stop_freq, num_points):
    freq_list = np.linspace(start_freq, stop_freq, num_points)
    steps_in_seq = 101
    num_averages = 100 # 10000
    ro_dur = 8000 #RO duration
    
    out_I, out_Q = np.zeros( (2, num_points, steps_in_seq))
    for i,freq in enumerate(freq_list):
        Fs.e_char(ge_coef=1,ssm_ef = freq) #replace with any sequence but feed in ssb_ge or ssb_ef
        wx_programs.wx_set_and_amplitude_and_offset(amp=[1.7, 1.7, 1.9, 1.9])
        daq_params, rec_readout_vs_pats, p_vs_pats,I,Q,a,b = run_daq(steps_in_seq, num_averages, ro_dur)
                
        out_Q[i] = Q
        out_I[i] = I
        print(i)
        #time.sleep(4)
        
    # make plots
#    save_dir = r"Z:\candle qubit/"
#    np.savetxt(save_dir+'4_18_3pm_' + 'chevron_Q',out_Q)
    plt.imshow(out_Q, extent=[0,steps_in_seq,stop_freq,start_freq],aspect='auto' )
    plt.show()
    plt.imshow(out_I, extent=[0,steps_in_seq,stop_freq,start_freq],aspect='auto' )
    return out_Q, out_I

def sweep_ssb_ef_freq(start_freq, stop_freq, num_points):
    freq_list = np.linspace(start_freq, stop_freq, num_points)
    steps_in_seq = 101
    num_averages = 100 # 10000
    ro_dur = 8000 #RO duration
    
    out_I, out_Q = np.zeros( (2, num_points, steps_in_seq))
    for i,freq in enumerate(freq_list):
        Fs.e_char(ge_coef=1,ssm_ef = freq) #replace with any sequence but feed in ssb_ge or ssb_ef
#        wx_programs.wx_set_and_amplitude_and_offset(amp=[1.7, 1.7, 1.9, 1.9])
        daq_params, rec_readout_vs_pats, p_vs_pats,I,Q,a,b = run_daq(steps_in_seq, num_averages, ro_dur)
        

        num_steps = 51
        sweep_time =100 #ns
        Fs.rabi_ef(num_steps,sweep_time,pi_ge_time,ssm_ge,ssm_ef)
        
        out_Q[i] = Q
        out_I[i] = I
        print(i)
        #time.sleep(4)
        
    # make plots
#    save_dir = r"Z:\candle qubit/"
#    np.savetxt(save_dir+'4_18_3pm_' + 'chevron_Q',out_Q)
    plt.imshow(out_Q, extent=[0,steps_in_seq,stop_freq,start_freq],aspect='auto' )
    plt.show()
    plt.imshow(out_I, extent=[0,steps_in_seq,stop_freq,start_freq],aspect='auto' )
    return out_Q, out_I
##END sweep_ssb_ef_freq

def sweep_rabi_J(start_val, stop_val, num_points,model,scale_matrix,Ax,By,num_steps,reps,ro_dur,sweep_time,pi_ge_time,pi_ef_time,ssb_ge,ssb_ef,IQ_angle):
    value_list = np.linspace(start_val, stop_val, num_points)

    out_y = np.zeros( (num_points, num_steps))
    for i,val in enumerate(value_list):
        print(val)
        print(i)
        ssm_ge = ssb_ge #-0.125
        ssm_ef = val#ssb_ef #-0.407

        ef_amp = 0.1 ##EP occurs at amp = 0.016
    

        Fs.rabi_J(num_steps,sweep_time,ef_amp,pi_ge_time,pi_ef_time,ssm_ge,ssm_ef)
        wx.wx_set_and_amplitude_and_offset(amp = [1, 1, 1,1])
        daq_params, t_histo, p,a,b,n_readout = daq_programs_homo.run_daq_cluster_threshold(model=model,Ax=Ax,By=By,num_patterns=num_steps,
                                                                                 num_records_per_pattern=reps,authr=0,fg=3,ro_dur=ro_dur,IQangle=IQ_angle)
        p_readout = p
#        p_readout= np.matmul(scale_matrix,p) # M*(g;e;f) 3x3 * 3xnum_steps
#
#        plt.plot(p_readout[0],label='|g>');plt.plot(p_readout[1],label='|e>');plt.plot(p_readout[2],label='|f>')
#        plt.legend();plt.title('data scaled with matrix');plt.show()
# 
#        p_post = analysis.p_readout_postselected(p_readout)
#        y = p_post[2]
        y = p_readout[2]
        
        plt.plot(y)
        plt.show()
        out_y[i] = y

        #time.sleep(4)
        
    # make plots
    plt.imshow(out_y, extent=[0,num_steps,stop_val,start_val],aspect='auto' )
    return out_y
##END sweep_rabi_J
    
def sweep_rabi_ef(start_val, stop_val, num_points,model,scale_matrix,Ax,By,num_steps,reps,ro_dur,sweep_time,pi_ge_time,ssb_ge,ssb_ef,IQ_angle,which_qubit):
    value_list = np.linspace(start_val, stop_val, num_points)

    out_y = np.zeros( (num_points, num_steps))
    for i,val in enumerate(value_list):
        print(val)
        print(i)
        ssm_ge = ssb_ge
        ssm_ef = val

        ef_amp = 0.08 ##change for chevron
    

        Fs.rabi_ef(num_steps,sweep_time,pi_ge_time,ssm_ge,ssm_ef,ef_amp,q=which_qubit)
        wx.wx_set_and_amplitude_and_offset(amp = [1, 1, 1,1])
        daq_params, t_histo, p,a,b,n_readout = daq_programs_homo.run_daq_cluster_threshold(model=model,Ax=Ax,By=By,num_patterns=num_steps,
                                                                                 num_records_per_pattern=reps,authr=0,fg=3,ro_dur=ro_dur,IQangle=IQ_angle)
        p_readout = p
        p_readout= np.matmul(scale_matrix,p) # M*(g;e;f) 3x3 * 3xnum_steps

        plt.plot(p_readout[0],label='|g>');plt.plot(p_readout[1],label='|e>');plt.plot(p_readout[2],label='|f>')
        plt.legend();plt.title('data scaled with matrix');plt.show()

        y = p_readout[2]
        
        plt.plot(y)
        plt.show()
        out_y[i] = y
        
    # make plots
    plt.imshow(out_y, extent=[0,num_steps,stop_val,start_val],aspect='auto' )
    return out_y
##END sweep_rabi_ef

def sweep_keithley(start_current, stop_current, num_points,steps_in_seq_,num_averages_,ro_dur_,ROIF1,ROIF2,deg_1,deg_2,qubit_1_thr,qubit_2_thr):
    sweep_vals = np.linspace(start_current, stop_current, num_points)
    steps_in_seq = steps_in_seq_
    num_steps = steps_in_seq_
    reps = num_averages_
    num_averages = num_averages_ # 10000
    ro_dur = ro_dur_ #RO duration
#    ro_freq = np.linspace(6.889,6.88888,num_points)
    out_I1, out_Q1, out_I2, out_Q2 = np.zeros( (4, num_points, steps_in_seq))
    for i,ss in enumerate(sweep_vals):
        keithley2401.set_current(ss, step_size_mA=0.01)
#        input_resonator_freq = 6.88888 #7.0137 #ro_freq[i]
#        bnc.set_bnc_output(input_resonator_freq, bnc_addr=target_bnc_qubit_address) #in GHz
#        keithley2401.set_voltage(ss, step_size_mV=10)
        
        #daq_params, rec_readout_vs_pats, p_vs_pats,I,Q,a,b,rec_readout,rec_avg_all,rec_all = daq_programs_homo.run_daq(steps_in_seq, num_averages, ro_dur)
        #above part was for homodyne
           
        n_vs_pats_1,n_vs_pats_2,rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1, deg_2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)

       
        out_Q1[i] = rec_avg_vs_pats_1[0];
        out_I1[i] = rec_avg_vs_pats_1[1];
        out_Q2[i] = rec_avg_vs_pats_2[0];
        out_I2[i] = rec_avg_vs_pats_2[1];
        plt.plot(out_Q1[i])
        plt.show()
        plt.plot(out_I1[i])
        plt.show()
        plt.plot(out_Q2[i])
        plt.show()
        plt.plot(out_I2[i])
        plt.show()
        print(i)
#        f1=-0.5#-0.05#
#        f2=0#-0.04#
#        freq = np.linspace(f1,f2,steps_in_seq)
#        Qrange = abs(np.max(Q)-np.min(Q))
#        Irange = abs(np.max(I)-np.min(I))
#        if Qrange>Irange:
#            freq_index = np.where(Q == np.amax(Q))
#        if Irange>Qrange:
#            freq_index = np.where(I == np.amax(I))
##    freq_index = np.where(y == np.max(y))
#        ssm_ge = freq[freq_index]
#        print(ssm_ge)
#        #time.sleep(4)
        
    # make plots
#    save_dir = r"Z:\candle qubit/"
#    np.savetxt(save_dir+'4_18_3pm_' + 'chevron_Q',out_Q)
    plt.imshow(out_Q1, extent=[0,steps_in_seq,stop_current,start_current],aspect='auto' )
    plt.show()
    plt.imshow(out_I1, extent=[0,steps_in_seq,stop_current,start_current],aspect='auto' )
    plt.show()
    plt.imshow(out_Q2, extent=[0,steps_in_seq,stop_current,start_current],aspect='auto' )
    plt.show()
    plt.imshow(out_I2, extent=[0,steps_in_seq,stop_current,start_current],aspect='auto' )
    plt.show()
    return out_Q1, out_I1,out_Q2, out_I2

def sweep_wx(start_amp, stop_amp, num_points):
    sweep_vals = np.linspace(start_amp, stop_amp, num_points)
    steps_in_seq = 81
    num_averages = 5000 # 10000
    ro_dur = 8000 #RO duration
    
    out_I, out_Q = np.zeros( (2, num_points, steps_in_seq))
    for i,ss in enumerate(sweep_vals):
        wx_programs.wx_set_and_amplitude_and_offset(amp=[ss,ss, .59,0.8])
        daq_params, rec_readout_vs_pats, p_vs_pats,I,Q,a,b = run_daq(steps_in_seq, num_averages, ro_dur)
                
        out_Q[i] = Q
        out_I[i] = I
        plt.plot(Q)
        plt.show()
        plt.plot(I)
        plt.show()
        print(i)
        #time.sleep(4)
        
    # make plots
#    save_dir = r"Z:\candle qubit/"
#    np.savetxt(save_dir+'4_18_3pm_' + 'chevron_Q',out_Q)
    plt.imshow(out_Q, extent=[0,steps_in_seq,stop_amp,start_amp],aspect='auto' )
    plt.show()
    plt.imshow(out_I, extent=[0,steps_in_seq,stop_amp,start_amp],aspect='auto' )
    return out_Q, out_I
#save_dir = r"C:\Data\2022\2022-04-08_SQUILL_nonHermitian/"

#sweep_keithley(0,0.08,1);np.savetxt(save_dir + 'keithley_sweep_0to0.08mA_Q',out_Q);np.savetxt(save_dir + 'keithley_sweep_0to0.08mA_I',out_I)
#for i in range(51):popt = fit_sine_decay(time,out_I[i],guess_vals=[10,10,2,1,157]);freq.append(popt[0][0]);gammas.append(popt[0][1])

##END sweep_ssb_freq
