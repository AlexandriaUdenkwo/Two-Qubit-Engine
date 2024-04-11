# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 14:47:59 2024

@author: crow104
"""

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
    
 if run_rabi_wx_sweep:
    num_steps = 81
    sweep_time =200 #ns
    Fs.rabi_ge(num_steps,sweep_time,ssm_ge,ROIF2,which_qubit)

if run_rabi_J: # 
    ef_amp = 0.01 #EP occurs at amp = 0.00227 
    num_steps = 81
    sweep_time =10000 #ns
    Fs.rabi_J(num_steps,sweep_time,ef_amp,pi_ge_time,pi_ef_time,ssm_ge,ssm_ef)
    
    
#    
if run_rabi_J:
    times = np.linspace(0,sweep_time/1000,num_steps)
    analysis.fit_sine_decay(times,y,guess_vals=[4,0.9,-0.4,-600,-100])#[10,2,2,1,-179]
    
    
    
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
    
if run_coupler_switch_sweep:
 
        sweep_steps_coupler_switch = 21
        
        coupler_switch_start = -0.5
        coupler_switch_stop = -1.5

        
        coupler_switch_steps = np.linspace(coupler_switch_start,coupler_switch_stop,sweep_steps_coupler_switch)
    
        
        guess_vals = [17,0.15,0.07,90,0.5]
    
        num_steps = 201
        sweep_time =200 #ns
        

        j=0
        out_I1, out_Q1, out_I2, out_Q2 = np.zeros( (4, sweep_steps_coupler_switch, num_steps))
        for j in range(len(coupler_switch_steps)):
            print('current coupler switch value:',coupler_switch_steps[j])
            print("j:",j)
            
            coupler_amp = coupler_switch_steps[j]
            if which_qubit == 0:
                Fs.vacuum_rabi(num_steps,sweep_time,ssm_geq2,ROIF1,ROIF2,which_qubit,pi_ge_time,coupler_amp=coupler_amp)
            if which_qubit == 1:
                Fs.vacuum_rabi(num_steps,sweep_time,ssm_geq1,ROIF1,ROIF2,which_qubit,pi_ge_time,coupler_amp=coupler_amp)

            
            wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs)
            n_vs_pats_1,n_vs_pats_2, rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
                                                                                                                                                                                                    

            out_Q1[j] = rec_avg_vs_pats_1[0];
            out_I1[j] = rec_avg_vs_pats_1[1];
            out_Q2[j] = rec_avg_vs_pats_2[0];
            out_I2[j] = rec_avg_vs_pats_2[1];
            plt.plot(out_Q1[j]);plt.show()
            plt.plot(out_I1[j]);plt.show()
            plt.plot(out_Q2[j]);plt.show()
            plt.plot(out_I2[j]);plt.show()

            print(j)
            
            plt.imshow(out_Q1, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()
            plt.imshow(out_I1, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()
            plt.imshow(out_Q2, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()
            plt.imshow(out_I2, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()


#        save_basename = '\vacuumrabiq2_5000_ns_bnc_4.32_ssmq1ge_'+str(ssm_geq1)+'_ssmq2ge_'+str(ssm_geq2)+'_wxch2amp_'+str(wx_start)+'_to_'+str(wx_stop)+'steps_'+str(steps_in_seq)+'_wxch2off_coupler'+str(wx_offs[1])+'_.txt'
        save_basename = '\sweep_vacuumrabi_coupler_switch_1.txt'
        np.savetxt(save_dir+save_basename+"_I1",out_I1)
        np.savetxt(save_dir+save_basename+"_I2",out_I2)
        np.savetxt(save_dir+save_basename+"_Q1",out_Q1)
        np.savetxt(save_dir+save_basename+"_Q2",out_Q2) 
        
if run_wx_coupler_switch_sweep:
        sweep_steps_amp = 31
        wx_amp_start = 0.5
        wx_amp_stop = 2.0
        
        sweep_steps_coupler_switch = 21
        coupler_switch_start = -0.5
        coupler_switch_stop = -1.5

        wx_ch2amp_steps = np.linspace(wx_amp_start,wx_amp_stop,sweep_steps_amp)
        coupler_switch_steps = np.linspace(coupler_switch_start,coupler_switch_stop,sweep_steps_coupler_switch)
    
        
    
        num_steps = 101
        sweep_time =200 #ns
        
        rabi_freq = np.zeros((sweep_steps_coupler_switch,sweep_steps_amp))
        rabi_amp = np.zeros((sweep_steps_coupler_switch,sweep_steps_amp))
        
        i=0
        j=0
        out_I1, out_Q1, out_I2, out_Q2, out_P1, out_P2 = np.zeros( (6, sweep_steps_amp, num_steps))
        for i in range(len(coupler_switch_steps)):
            print('current coupler switch value:',coupler_switch_steps[i])
            print("i:",i)
            
            coupler_amp = coupler_switch_steps[i]
            if which_qubit == 0:
                Fs.vacuum_rabi(num_steps,sweep_time,ssm_geq2,ROIF1,ROIF2,which_qubit,pi_ge_time,coupler_amp=coupler_amp)
            if which_qubit == 1:
                Fs.vacuum_rabi(num_steps,sweep_time,ssm_geq1,ROIF1,ROIF2,which_qubit,pi_ge_time,coupler_amp=coupler_amp)

            for j in range(len(wx_ch2amp_steps)):
                print('current coupler switch value:',coupler_switch_steps[i])
                print("i:",i)
                print('current amp:',wx_ch2amp_steps[j])
                print("amp j:",j)
    
                
                wx_amps[1] = wx_ch2amp_steps[j]
                wx.wx_set_and_amplitude_and_offset(amp=wx_amps,offset=wx_offs) #controls awg to give abitrary wave
            
                n_vs_pats_1,n_vs_pats_2, rec_avg_all, rec_all, rec_readout_1, rec_readout_2, rec_avg_vs_pats_1, rec_avg_vs_pats_2 , rec_all_het_1, rec_all_het_2, bins_1, bins_2, counts_1, counts_2,prob_vs_pats_1,prob_vs_pats_2,n_readout_1,n_readout_2 = daq_programs_homo.run_daq_het_2q(ROIF1,ROIF2, deg_1 = IQ_angle_q1, deg_2 = IQ_angle_q2,num_patterns=num_steps, num_records_per_pattern=reps,ro_dur=ro_dur,qubit_1_thr=qubit_1_thr,qubit_2_thr=qubit_2_thr, verbose=True)
                                                                                            
                out_Q1[j] = rec_avg_vs_pats_1[0];plt.plot(out_Q1[j]);plt.title('Q Q1');plt.show()
                out_I1[j] = rec_avg_vs_pats_1[1];plt.plot(out_I1[j]);plt.title('I Q1');plt.show()
                out_P1[j] = prob_vs_pats_1[0];plt.plot(out_P1[j]);plt.title('Q1 thresholded');plt.show()
                out_Q2[j] = rec_avg_vs_pats_2[0]; plt.plot(out_Q2[j]);plt.title('Q Q2');plt.show()
                out_I2[j] = rec_avg_vs_pats_2[1]; plt.plot(out_I2[j]);plt.title('I Q2');plt.show()
                out_P2[j] = prob_vs_pats_2[0];plt.plot(out_P2[j]);plt.title('Q2 thresholded');plt.show()
                
                
                if which_qubit == 0:
                    I = out_I2[j]; Q = out_Q2[j]
                    y2 = out_P2[j];
                if which_qubit == 1:
                    I = out_I1[j]; Q = out_Q1[j]
                    y1 = out_P1[j]
                    
                times = np.linspace(0,sweep_time/1000,num_steps)
                    
#                Qrange = abs(np.max(Q)-np.min(Q))
#                Irange = abs(np.max(I)-np.min(I))
#                if Qrange>Irange:
#                    pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,Q,guess_vals=[17,0.15,0.07,90,0.5])
#                if Irange>Qrange:
#                    pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,I,guess_vals=[17,0.15,0.07,90,0.5])
                 
                pi_ge_fit_vals,_,_,_ = analysis.fit_sine_decay(times,y2,guess_vals=[19,0.3,0.08,38,0.5])
                #guess_vals =  pi_ge_fit_vals
                pi_ge = abs((1/2/pi_ge_fit_vals[0])*1000)
                print("\u03C0_ge time = {} ns".format(pi_ge))
                
                rabi_freq[i][j] = pi_ge_fit_vals[0]
                rabi_amp[i][j] =  pi_ge_fit_vals[2]
                
               
                print(j)
                
                plt.imshow(out_Q1, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()
                plt.imshow(out_I1, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()
                plt.imshow(out_P1, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()
                plt.imshow(out_Q2, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()
                plt.imshow(out_I2, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()
                plt.imshow(out_P2, extent=[0,sweep_time,coupler_switch_stop,coupler_switch_start],aspect='auto' );plt.show()
                
            save_basename = '\sweep_vacuumrabi_coupler_switch_i_'+str(i)+'_file3.txt'
            np.savetxt(save_dir+save_basename+"_I1",out_I1)
            np.savetxt(save_dir+save_basename+"_I2",out_I2)
            np.savetxt(save_dir+save_basename+"_Q1",out_Q1)
            np.savetxt(save_dir+save_basename+"_Q2",out_Q2) 
            np.savetxt(save_dir+save_basename+"_P1",out_P1)
            np.savetxt(save_dir+save_basename+"_P2",out_P2)
        plt.imshow(rabi_freq,extent=[coupler_switch_start,coupler_switch_stop,wx_amp_start,wx_amp_stop],aspect='auto',cmap='gray', vmin=-20, vmax=20 );plt.show()
        plt.imshow(rabi_amp,extent=[coupler_switch_start,coupler_switch_stop,wx_amp_start,wx_amp_stop],aspect='auto',cmap='gray', vmin=-0.2, vmax=0.2 );plt.show()
        save_basename2 = '\sweep_vacuumrabi_characteristics_file3.txt'
        np.savetxt(save_dir+save_basename2+"_rabi_freq",rabi_freq)
        np.savetxt(save_dir+save_basename2+"_rabi_amp",rabi_amp)
  