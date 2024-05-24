#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 11:34:25 2024

@author: alexandriaudenkwo
"""
#Detunings
# 0.02 GHz
# 0.006 GHz
# 0.013 GHz

#Notations
#Q1 is QB
#Q2 is QA

import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.optimize import curve_fit
import numpy as np
import fit_functions as fitfun
import analysis as analy

#Directory where files are
importdir = r"/Users/alexandriaudenkwo/Desktop/Research/two_qubit_engine/"
#Directory to save data
savedir = r"/Users/alexandriaudenkwo/Desktop/Research/two_qubit_engine/13_mhz/"

#Read in file to array
data_list= [[]]*6
name_list = ['P1', 'P2', 'I1', 'I2', 'Q1', 'Q2']

thirteen = '13_mhz/wx_sweep4_vacuumrabi_500_detuning_13_mhz.txt_'
six = '6_mhz/wx_sweep1_rabi_500.txt_'
twenty = '20_mhz/wx_sweep2_vacuumrabi_500_detuning_20_mhz.txt_'

whichfile = thirteen
#denote as string to save files later on
whichfile_str = 'thirteen'

for i in range(6):
    data_list[i] = np.loadtxt(importdir+whichfile+ name_list[i], np.float)
  

P1 =data_list[0]; P2 = data_list[1]; I1 = data_list[2]; I2 = data_list[3]; Q1 = data_list[4]; Q2 = data_list[5];

plt.figure(1)
fig1 = plt.figure(1)
plt.imshow(P1, extent =[0, 500, 0.5, 2.0], aspect = 'auto' , origin = 'lower')
plt.xlabel("time (ns)")
plt.ylabel("FF line Amplitude (V)")
#fig1.savefig('ffl_amp.pdf') 



plt.figure(2)
plt.imshow(P2,extent =[0, 500, 0.5, 2.0], aspect = 'auto' , origin = 'lower')
plt.xlabel("time (ns)")
plt.ylabel("FF line Amplitude (V)")

#Raw Data shows ground state, we want excited state
y1= 1-P1
y2= 1-P2

#j denotes the index when oscillations start
j = 5

# np.savetxt(savedir+"thirteen_qb.txt", y1[j:])
# np.savetxt(savedir+"thirteen_qa.txt", y2[j:])



plt.figure(3)
plt.plot(1-P1[j], label = 'QB' );plt.plot(1-P2[j], label = 'QA'); plt.legend()


#Parameters
steps = len(P1)
wx_start = 0.5 
wx_stop = 2 
wx_steps = np.linspace(wx_start,wx_stop,steps) 

num_steps =101
sweep_time =500 #ns

###for 13 and 20 mhz
QA_min = np.zeros(steps)
QB_max = np.zeros(steps)

### for 6 mhz 
# QA_max = np.zeros(steps)
# QB_min = np.zeros(steps)


for i in range(steps):
    guess_amp_y1 = np.max(y1[i]) - np.min(y1[i])
    guess_amp_y2 = np.max(y2[i]) - np.min(y2[i])
    print(str(i))
    print("For QB")
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ge_fit_vals1,_,y1_vals_fit,_ = analy.fit_sine_decay(times,y1[i],guess_vals =[9,0.2,guess_amp_y1,0,y1[i][0]])
    
    ###for 13 and 20 mhz
    QB_max[i] = np.max(y1_vals_fit)
    print("QB_max = ", QB_max[i])
    
    ###for 6 mhz
    # QB_min[i] = np.min(y1_vals_fit)
    # print("QB_min = ", QB_min[i])
    
   
    print(str(i))
    print("For QA")
    times = np.linspace(0,sweep_time/1000,num_steps)
    pi_ge_fit_vals2,_,y2_vals_fit,_ = analy.fit_sine_decay(times,y2[i],guess_vals=[9,0.3,guess_amp_y2,0,y2[i][0]])
    
    
    # ###for 13 and 20 mhz
    QA_min[i] = np.min(y2_vals_fit)
    print("QA_min = ", QA_min[i])
    
    ###for 6 mhz
    # QA_max[i] = np.max(y2_vals_fit)
    # print("QA_max = ", QA_max[i])

    #Grab example of vaccuum Rabi oscillations
    if i == 20:
        plt.figure(6)
        plt.plot(times*1000, y1[i], label = "P(01)"); plt.plot(times*1000, y2[i], label = "P(10)"); plt.legend(); plt.title("Raw data")
        plt.figure(7)
        current_wx_val = wx_start +  i*((wx_stop - wx_start)/steps)
        plt.plot(times*1000, y1_vals_fit, label = "P(01)"); plt.plot(times*1000, y2_vals_fit, label = "P(10)")
        plt.legend(); plt.title("Vacuum Rabi Oscillations for Detuning 13 MHz, wx amp = "+str(current_wx_val))
        plt.xlabel(xlabel = "time (ns)")
 
#Calculate theta          
theta_qb = np.zeros(steps)
theta_qa = np.zeros(steps)
for i in range(j, steps):
    
    # ###for 13 and 20 mhz
    theta_qb[i] = np.arcsin(np.sqrt(QB_max[i]))
    theta_qa[i] = np.arccos(np.sqrt(QA_min[i]))
    
    ### for 6 mhz
    # theta_qb[i] = np.arccos(np.sqrt(QB_min[i]))
    # theta_qa[i] = np.arcsin(np.sqrt(QA_max[i]))

#Plot theta vs wx_amplitude
plt.figure(8)
plt.plot(wx_steps[j:], theta_qb[j:], label = r"$\theta$"+"_QB"); plt.plot(wx_steps[j:], theta_qa[j:], label = r"$\theta$"+"_QA");plt.legend()
plt.title(label = "Detuning: 13 MHz")
plt.xlabel("wx amplitude")
plt.ylabel(r"$\theta$")


#to deal with formating issues in save file
new_wx_start = wx_start + j*((wx_stop - wx_start)/steps)
dummy_wx_str = str(new_wx_start)
for i in range(len(dummy_wx_str)):
    if "." == dummy_wx_str[i]:
        wx_str = dummy_wx_str[:i] + 'p' + dummy_wx_str[i+1:]
#save
# np.savetxt(savedir+whichfile_str+"wx_amplitude_steps_"+ wx_str+"to2p0.txt", wx_steps[j:])
# np.savetxt(savedir+whichfile_str+"theta_qb.txt",theta_qb[j:] )
# np.savetxt(savedir+whichfile_str+"theta_qa.txt",theta_qa[j:] )

#Calculate tan(theta)
tan_qb = np.zeros(steps)
tan_qa = np.zeros(steps)
for i in range(j, steps):
    
    tan_qb[i] = np.tan(theta_qb[i])
    tan_qa[i] = np.tan(theta_qa[i])
    
plt.figure(9)
plt.plot(wx_steps[j:], tan_qb[j:], label = "QB"); plt.plot(wx_steps[j:], tan_qa[j:], label = "QA");plt.legend() 
plt.title(label = "Detuning: 13 MHz")   
plt.xlabel("wx amplitude")
plt.ylabel("g/"+"$\delta$")

#save
# np.savetxt(savedir+whichfile_str+"tan_qb.txt",tan_qb[j:] )
# np.savetxt(savedir+whichfile_str+"tan_qa.txt",tan_qa[j:] )


