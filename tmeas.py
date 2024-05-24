#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 08:33:27 2024

@author: alexandriaudenkwo
"""

import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.optimize import curve_fit
import numpy as np


#Directory where files are
savedir = r"/Users/alexandriaudenkwo/Desktop/Research/two_qubit_engine/"

thirteen = '13_mhz/'
twenty = '20_mhz/'

#let hbar =1 mhz 
#Steps are all the same so doesn't matter which detuning or qubit is loaded
wxsteps = np.loadtxt(savedir+thirteen+ 'thirteenwx_amplitude_steps_0p8to2p0.txt',np.float)

#Load data for theta
theta_qb_13 = np.loadtxt(savedir+thirteen+ 'thirteentheta_qb.txt',np.float)
theta_qa_13 = np.loadtxt(savedir+thirteen+ 'thirteentheta_qa.txt',np.float)
theta_qb_20 = np.loadtxt(savedir+twenty+ 'twentytheta_qb.txt',np.float)
theta_qa_20 = np.loadtxt(savedir+twenty+ 'twentytheta_qa.txt',np.float)

 #Load data for tan(theta) = g/delta   
tan_qb_13 = np.loadtxt(thirteen+ 'thirteentan_qb.txt',np.float)
tan_qa_13 = np.loadtxt(thirteen+ 'thirteentan_qa.txt',np.float)
tan_qb_20 = np.loadtxt(twenty+ 'twentytan_qb.txt',np.float)
tan_qa_20 = np.loadtxt(twenty+ 'twentytan_qa.txt',np.float)



#Intialize arrays to calculate Tmeas ratio
entropy_qb_13 = np.zeros(len(wxsteps))
entropy_qb_20 = np.zeros(len(wxsteps))
entropy_qa_13 = np.zeros(len(wxsteps))
entropy_qa_20 = np.zeros(len(wxsteps))

energy_qb_13 = np.zeros(len(wxsteps))
energy_qb_20 = np.zeros(len(wxsteps))
energy_qa_13 = np.zeros(len(wxsteps))
energy_qa_20 = np.zeros(len(wxsteps))

t_qb_13 = np.zeros(len(wxsteps))
t_qb_20 = np.zeros(len(wxsteps))
t_qa_13 = np.zeros(len(wxsteps))
t_qa_20 = np.zeros(len(wxsteps))

#Calculate tmeas ratio
for i in range(len(wxsteps)):
    entropy_qb_13[i] = -((np.cos(theta_qb_13[i]))**2)*np.log2((np.cos(theta_qb_13[i]))**2) -((np.sin(theta_qb_13[i]))**2)*np.log2((np.sin(theta_qb_13[i]))**2)
    entropy_qb_20[i] = -((np.cos(theta_qb_20[i]))**2)*np.log2((np.cos(theta_qb_20[i]))**2) -((np.sin(theta_qb_20[i]))**2)*np.log2((np.sin(theta_qb_20[i]))**2)

    entropy_qa_13[i] = -((np.cos(theta_qa_13[i]))**2)*np.log2((np.cos(theta_qa_13[i]))**2) -((np.sin(theta_qa_13[i]))**2)*np.log2((np.sin(theta_qa_13[i]))**2)
    entropy_qa_20[i] = -((np.cos(theta_qa_20[i]))**2)*np.log2((np.cos(theta_qa_20[i]))**2) -((np.sin(theta_qa_20[i]))**2)*np.log2((np.sin(theta_qa_20[i]))**2)

    energy_qb_13[i] = 13*np.sin(theta_qb_13[i])**2
    energy_qb_20[i] = 20*np.sin(theta_qb_20[i])**2
    
    energy_qa_13[i] = 13*np.sin(theta_qa_13[i])**2
    energy_qa_20[i] = 20*np.sin(theta_qa_20[i])**2
    
    t_qb_13[i] = energy_qb_13[i]/entropy_qb_13[i]
    t_qb_20[i] = energy_qb_20[i]/entropy_qb_20[i]
    t_qa_13[i] = energy_qa_13[i]/entropy_qa_13[i]
    t_qa_20[i] = energy_qa_20[i]/entropy_qa_20[i]


#Plot Tmeas
plt.figure(1)
fig1 = plt.figure(1)
plt.plot(tan_qb_20[:len(wxsteps)-1], t_qb_20[:len(wxsteps)-1], label = r"$\delta$ = 20 MHz"); #plt.plot(wxsteps, t_qa_20, label = r"QA 20 MHz")
plt.plot(tan_qb_13[:len(wxsteps)-1], t_qb_13[:len(wxsteps)-1], label = r"$\delta$ = 13 MHz"); plt.legend() #plt.plot(wxsteps, t_qa_13, label = r"QA 13 MHz");plt.legend()
plt.title(label = "$T^{meas}$ vs g/"+r"$\delta$")
plt.xlabel("g/"+r"$\delta$")
plt.ylabel("$T^{meas}$ (MHz)") 
#fig1.savefig('tmeas_vs_gdelta_paper.pdf')   
    