#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 07:48:31 2024

@author: alexandriaudenkwo
"""

import matplotlib.pyplot as plt
import numpy as np


#Directory where files are
# savedir = r"/Users/alexandriaudenkwo/Desktop/Research/two_qubit_engine/"

#Directory where files are
# importdir = r"/Users/alexandriaudenkwo/Desktop/Research/two_qubit_engine/"
importdir = r"/Users/katermurch/Documents/Python/entanglement_engine_data/"
#Directory to save data
savedir = r"/Users/katermurch/Documents/Python/entanglement_engine_data/"

# thirteen = 'detuning_13_mhz/'
# twenty = 'detuning_20_mhz/'
thirteen = ''
twenty = ''

#let hbar =1 mhz 

#Steps are all the same so doesn't matter which detuning or qubit is loaded
wxsteps = np.loadtxt(savedir+thirteen+ 'thirteenwx_amplitude_steps_0p8to2p0.txt',float)

#Load data for theta
theta_qb_13 = np.loadtxt(savedir+thirteen+ 'thirteentheta_qb.txt',float)
theta_qa_13 = np.loadtxt(savedir+thirteen+ 'thirteentheta_qa.txt',float)
theta_qb_20 = np.loadtxt(savedir+twenty+ 'twentytheta_qb.txt',float)
theta_qa_20 = np.loadtxt(savedir+twenty+ 'twentytheta_qa.txt',float)

#Load data for tan(theta)
tan_qb_13 = np.loadtxt(savedir+thirteen+ 'thirteentan_qb.txt',float)
tan_qa_13 = np.loadtxt(savedir+thirteen+ 'thirteentan_qa.txt',float)
tan_qb_20 = np.loadtxt(savedir+twenty+ 'twentytan_qb.txt',float)
tan_qa_20 = np.loadtxt(savedir+twenty+ 'twentytan_qa.txt',float)

#Work bound
work13 = np.zeros(len(wxsteps))
work20 = np.zeros(len(wxsteps))
work13_fig2 = np.zeros(len(wxsteps))
work20_fig2 = np.zeros(len(wxsteps))
for i in range(len(wxsteps)):
    work13[i] = 13
    work20[i] = 20
    work13_fig2[i] = (13/(4.1462*1000))
    work20_fig2[i] = (20/(4.1462*1000))

#%%
#Updated by Kater with ChatGPT
#Plot work vs g/delta
plt.figure(figsize=(3, 1.5))

# Plotting the data
plt.plot(tan_qb_20[:len(wxsteps)-1], 20*np.sin(theta_qb_20[:len(wxsteps)-1])**2, label=r"$\delta/2\pi$ = 20 MHz", color='black')
plt.plot(tan_qb_13[:len(wxsteps)-1], 13*np.sin(theta_qb_13[:len(wxsteps)-1])**2, label=r"$\delta/2\pi$ = 13 MHz", color ='purple')
plt.plot(tan_qb_13[:len(wxsteps)-1], work13[:len(wxsteps)-1],  linestyle="dashed", color ='purple')
plt.plot(tan_qb_13[:len(wxsteps)-1], work20[:len(wxsteps)-1],  linestyle="dashed", color = 'black')

# Adding legend and labels
plt.legend( loc="lower right")
#plt.title(label="$E^{meas}$ vs g/"+r"$\delta$")
plt.xlim(0, 5)
plt.ylim(0, 21)
plt.xlabel("g/"+r"$\delta$")
plt.ylabel("$E^{\mathrmmeas}$ (MHz)")

# Saving the figure
plt.savefig('emeas_vs_gdelta_paper.pdf')

# Displaying the plot
plt.show()

#%%

#Initialize entropy arrays for each detuning and qubit
entropy_qb_13 = np.zeros(len(wxsteps))
entropy_qb_20 = np.zeros(len(wxsteps))
entropy_qa_13 = np.zeros(len(wxsteps))
entropy_qa_20 = np.zeros(len(wxsteps))

#Calculate entropy
for i in range(len(wxsteps)):
    entropy_qb_13[i] = -((np.cos(theta_qb_13[i]))**2)*np.log2((np.cos(theta_qb_13[i]))**2) -((np.sin(theta_qb_13[i]))**2)*np.log2((np.sin(theta_qb_13[i]))**2)
    entropy_qb_20[i] = -((np.cos(theta_qb_20[i]))**2)*np.log2((np.cos(theta_qb_20[i]))**2) -((np.sin(theta_qb_20[i]))**2)*np.log2((np.sin(theta_qb_20[i]))**2)

    entropy_qa_13[i] = -((np.cos(theta_qa_13[i]))**2)*np.log2((np.cos(theta_qa_13[i]))**2) -((np.sin(theta_qa_13[i]))**2)*np.log2((np.sin(theta_qa_13[i]))**2)
    entropy_qa_20[i] = -((np.cos(theta_qa_20[i]))**2)*np.log2((np.cos(theta_qa_20[i]))**2) -((np.sin(theta_qa_20[i]))**2)*np.log2((np.sin(theta_qa_20[i]))**2)
    
#%%
#Plot entropy vs g/delta    

plt.figure(figsize=(3,1.5))
plt.plot(tan_qb_13[:len(wxsteps)-1], entropy_qb_13[:len(wxsteps)-1], label = r"$\delta$"+r" = 13 MHz", color='purple')
plt.plot(tan_qb_20[:len(wxsteps)-1], entropy_qb_20[:len(wxsteps)-1], label = r"$\delta$"+" = 20 MHz", color = 'black'); #plt.plot(wxsteps, entropy_qa_20, label = r"QA 20 MHz")


plt.legend()
plt.xlim(0, 5)
plt.ylim(0,1.1)
plt.xlabel("g/"+r"$\delta$")
plt.ylabel("S"+"$^\mathrm{meas}$"+" (bits)")

# Saving the figure
plt.savefig('s_vs_gdelta_paper.pdf')

# Displaying the plot
plt.show()

#%%
#Plot efficiency vs g/delta    
t_meas_13 = work13/entropy_qb_13
t_meas_20 = work20/entropy_qb_20

plt.figure(figsize=(3,1.5))
plt.plot(tan_qb_13[:len(wxsteps)-1], t_meas_13[:len(wxsteps)-1], label = r"$\delta$"+r" = 13 MHz", color='purple')
plt.plot(tan_qb_20[:len(wxsteps)-1], t_meas_20[:len(wxsteps)-1], label = r"$\delta$"+r" = 20 MHz", color='black')


plt.legend()
# plt.xlim(0, 5)
# plt.ylim(0,1.1)
plt.xlabel("g/"+r"$\delta$")
plt.ylabel("$T$"+"$^\mathrm{meas}$")

# Saving the figure
plt.savefig('t_vs_gdelta_paper.pdf')

# Displaying the plot
plt.show()

#Let's look at rabi oscillations at the max entropy and near the max

# plt.figure(2)
# fig3 = plt.figure(2)
# plt.figure(figsize=(3,1.5))
# plt.plot(timestep, thirteen_file_qb[cc],label = "P(01)", color = 'blue');plt.plot(timestep,thirteen_file_qa[cc], label ="P(10)", color = 'black')
# plt.legend()
# #plt.title("Rabi Oscillations at " +"S"+"$^{meas}$" + "= 0.99, "+"g/"+r"$\delta$"+"= 1.0")
# plt.xlabel("time (ns)")
# #fig3.savefig('rabi_touch_paper.pdf')

# plt.figure(3)
# fig4 = plt.figure(3)
# plt.figure(figsize=(3,1.5))
# plt.plot(timestep,thirteen_file_qb[cc+2],label = "P(01)", color = 'blue');plt.plot(timestep,thirteen_file_qa[cc+2], label ="P(10)", color = 'black')
# plt.legend()
# #plt.title("Rabi Oscillations at " +"S"+"$^{meas}$" + "= 0.93, "+"g/"+r"$\delta$"+"= 1.4")
# plt.xlabel("time (ns)")
# #fig4.savefig('rabi_intersect_paper.pdf')

#%%