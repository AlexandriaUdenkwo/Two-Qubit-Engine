#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 07:48:31 2024

@author: alexandriaudenkwo
"""

import matplotlib.pyplot as plt
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

#Load data for tan(theta)
tan_qb_13 = np.loadtxt(savedir+thirteen+ 'thirteentan_qb.txt',np.float)
tan_qa_13 = np.loadtxt(savedir+thirteen+ 'thirteentan_qa.txt',np.float)
tan_qb_20 = np.loadtxt(savedir+twenty+ 'twentytan_qb.txt',np.float)
tan_qa_20 = np.loadtxt(savedir+twenty+ 'twentytan_qa.txt',np.float)

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


#Plot work vs g/delta
plt.figure(1)
fig1 = plt.figure(1)
plt.figure(figsize=(3,1.5))
plt.plot(tan_qb_20[:len(wxsteps)-1], 20*np.sin(theta_qb_20[:len(wxsteps)-1])**2, label = r"$\delta$"+" = 20 MHz"); #plt.plot(wxsteps, 20*np.sin(theta_qa_20)**2, label = r"QA 20 MHz")
plt.plot(tan_qb_13[:len(wxsteps)-1], 13*np.sin(theta_qb_13[:len(wxsteps)-1])**2, label = r"$\delta$"+r" = 13 MHz"); #plt.plot(wxsteps, 13*np.sin(theta_qa_13)**2, label = r"QA 13 MHz")
plt.plot(tan_qb_13[:len(wxsteps)-1], work13[:len(wxsteps)-1], label = "W = 13 MHz", linestyle = "dashed");plt.plot(tan_qb_13[:len(wxsteps)-1], work20[:len(wxsteps)-1], label = "W = 20 MHz", linestyle = "dashed")

plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")


plt.title(label = "$E^{meas}$ vs g/"+r"$\delta$")
plt.xlim(0, 5)
plt.ylim(0,21)
plt.xlabel("g/"+r"$\delta$")
plt.ylabel("$E^{meas}$ (MHz)")
fig1.savefig('emeas_vs_gdelta_paper.pdf') 

#Below I played around with plotting work vs other parameters

# plt.figure(2)
# plt.plot((20/(4.1462*1000))*tan_qb_20[:len(wxsteps)-1], (20/(4.1462*1000))*np.sin(theta_qb_20[:len(wxsteps)-1])**2, label = r"$\delta$/"+r"$\omega_A$"+" = 0.003 MHz"); #plt.plot(wxsteps, 20*np.sin(theta_qa_20)**2, label = r"QA 20 MHz")
# plt.plot((13/(4.1462*1000))*tan_qb_13[:len(wxsteps)-1], (13/(4.1462*1000))*np.sin(theta_qb_13[:len(wxsteps)-1])**2, label = r"$\delta$/"+r"$\omega_A$"+r" = 0.005 MHz"); #plt.plot(wxsteps, 13*np.sin(theta_qa_13)**2, label = r"QA 13 MHz")
# plt.plot((13/(4.1462*1000))*tan_qb_13[:len(wxsteps)-1], work13_fig2[:len(wxsteps)-1], label = "W = 13 MHz", linestyle = "dashed");plt.plot((13/(4.1462*1000))*tan_qb_13[:len(wxsteps)-1], work20_fig2[:len(wxsteps)-1], label = "W = 20 MHz", linestyle = "dashed");plt.legend()
# plt.title(label = "$E^{meas}$ vs g/"+r"$\omega_A$")
# plt.xlabel("g/"+r"$\omega_A$")
# plt.ylabel("$E^{meas}$/" +r"$\omega_A$"+" (MHz)")

# plt.figure(3)
# plt.plot((20)*tan_qb_20[:len(wxsteps)-1], (20)*np.sin(theta_qb_20[:len(wxsteps)-1])**2, label = r"$\delta$"+" = 13 MHz"); #plt.plot(wxsteps, 20*np.sin(theta_qa_20)**2, label = r"QA 20 MHz")
# plt.plot((13)*tan_qb_13[:len(wxsteps)-1], (13)*np.sin(theta_qb_13[:len(wxsteps)-1])**2, label = r"$\delta$"+r" = 20 MHz"); #plt.plot(wxsteps, 13*np.sin(theta_qa_13)**2, label = r"QA 13 MHz")
# plt.plot((13)*tan_qb_13[:len(wxsteps)-1], work13[:len(wxsteps)-1], label = "W = 13 MHz", linestyle = "dashed");plt.plot((13)*tan_qb_13[:len(wxsteps)-1], work20[:len(wxsteps)-1], label = "W = 20 MHz", linestyle = "dashed");plt.legend()
# plt.title(label = "$E^{meas}$ vs g")
# plt.xlabel("g (MHz)")
# plt.ylabel("$E^{meas}$ (MHz)")

# plt.figure(4)
# plt.plot(1/tan_qb_20[:len(wxsteps)-1], 20*np.sin(theta_qb_20[:len(wxsteps)-1])**2, label = r"$\delta$"+" = 20 MHz"); #plt.plot(wxsteps, 20*np.sin(theta_qa_20)**2, label = r"QA 20 MHz")
# plt.plot(1/tan_qb_13[:len(wxsteps)-1], 13*np.sin(theta_qb_13[:len(wxsteps)-1])**2, label = r"$\delta$"+r" = 13 MHz"); #plt.plot(wxsteps, 13*np.sin(theta_qa_13)**2, label = r"QA 13 MHz")
# plt.plot(1/tan_qb_13[:len(wxsteps)-1], work13[:len(wxsteps)-1], label = "W = 13 MHz", linestyle = "dashed");plt.plot(1/tan_qb_13[:len(wxsteps)-1], work20[:len(wxsteps)-1], label = "W = 20 MHz", linestyle = "dashed");plt.legend()


# plt.title(label = "$E^{meas}$ vs "+r"$\delta$" +"/g")
# plt.xlabel(r"$\delta$"+"/g")
# plt.ylabel("$E^{meas}$ (MHz)")

# plt.figure(5)
# plt.plot((20/(4.1462*1000))*(1/tan_qb_20[:len(wxsteps)-1]), (20/(4.1462*1000))*np.sin(theta_qb_20[:len(wxsteps)-1])**2, label = r"$\delta$"+" = 20 MHz"); #plt.plot(wxsteps, 20*np.sin(theta_qa_20)**2, label = r"QA 20 MHz")
# plt.plot((13/(4.1462*1000))*(1/tan_qb_13[:len(wxsteps)-1]), (13/(4.1462*1000))*np.sin(theta_qb_13[:len(wxsteps)-1])**2, label = r"$\delta$"+r" = 13 MHz"); #plt.plot(wxsteps, 13*np.sin(theta_qa_13)**2, label = r"QA 13 MHz")
# plt.plot((13/(4.1462*1000))*(1/tan_qb_13[:len(wxsteps)-1]), work13_fig2[:len(wxsteps)-1], label = "W = 13 MHz", linestyle = "dashed");plt.plot((13/(4.1462*1000))*(1/tan_qb_13[:len(wxsteps)-1]), work20_fig2[:len(wxsteps)-1], label = "W = 20 MHz", linestyle = "dashed");plt.legend()


# plt.title(label = "$E^{meas}$ vs "+r"$\delta$" +"/g")
# plt.xlabel(r"$\delta$"+"/g")
# plt.ylabel("$E^{meas}$ (MHz)")

