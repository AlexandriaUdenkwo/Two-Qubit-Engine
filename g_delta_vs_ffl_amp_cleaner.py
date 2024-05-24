#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 07:28:30 2024

@author: alexandriaudenkwo
"""

import matplotlib.pyplot as plt
import numpy as np


#Directory where files are
savedir = r"/Users/alexandriaudenkwo/Desktop/Research/two_qubit_engine/"

thirteen = '13_mhz/'
twenty = '20_mhz/'

#Steps are all the same so doesn't matter which detuning or qubit is loaded
wxsteps = np.loadtxt(savedir+thirteen+ 'thirteenwx_amplitude_steps_0p8to2p0.txt',np.float)

#Load data for theta
theta_qb_13 = np.loadtxt(savedir+thirteen+ 'thirteentheta_qb.txt',np.float)
theta_qa_13 = np.loadtxt(savedir+thirteen+ 'thirteentheta_qa.txt',np.float)
theta_qb_20 = np.loadtxt(savedir+twenty+ 'twentytheta_qb.txt',np.float)
theta_qa_20 = np.loadtxt(savedir+twenty+ 'twentytheta_qa.txt',np.float)

plt.figure(1)
plt.plot(wxsteps[:len(wxsteps)-1], theta_qb_20[:len(wxsteps)-1], label = r"$\delta$"+"= 20 MHz"); #plt.plot(wxsteps, theta_qa_20, label = r"QA 20 MHz")
plt.plot(wxsteps[:len(wxsteps)-1], theta_qb_13[:len(wxsteps)-1], label = r"$\delta$"+r"= 13 MHz");plt.legend() #plt.plot(wxsteps, theta_qa_13, label = r"QA 13 MHz")
plt.title(label = r"$\theta$"+" vs wx amplitude setting")
plt.xlabel("wx amplitude")
plt.ylabel(r"$\theta$")

#Load data for tan(theta) = g/delta
tan_qb_13 = np.loadtxt(savedir+thirteen+ 'thirteentan_qb.txt',np.float)
tan_qa_13 = np.loadtxt(savedir+thirteen+ 'thirteentan_qa.txt',np.float)
tan_qb_20 = np.loadtxt(savedir+twenty+ 'twentytan_qb.txt',np.float)
tan_qa_20 = np.loadtxt(savedir+twenty+ 'twentytan_qa.txt',np.float)

plt.figure(2)
fig2 = plt.figure(2)
plt.plot(wxsteps[:len(wxsteps)-1], tan_qb_20[:len(wxsteps)-1], label = r"$\delta$"+"= 20 MHz"); #plt.plot(wxsteps, tan_qa_20, label = "QA 20 MHz")
plt.plot(wxsteps[:len(wxsteps)-1], tan_qb_13[:len(wxsteps)-1], label = r"$\delta$"+r"= 13 MHz"); plt.legend()#plt.plot(wxsteps, tan_qa_13, label = "QA 13 MHz");plt.legend()
plt.title(label = "g/"+r"$\delta$"+" vs FF line amplitude")
plt.xlabel("FF line amplitude (V)")
plt.ylabel("g/"+r"$\delta$")
#fig2.savefig('gdelta_vs_ffl_V_paper.pdf')
