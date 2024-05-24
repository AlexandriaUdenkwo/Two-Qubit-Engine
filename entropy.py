#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: alexandriaudenkwo
"""

import matplotlib.pyplot as plt
import numpy as np


#Directory where files are
savedir = r"/Users/alexandriaudenkwo/Desktop/Research/two_qubit_engine/"

thirteen = '13_mhz/'
twenty = '20_mhz/'

#let hbar =1 mhz 

#Load analyzed data
wxsteps = np.loadtxt(thirteen+ 'thirteenwx_amplitude_steps_0p8to2p0.txt',np.float)

theta_qb_13 = np.loadtxt(savedir+thirteen+ 'thirteentheta_qb.txt',np.float)
theta_qa_13 = np.loadtxt(savedir+thirteen+ 'thirteentheta_qa.txt',np.float)
theta_qb_20 = np.loadtxt(savedir+twenty+ 'twentytheta_qb.txt',np.float)
theta_qa_20 = np.loadtxt(savedir+twenty+ 'twentytheta_qa.txt',np.float)


tan_qb_13 = np.loadtxt(savedir+thirteen+ 'thirteentan_qb.txt',np.float)
tan_qa_13 = np.loadtxt(savedir+thirteen+ 'thirteentan_qa.txt',np.float)
tan_qb_20 = np.loadtxt(savedir+twenty+ 'twentytan_qb.txt',np.float)
tan_qa_20 = np.loadtxt(savedir+twenty+ 'twentytan_qa.txt',np.float)


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
    
    
#Plot entropy vs g/delta    
plt.figure(1)
fig1 = plt.figure(1)
plt.plot(tan_qb_20[:len(wxsteps)-1], entropy_qb_20[:len(wxsteps)-1], label = r"$\delta$"+" = 20 MHz"); #plt.plot(wxsteps, entropy_qa_20, label = r"QA 20 MHz")
plt.plot(tan_qb_13[:len(wxsteps)-1], entropy_qb_13[:len(wxsteps)-1], label = r"$\delta$"+r" = 13 MHz"); plt.legend() #plt.plot(wxsteps, entropy_qa_13, label = r"QA 13 MHz");plt.legend()
plt.title(label = "S"+"$^{meas}$"+" vs g/"+r"$\delta$")
plt.xlabel("g/"+r"$\delta$")
plt.ylabel("S"+"$^{meas}$")
#fig1.savefig('smeas_vs_gdelta_paper.pdf')

#Let's look at rabi oscillations at the max entropy and near the max
print("max at (for detuning 13 mhz): ", entropy_qb_13.argmax())
print("max at (for detuning 20 mhz): ", entropy_qb_20.argmax())
    
cc = entropy_qb_13.argmax()
dd = entropy_qb_20.argmax()
#import rabi oscillations at max
thirteen_file_qb = np.loadtxt(savedir+thirteen+ 'thirteen_qb.txt',np.float)
thirteen_file_qa = np.loadtxt(savedir+thirteen+ 'thirteen_qa.txt',np.float)
timestep = np.linspace(0,500, 101)


plt.figure(3)
fig3 = plt.figure(3)
plt.plot(timestep, thirteen_file_qb[cc],label = "P(01)", color = 'blue');plt.plot(timestep,thirteen_file_qa[cc], label ="P(10)", color = 'black');plt.legend()
plt.title("Rabi Oscillations at " +"S"+"$^{meas}$" + "= 0.99, "+"g/"+r"$\delta$"+"= 1.0")
plt.xlabel("time (ns)")
#fig3.savefig('rabi_touch_paper.pdf')

plt.figure(4)
fig4 = plt.figure(4)
plt.plot(timestep,thirteen_file_qb[cc+2],label = "P(01)", color = 'blue');plt.plot(timestep,thirteen_file_qa[cc+2], label ="P(10)", color = 'black');plt.legend()
plt.title("Rabi Oscillations at " +"S"+"$^{meas}$" + "= 0.93, "+"g/"+r"$\delta$"+"= 1.4")
plt.xlabel("time (ns)")
#fig4.savefig('rabi_intersect_paper.pdf')

# plt.figure(5)
# plt.plot(timestep, thirteen_file_qb[cc-2],label = "P(01)", color = 'blue');plt.plot(timestep, thirteen_file_qa[cc-2], label ="P(10)", color = 'black')
# plt.title("Rabi Oscillations at " +"S"+"$^{meas}$" + "= 0.92, " +"g/"+r"$\delta$"+"= 0.7");plt.legend()
# plt.xlabel("time (ns)")

# plt.figure(6)
# plt.plot(timestep, thirteen_file_qb[10],label = "P(01)", color = 'blue');plt.plot(timestep, thirteen_file_qa[10], label ="P(10)", color = 'black');plt.legend()
# plt.title("Rabi Oscillations at " +r"$\delta$"+"= 13 MHz, g = 9.2 MHz")
# plt.xlabel("time (ns)")


# print(len(wxsteps))
# print(wxsteps[cc+4])



# plt.figure(7)
# plt.plot(1/tan_qb_20[:len(wxsteps)-1], entropy_qb_20[:len(wxsteps)-1], label = "QB 20 MHz"); #plt.plot(wxsteps, entropy_qa_20, label = r"QA 20 MHz")
# plt.plot(1/tan_qb_13[:len(wxsteps)-1], entropy_qb_13[:len(wxsteps)-1], label = r"QB 13 MHz"); plt.legend() #plt.plot(wxsteps, entropy_qa_13, label = r"QA 13 MHz");plt.legend()
# plt.title(label = "S"+"$^{meas}$ "+"vs "+r"$\delta$"+"/g")
# plt.xlabel(r"$\delta$"+"/g")
# plt.ylabel("S"+"$^{meas}$")