# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 14:58:16 2023

@author: crow104
"""

class qubit1():
    self.ssb = -0.12
    self.ch = 3
    self.ro_pulse_dur=5000
    self.pi_ge_time =  34 #Keithley  #41#Hypres        #Keithley#33#41 #40.6 #37 #38 #q1
    wx_amps =[0.7,1.1,.55,.5] # [0.7,1.1,.55,.5]
    wx_offs = [0,0.01,0,0] # [0,0.2,0,0]
    IQ_angle =250 #225    
    RO_LO = 6.77#6.77#6.85 #6.6
    ROq1 = 6.872626 #6.872870#6.873020     
    ssm_geq1 =-0.1985 #-0.0767 #-0.0767#-0.0564  #   #-0.0796#Hypres    #-0.0796#Keithley                  #Keithley#-0.0564#-.2 #-0.077737#-.090+0.00025#-.126#-0.3275  
    ssm_efq1 = -.2375
    qubit_1_thr= [-1500,-400]#[-800,550]
    IQ_angle_q1 =  110 #130



class qubit2():
    ssm_efq2 = -.249
    ssm_geq2 = -0.1468
    self.ro_pulse_dur=5000
    ROq2 = 6.68862
    pi_ge_time = 44#50#49#50 #q2
    qubit_2_thr=[-2500,-1100]#[500,1700]#
    IQ_angle_q2 = 30
    self.ch = 4
    self.ro_dur=8000
    wx_amps =[0.7,1.1,.55,.5] # [0.7,1.1,.55,.5]
    wx_offs = [0,0.01,0,0] # [0,0.2,0,0]
    IQ_angle =250 #225    
    RO_LO = 6.77#6.77#6.85 #6.6
    ROq1 = 6.872626 #6.872870#6.873020  