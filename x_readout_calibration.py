import numpy as np
import sys
file_dir= r"C:\Users\crow104\Documents\Python Scripts\sequence_generator"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Users\crow104\Documents\Python Scripts\data_acquisition_scripts"
if file_dir not in sys.path: sys.path.append(file_dir)
file_dir = r"C:\Data\2020\201021_thermal_readout\Ramsey vs thermal noise"
if file_dir not in sys.path: sys.path.append(file_dir)
from generator import *
import os
pi = np.pi
import pyvisa
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
import daq_programs as dp
from tkinter import Tk
import tkinter as tki
from tkinter.filedialog import askopenfilename, asksaveasfilename, asksaveasfile
import wx_programs
import scipy.io as io
import data_analysis as da
import time
from IPython import get_ipython
import data_analysis as da
import thermal_readout as tr
from scipy import optimize as opt
from matplotlib.patches import Rectangle


 # Readout Integration Parameters 
startindex = 300
 
#angle = 205 # for coherent readout
#duration = 3000 # for coherent readout
#angle = 100 # for high-power readout JPA 0 dB
angle = 229 # for high-power readout JPA 20 dB
duration = 1000 # for high-power readout

#angle = 150 # for thermal readout
#duration = 3000 # for thermal readout

#angle = 168 # for thermal readout
#duration = 3000 # for FT readout


endindex = startindex + duration
bgstartindex = 5000
bgendindex = bgstartindex + duration


# cutoff value to determine qubit in g or e state (run compare_ge_readout to get this value)
cutoff = 7
#cutoff = 17


def compare_ge_readout(numrec=1000):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    tr.check_readout(state = 0, measure = 1, num_steps = 5)
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=5, num_records_per_pattern=numrec))
    
    
    rec_rot,_ = process_readout_raw(rec_all_raw_ave, ifplot=1)
    
    ch_a_g = rec_rot[0]
    ch_b_g = rec_rot[1]
    
    tr.check_readout(state = 1, measure = 1, num_steps = 5)
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=5, num_records_per_pattern=numrec))
    
    
    rec_rot,_ = process_readout_raw(rec_all_raw_ave, ifplot=1)
    
    ch_a_e = rec_rot[0]
    ch_b_e = rec_rot[1]
    
    
    plt.figure()
    ax_scatter = plt.axes()
    ax_scatter.scatter(ch_a_g, ch_b_g, s=0.5, c='b', marker='.', alpha=0.2)
    ax_scatter.scatter(ch_a_e, ch_b_e, s=0.5, c='r', marker='.', alpha=0.2)
    
    lim = np.max(np.concatenate((np.abs(ch_a_g), np.abs(ch_a_e),
                                 np.abs(ch_b_g), np.abs(ch_b_e))))*1.2
    lim_x = np.max(np.concatenate((np.abs(ch_a_g), np.abs(ch_a_e))))*1.2
    lim_y = np.max(np.concatenate((np.abs(ch_b_g), np.abs(ch_b_e))))*1.2                             
    
#    lim = 5
    ax_scatter.set_xlim([-lim, lim])
    ax_scatter.set_ylim([-lim, lim])
    
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    
    plt.xlabel('ch1')
    plt.ylabel('ch2')
    
    plt.gca().set_aspect('equal')
    plt.show()
    
    
    plt.figure()
    plt.hist(ch_b_g, bins=bins_y, color = 'b', histtype='step', orientation='vertical')

    plt.hist(ch_b_e, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
    
    plt.xlabel('ch2 reading')
    plt.ylabel('Count')
    plt.legend(['g state', 'e state'],loc='upper left')    
    
    plt.show()
    
    plt.figure()
    plt.hist(ch_a_g, bins=bins_y, color = 'b', histtype='step', orientation='vertical')

    plt.hist(ch_a_e, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
    
    plt.xlabel('ch1 reading')
    plt.ylabel('Count')
    plt.legend(['g state', 'e state'],loc='upper left') 
    
    plt.show()
    
    plt.figure()
    binwidth = 1
    bins_y = np.arange(-lim, lim, binwidth)
    hist_g, bins = np.histogram(ch_b_g, bins=bins_y)
    hist_e, bins = np.histogram(ch_b_e, bins=bins_y)
    intcount_g = np.zeros(np.size(hist_g))
    intcount_e = np.zeros(np.size(hist_g))
    
    for k in np.arange(np.size(hist_g)):
        intcount_g[k] = np.sum(hist_g[0:k+1])
        intcount_e[k] = np.sum(hist_e[0:k+1])
        
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g), 'b-')
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_e/np.max(intcount_e), 'r-')
    plt.xlabel('ch1 read')
    plt.ylabel('Integrated count')
    plt.legend(['g state', 'e state'])  
    plt.show()
    
     
    plt.figure()
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e), 'b-')
    plt.xlabel('ch1 read')
    plt.ylabel('readout fidelity')
    fid = intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e)
    voltage = bins_y[0:-1]+binwidth/2
    index = np.argmax(fid)
    print('Cutoff voltage is ' + str(voltage[index]))
    print('Readout fidelity is ' + str(fid[index]))
    plt.plot([voltage[index], voltage[index]], plt.gca().get_ylim(), 'r')
    maxfid = np.max(fid)
    fid_g = intcount_g[index]/np.max(intcount_g)
    fid_e = 1-intcount_e[index]/np.max(intcount_e)
    plt.show()
    
    return maxfid, fid_g, fid_e


def compare_ge_readout_simultaneously(measure=1, num_steps=5, numrec=1000):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    tr.check_readout_multiplestate(measure=measure, num_steps=5, measTime = 2000)
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=5, num_records_per_pattern=numrec))
    
    
    rec_rot,_ = process_readout_raw(rec_all_raw_ave, ifplot=1)
    
    ch_a = rec_rot[0]
    ch_b = rec_rot[1]
    
    ch_a_reshape = np.reshape(ch_a,[np.int(np.size(ch_a)/(2*num_steps)),2*num_steps])
    ch_b_reshape = np.reshape(ch_b,[np.int(np.size(ch_a)/(2*num_steps)),2*num_steps])
    
    ch_a_g = ch_a_reshape[:,0:num_steps-2].flatten()
    ch_b_g = ch_b_reshape[:,0:num_steps-2].flatten()
    
    ch_a_e = ch_a_reshape[:,num_steps:2*num_steps-2].flatten()
    ch_b_e = ch_b_reshape[:,num_steps:2*num_steps-2].flatten()
    
    
    plt.figure()
    ax_scatter = plt.axes()
    ax_scatter.scatter(ch_a_g, ch_b_g, s=2, c='b', marker='.', alpha=0.2)
    ax_scatter.scatter(ch_a_e, ch_b_e, s=2, c='r', marker='.', alpha=0.2)
    
    lim = np.max(np.concatenate((np.abs(ch_a_g), np.abs(ch_a_e),
                                 np.abs(ch_b_g), np.abs(ch_b_e))))*1.2
    lim_x = np.max(np.concatenate((np.abs(ch_a_g), np.abs(ch_a_e))))*1.2
    lim_y = np.max(np.concatenate((np.abs(ch_b_g), np.abs(ch_b_e))))*1.2                             
    
#    lim = 5
    ax_scatter.set_xlim([-lim, lim])
    ax_scatter.set_ylim([-lim, lim])
    
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    
    plt.xlabel('ch1')
    plt.ylabel('ch2')
    
    plt.gca().set_aspect('equal')
    plt.show()
    
    
    plt.figure()
    plt.hist(ch_b_g, bins=bins_y, color = 'b', histtype='step', orientation='vertical')

    plt.hist(ch_b_e, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
    
    plt.xlabel('ch2 reading')
    plt.ylabel('Count')
    plt.legend(['g state', 'e state'],loc='upper left')    
    
    plt.show()
    
    plt.figure()
    plt.hist(ch_a_g, bins=bins_y, color = 'b', histtype='step', orientation='vertical')

    plt.hist(ch_a_e, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
    
    plt.xlabel('ch1 reading')
    plt.ylabel('Count')
    plt.legend(['g state', 'e state'],loc='upper left') 
    
    plt.show()
    
    plt.figure()
    binwidth = 1
    bins_y = np.arange(-lim, lim, binwidth)
    hist_g, bins = np.histogram(ch_b_g, bins=bins_y)
    hist_e, bins = np.histogram(ch_b_e, bins=bins_y)
    intcount_g = np.zeros(np.size(hist_g))
    intcount_e = np.zeros(np.size(hist_g))
    
    for k in np.arange(np.size(hist_g)):
        intcount_g[k] = np.sum(hist_g[0:k+1])
        intcount_e[k] = np.sum(hist_e[0:k+1])
        
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g), 'b-')
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_e/np.max(intcount_e), 'r-')
    plt.xlabel('ch1 read')
    plt.ylabel('Integrated count')
    plt.legend(['g state', 'e state'])  
    plt.show()
    
     
    plt.figure()
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e), 'b-')
    plt.xlabel('ch1 read')
    plt.ylabel('readout fidelity')
    fid = intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e)
    voltage = bins_y[0:-1]+binwidth/2
    index = np.argmax(fid)
    print('Cutoff voltage is ' + str(voltage[index]))
    print('Readout fidelity is ' + str(fid[index]))
    plt.plot([voltage[index], voltage[index]], plt.gca().get_ylim(), 'r')
    maxfid = np.max(fid)
    fid_g = intcount_g[index]/np.max(intcount_g)
    fid_e = 1-intcount_e[index]/np.max(intcount_e)
    plt.show()
    
    return maxfid, fid_g, fid_e

def compare_ge_readout_simultaneously_noSequence(num_steps=5, numrec=1000, 
                                                 ifplot=1,ifplotSNR=0,ifplotamp=0,ifreturn=0):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=num_steps*2, num_records_per_pattern=numrec))
    
    rec_rot,_ = process_readout_raw(rec_all_raw_ave, ifplot=0)
    
    
    if ifplot:
        
        shape1 = np.shape(rec_all_raw_ave)[1]
        shape2 = np.shape(rec_all_raw_ave)[2]
        cha_raw_reshape = np.reshape(rec_all_raw_ave[0], [np.int(shape1/(2*num_steps)),
                                     2*num_steps, shape2])
        chb_raw_reshape = np.reshape(rec_all_raw_ave[1], [np.int(shape1/(2*num_steps)),
                                     2*num_steps, shape2])
        plt.figure()
        plt.plot(np.mean(cha_raw_reshape[:,0:num_steps,:],axis=(0,1)),'b')
        plt.plot(np.mean(cha_raw_reshape[:,num_steps:2*num_steps,:],axis=(0,1)),'r')
        plt.xlabel('Time (ns)')
        plt.ylabel('Ch1 average read')
        plt.show()
        
        plt.figure()
        plt.plot(cha_raw_reshape[0,0,:],'b')
        plt.plot(cha_raw_reshape[0,num_steps,:],'r')
        plt.xlabel('Time (ns)')
        plt.ylabel('Ch1 single read')
        plt.show()
        
        plt.figure()
        plt.plot(np.mean(chb_raw_reshape[:,0:num_steps,:],axis=(0,1)),'b')
        plt.plot(np.mean(chb_raw_reshape[:,num_steps:2*num_steps,:],axis=(0,1)),'r')
        plt.xlabel('Time (ns)')
        plt.ylabel('Ch2 average read')
        plt.xlim([0,1000])
        plt.show()
    
    ch_a = rec_rot[0]
    ch_b = rec_rot[1]
    
    
    ch_a_reshape = np.reshape(ch_a,[np.int(np.size(ch_a)/(2*num_steps)),2*num_steps])
    ch_b_reshape = np.reshape(ch_b,[np.int(np.size(ch_a)/(2*num_steps)),2*num_steps])
    
    ch_a_g = ch_a_reshape[:,0:num_steps-0].flatten()
    ch_b_g = ch_b_reshape[:,0:num_steps-0].flatten()
    
    ch_a_e = ch_a_reshape[:,num_steps:2*num_steps-0].flatten()
    ch_b_e = ch_b_reshape[:,num_steps:2*num_steps-0].flatten()
    
    
    lim = np.max(np.concatenate((np.abs(ch_a_g), np.abs(ch_a_e),
                                 np.abs(ch_b_g), np.abs(ch_b_e))))*1.2
    lim_x = np.max(np.concatenate((np.abs(ch_a_g), np.abs(ch_a_e))))*1.2
    lim_y = np.max(np.concatenate((np.abs(ch_b_g), np.abs(ch_b_e))))*1.2   
    
    
    if ifplot:
#        ax = plt.figure()
#        plt.plot(rec_all_raw_ave[0][0,:],'r')
#        plt.plot(rec_all_raw_ave[1][0,:],'b')
        
        plt.figure()
        ax_scatter = plt.axes()
        ax_scatter.scatter(ch_a_g, ch_b_g, s=3, c='b', marker='.', alpha=0.2)
        ax_scatter.scatter(ch_a_e, ch_b_e, s=3, c='r', marker='.', alpha=0.2)
    
        ax_scatter.set_xlim([-lim, lim])
        ax_scatter.set_ylim([-lim, lim])
        
        
        ax_scatter.plot(np.mean(ch_a_g),np.mean(ch_b_g), 'bo')
        ax_scatter.plot(np.mean(ch_a_e),np.mean(ch_b_e), 'kx')
        
        
        vec = 1j * (np.mean(ch_b_e)-np.mean(ch_b_g)) +  (np.mean(ch_a_e)-np.mean(ch_a_g))
        angle_c = 180/np.pi * np.angle(vec)
        
        
        print('Correct angle is {}'.format(np.mod(angle+90-angle_c,360)))
        
        plt.xlabel('ch1')
        plt.ylabel('ch2')
        
        plt.gca().set_aspect('equal')
        plt.show()
    
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    
    
    
    
    hist_g, bins = np.histogram(ch_b_g, bins=bins_y)
    hist_e, bins = np.histogram(ch_b_e, bins=bins_y)
    binwidth = bins[1]-bins[0]
    
    
    intcount_g = np.zeros(np.size(hist_g))
    intcount_e = np.zeros(np.size(hist_g))
    
    for k in np.arange(np.size(hist_g)):
        intcount_g[k] = np.sum(hist_g[0:k+1])
        intcount_e[k] = np.sum(hist_e[0:k+1])
    
    
    fid = intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e)
    voltage = bins_y[0:-1]+binwidth/2
    index = np.argmax(fid)
    
    print('Readout fidelity is ' + str(fid[index]))
    print('Cutoff voltage is ' + str(voltage[index]))
    
    if ifplot:
        plt.figure()
        plt.hist(ch_b_g, bins=bins_y, color = 'b', histtype='step', orientation='vertical')
    
        plt.hist(ch_b_e, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
        
        plt.xlabel('ch2 reading')
        plt.ylabel('Count')
        plt.legend(['g state', 'e state'],loc='upper left')    
        
        plt.title('Fidelity is {}'.format(fid[index]))
        
        plt.show()
    
    
    
    
    xdata = voltage
    ydata_g = hist_g
    ydata_e = hist_e
    
    
    if ifplotSNR:
        index = np.argmax(ydata_g)
        startpoint = [np.max(ydata_g), xdata[index], np.ptp(xdata)/6]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        index = np.argmax(ydata_e)
        startpoint = [np.max(ydata_e), xdata[index], np.ptp(xdata)/6]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.figure()
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        
        print(startpoint)
        print(para_e)
        plt.title('SNR is {}'.format(SNR))
        plt.show()
    
    bins_amp = np.linspace(-lim/3, lim, binnum)
    hist_g_amp, bins = np.histogram(np.sqrt(ch_b_g**2 + ch_a_g**2), bins=bins_amp)
    hist_e_amp, bins = np.histogram(np.sqrt(ch_b_e**2 + ch_a_e**2), bins=bins_amp)
    
    binwidth = bins[1]-bins[0]
    xdata = bins_amp[0:-1]+binwidth/2
    voltage_amp = xdata    
        
    if ifplotamp:
        
        plt.figure()
#        plt.plot(xdata, hist_g_amp, 'b.-')
#        plt.plot(xdata, hist_e_amp, 'r.-')
        plt.xlabel('Amplitude')
        plt.ylabel('Count')
        
        ydata_g = hist_g_amp
        ydata_e = hist_e_amp
        
        startpoint = [np.max(ydata_g), 0, np.ptp(xdata)/4]
        
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        index = np.argmax(ydata_e)
        startpoint = [np.max(ydata_e), xdata[index], np.ptp(xdata)/4]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        
        plt.title('SNR is {}'.format(SNR))
        
        plt.show()
    
    rawdata = {'ch_b_g':ch_b_g,'ch_b_e':ch_b_e,'ch_a_g':ch_a_g,'ch_a_e':ch_a_e}
    
    if ifreturn == 1:
        return bins_y, voltage, hist_g, hist_e, fid[index]
    elif ifreturn == 2:
        return voltage, hist_g, hist_e, voltage_amp, hist_g_amp, hist_e_amp, rawdata

    
    
    

    
    
    
    

def compare_ge_readout_withNoise(noise=1):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    tr.check_readout_strongnoise(state = 0, measure = 1, noise=noise, num_steps = 5)
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=5, num_records_per_pattern=1000))
    
    
    rec_rot,_ = process_readout_raw(rec_all_raw_ave, ifplot=1)
    
    ch_a_g = rec_rot[0]
    ch_b_g = rec_rot[1]
    
    tr.check_readout_strongnoise(state = 1, measure = 1, noise=noise, num_steps = 5)
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=5, num_records_per_pattern=1000))
    
    
    rec_rot,_ = process_readout_raw(rec_all_raw_ave)
    
    ch_a_e = rec_rot[0]
    ch_b_e = rec_rot[1]
    
    
    plt.figure()
    ax_scatter = plt.axes()
    ax_scatter.scatter(ch_a_g, ch_b_g, s=0.5, c='b', marker='.', alpha=0.2)
    ax_scatter.scatter(ch_a_e, ch_b_e, s=0.5, c='r', marker='.', alpha=0.2)
    
    lim = np.max(np.concatenate((np.abs(ch_a_g), np.abs(ch_a_e),
                                 np.abs(ch_b_g), np.abs(ch_b_e))))*1.2
    lim_x = np.max(np.concatenate((np.abs(ch_a_g), np.abs(ch_a_e))))*1.2
    lim_y = np.max(np.concatenate((np.abs(ch_b_g), np.abs(ch_b_e))))*1.2                             
    
#    lim = 5
    ax_scatter.set_xlim([-lim, lim])
    ax_scatter.set_ylim([-lim, lim])
    
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    
    plt.xlabel('ch1')
    plt.ylabel('ch2')
    
    plt.gca().set_aspect('equal')
    plt.show()
    
    
    plt.figure()
    plt.hist(ch_b_g, bins=bins_y, color = 'b', histtype='step', orientation='vertical')

    plt.hist(ch_b_e, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
    
    plt.xlabel('ch2 reading')
    plt.ylabel('Count')
    plt.legend(['g state', 'e state'],loc='upper left')    
    
    plt.show()
    
    plt.figure()
    plt.hist(ch_a_g, bins=bins_y, color = 'b', histtype='step', orientation='vertical')

    plt.hist(ch_a_e, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
    
    plt.xlabel('ch1 reading')
    plt.ylabel('Count')
    plt.legend(['g state', 'e state'],loc='upper left') 
    
    plt.show()
    
    plt.figure()
    binwidth = 1
    bins_y = np.arange(-lim, lim, binwidth)
    hist_g, bins = np.histogram(ch_b_g, bins=bins_y)
    hist_e, bins = np.histogram(ch_b_e, bins=bins_y)
    intcount_g = np.zeros(np.size(hist_g))
    intcount_e = np.zeros(np.size(hist_g))
    
    for k in np.arange(np.size(hist_g)):
        intcount_g[k] = np.sum(hist_g[0:k+1])
        intcount_e[k] = np.sum(hist_e[0:k+1])
        
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g), 'b-')
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_e/np.max(intcount_e), 'r-')
    plt.xlabel('ch1 read')
    plt.ylabel('Integrated count')
    plt.legend(['g state', 'e state'])  
    plt.show()
    
     
    plt.figure()
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e), 'b-')
    plt.xlabel('ch1 read')
    plt.ylabel('readout fidelity')
    fid = intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e)
    voltage = bins_y[0:-1]+binwidth/2
    index = np.argmax(fid)
    print('Cutoff voltage is ' + str(voltage[index]))
    print('Readout fidelity is ' + str(fid[index]))
    plt.plot([voltage[index], voltage[index]], plt.gca().get_ylim(), 'r')
    maxfid = np.max(fid)
    fid_g = intcount_g[index]/np.max(intcount_g)
    fid_e = 1-intcount_e[index]/np.max(intcount_e)
    plt.show()
    
    return maxfid, fid_g, fid_e


def compare_ge_readout_FT(numrec=1000, ifplotSNR=0):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    tr.check_readout(state = 0, measure = 0, num_steps = 5)
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=5, num_records_per_pattern=numrec))
    
    
    sig_subBg_gpeak_gstate, sig_subBg_epeak_gstate, spec_gstate, freq = process_thermalreadout_raw(rec_all_raw_ave, ifplot=1)
    sig_gstate = sig_subBg_epeak_gstate - sig_subBg_gpeak_gstate
    
    
    tr.check_readout(state = 1, measure = 0, num_steps = 5)
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=5, num_records_per_pattern=numrec))
    
    
    sig_subBg_gpeak_estate, sig_subBg_epeak_estate, spec_estate, freq = process_thermalreadout_raw(rec_all_raw_ave, ifplot=1)
    sig_estate = sig_subBg_epeak_estate - sig_subBg_gpeak_estate
    
    print(np.shape(sig_estate))
    print(np.shape(sig_gstate))
#    lim = np.max(np.concatenate((sig_estate,sig_gstate)))*1.2
    lim = np.max(sig_estate)*1.2   
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    
    
    
    plt.figure()
    plt.plot(freq, np.mean(spec_gstate[:,:],axis=0), 'b-')
    plt.plot(freq, np.mean(spec_estate[:,:],axis=0), 'r-')
    plt.legend(['g state','e state'])
    plt.xlabel('Freq (GHz)')
    plt.ylabel('FT signal')
    plt.xlim([0,0.02])
    autoscale_y(plt.gca())
    plt.show()
    

    
    plt.figure()
    bins_y = np.linspace(-lim, lim, binnum)
    binwidth = bins_y[1] - bins_y[0]
    hist_g, bins = np.histogram(sig_gstate, bins=bins_y)
    hist_e, bins = np.histogram(sig_estate, bins=bins_y)
    intcount_g = np.zeros(np.size(hist_g))
    intcount_e = np.zeros(np.size(hist_g))
    
    for k in np.arange(np.size(hist_g)):
        intcount_g[k] = np.sum(hist_g[0:k+1])
        intcount_e[k] = np.sum(hist_e[0:k+1])
        
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g), 'b-')
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_e/np.max(intcount_e), 'r-')
    plt.xlabel('ch1 read')
    plt.ylabel('Integrated count')
    plt.legend(['g state', 'e state'])  
    plt.show()
    
     
    plt.figure()
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e), 'b-')
    plt.xlabel('ch1 read')
    plt.ylabel('readout fidelity')
    fid = intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e)
    voltage = bins_y[0:-1]+binwidth/2
    index = np.argmax(fid)
    print('Cutoff voltage is ' + str(voltage[index]))
    print('Readout fidelity is ' + str(fid[index]))
    plt.plot([voltage[index], voltage[index]], plt.gca().get_ylim(), 'r')
    maxfid = np.max(fid)
    fid_g = intcount_g[index]/np.max(intcount_g)
    fid_e = 1-intcount_e[index]/np.max(intcount_e)
    plt.show()
    
    
    
    plt.figure()
    bins_y = np.linspace(-lim, lim, binnum)
    plt.hist(sig_gstate, bins=bins_y, color = 'b', histtype='step', orientation='vertical')

    plt.hist(sig_estate, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
    
    plt.xlabel('FT signal')
    plt.ylabel('Count')
    plt.legend(['g state', 'e state'],loc='upper left')    
    plt.title('Readout fidelity is ' + str(fid[index]))
    plt.show()
    
    
    
    xdata = voltage
    ydata_g = hist_g
    ydata_e = hist_e
    
    if ifplotSNR:
        startpoint = [np.max(ydata_g), 0, np.ptp(xdata)/4]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        startpoint = [np.max(ydata_e), 0, np.ptp(xdata)/4]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.figure()
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        
        
        plt.title('SNR is {}'.format(SNR))
        plt.show()
#    return maxfid, fid_g, fid_e

def compare_ge_readout_FT_simultaneously_nosequence(num_steps=5, numrec=2000, ifplot=1, ifplotSNR=1,ifsave=0):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    fg = 5.6185 - 5.604
    fe = 5.606 - 5.604
    
    gpeak = fg  # 10 dB
    epeak = fe    
    
    linewidth = 0.0005 
        
        
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=num_steps, num_records_per_pattern=numrec))
    
    
    (sig_subBg_gpeak, sig_subBg_epeak,
     spec, specbg, freq) = process_thermalreadout_raw(rec_all_raw_ave,
                       gpeak=gpeak, epeak=epeak, linewidth=linewidth, ifplot=0)
    sig = sig_subBg_epeak - sig_subBg_gpeak
#    sig = sig_subBg_epeak 
    
    
    sig_reshape = np.reshape(sig,[np.int(np.size(sig)/(2*num_steps)),2*num_steps])
    
#    sig_gstate = sig_reshape[:,0:num_steps-0].flatten()
#    sig_estate = sig_reshape[:,num_steps:2*num_steps-0].flatten()
    
    sig_gstate = sig_reshape[:,0:num_steps-1].flatten()
    sig_estate = sig_reshape[:,num_steps:2*num_steps-1].flatten()
    
    spec_reshape = np.reshape(spec,[np.int(np.shape(spec)[0]/(2*num_steps)),2*num_steps,np.shape(spec)[1]])
    spec_gstate = spec_reshape[:,0:num_steps-1,:]
    spec_estate = spec_reshape[:,num_steps:2*num_steps-1,:]
    
    specbg_reshape = np.reshape(specbg,[np.int(np.shape(spec)[0]/(2*num_steps)),2*num_steps,np.shape(spec)[1]])
    specbg_gstate = specbg_reshape[:,0:num_steps-1,:]
    specbg_estate = specbg_reshape[:,num_steps:2*num_steps-1,:]
    
    mean_spec_gstate_subbg = np.mean(spec_gstate, axis = (0,1)) - np.mean(specbg_gstate, axis = (0,1))
    mean_spec_estate_subbg = np.mean(spec_estate, axis = (0,1)) - np.mean(specbg_estate, axis = (0,1))
    
    mean_spec_gstate = np.mean(spec_gstate, axis = (0,1)) 
    mean_spec_estate = np.mean(spec_estate, axis = (0,1)) 
    
    mean_spec_gstate_bg = np.mean(specbg_gstate, axis = (0,1)) 
    mean_spec_estate_bg = np.mean(specbg_estate, axis = (0,1)) 

    lim = np.max(sig_estate)*1.2   
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    
    
    
#    fg = 5.6185 - 5.6205
#    fe = 5.606 - 5.6205
    
    if ifplot:
        fig, ax = plt.subplots()
        plt.plot(1e3*freq, mean_spec_gstate_subbg, 'b-')
        plt.plot(1e3*freq, mean_spec_estate_subbg, 'r-')
        plt.legend(['g state','e state'])
        plt.xlabel('Freq (MHz)')
        plt.ylabel('FT signal')
        plt.xlim([-20,20])
        plt.ylim([-1000, 3000])
#        autoscale_y(plt.gca())
        
        
        plt.plot(1e3*np.array([fg, fg]), plt.gca().get_ylim(), 'b--')
        plt.plot(1e3*np.array([-fg, -fg]), plt.gca().get_ylim(), 'b--')
        
        
        plt.plot(1e3*np.array([fe, fe]), plt.gca().get_ylim(), 'r--')
        plt.plot(1e3*np.array([-fe, -fe]), plt.gca().get_ylim(), 'r--')
        
        
            
        
        plt.show()
        
    if ifplot:
        fig, ax = plt.subplots()
        plt.plot(1e3*freq, mean_spec_gstate, 'b-')
        plt.plot(1e3*freq, mean_spec_estate, 'r-')
#        plt.plot(1e3*freq, mean_spec_gstate_bg, 'k--')
        
        
        plt.legend(['g state','e state'])
        
#        plt.plot(freq, spec_gstate[0,0,:], 'b--')
#        plt.plot(freq, spec_estate[0,0,:], 'r--')
        
        plt.xlabel('Freq (MHz)')
        plt.ylabel('FT signal')
        plt.xlim([-20,20])
        plt.ylim([-50, 10000])
        
        
        plt.plot([1e3*fg, 1e3*fg], plt.gca().get_ylim(), 'b--')
        plt.plot([-1e3*fg, -1e3*fg], plt.gca().get_ylim(), 'b--')
        
        
        plt.plot([1e3*fe, 1e3*fe], plt.gca().get_ylim(), 'r--')
        plt.plot([-1e3*fe, -1e3*fe], plt.gca().get_ylim(), 'r--')
#        autoscale_y(plt.gca())
        
        
        gpeak = gpeak * 1e3  # 10 dB
        epeak = epeak * 1e3    
        
        linewidth = linewidth * 1e3
        ylimp = plt.gca().get_ylim()
        ax.add_patch(Rectangle((epeak-linewidth,ylimp[0]),2*linewidth,np.ptp(ylimp)
                               ,alpha=0.1,color='r',edgecolor=None))
        ax.add_patch(Rectangle((-epeak-linewidth,ylimp[0]),2*linewidth,np.ptp(ylimp)
                               ,alpha=0.1,color='r',edgecolor=None))
        ax.add_patch(Rectangle((gpeak-linewidth,ylimp[0]),2*linewidth,np.ptp(ylimp)
                               ,alpha=0.1,color='b',edgecolor=None))
        ax.add_patch(Rectangle((-gpeak-linewidth,ylimp[0]),2*linewidth,np.ptp(ylimp)
                               ,alpha=0.1,color='b',edgecolor=None))
        
        
        plt.show()
        
        
        
    

    
    
    bins_y = np.linspace(-lim, lim, binnum)
    binwidth = bins_y[1] - bins_y[0]
    hist_g, bins = np.histogram(sig_gstate, bins=bins_y)
    hist_e, bins = np.histogram(sig_estate, bins=bins_y)
    intcount_g = np.zeros(np.size(hist_g))
    intcount_e = np.zeros(np.size(hist_g))
    
    for k in np.arange(np.size(hist_g)):
        intcount_g[k] = np.sum(hist_g[0:k+1])
        intcount_e[k] = np.sum(hist_e[0:k+1])
    
    if ifplot:
        plt.figure()    
        plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g), 'b-')
        plt.plot(bins_y[0:-1]+binwidth/2, intcount_e/np.max(intcount_e), 'r-')
        plt.xlabel('ch1 read')
        plt.ylabel('Integrated count')
        plt.legend(['g state', 'e state'])  
        plt.show()
    
    
    fid = intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e)
    voltage = bins_y[0:-1]+binwidth/2
    index = np.argmax(fid)
    if ifplot: 
        plt.figure()
        plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e), 'b-')
        plt.xlabel('ch1 read')
        plt.ylabel('readout fidelity')
        
        print('Cutoff voltage is ' + str(voltage[index]))
        print('Readout fidelity is ' + str(fid[index]))
        plt.plot([voltage[index], voltage[index]], plt.gca().get_ylim(), 'r')
        maxfid = np.max(fid)
        fid_g = intcount_g[index]/np.max(intcount_g)
        fid_e = 1-intcount_e[index]/np.max(intcount_e)
        plt.show()
    
    
    if ifplot:
        plt.figure()
        bins_y = np.linspace(-lim, lim, binnum)
        plt.hist(sig_gstate, bins=bins_y, color = 'b', histtype='step', orientation='vertical')
    
        plt.hist(sig_estate, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
        
        plt.xlabel('FT signal')
        plt.ylabel('Count')
        plt.legend(['g state', 'e state'],loc='upper left')    
        plt.title('Readout fidelity is ' + str(fid[index]))
        plt.show()
    
    
    
    xdata = voltage
    ydata_g = hist_g
    ydata_e = hist_e
    
    if ifplotSNR:
        startpoint = [np.max(ydata_g), 0, np.ptp(xdata)/4]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        startpoint = [np.max(ydata_e), 0, np.ptp(xdata)/4]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.figure()
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        plt.xlabel('FT signal')
        plt.ylabel('Count')
        
        plt.title('SNR is {}'.format(SNR))
        plt.show()
        
    if ifsave:
        data = {'freq':freq, 'mean_spec_gstate_subbg':mean_spec_gstate_subbg,
                'mean_spec_estate_subbg':mean_spec_gstate_subbg,
                'mean_spec_gstate':mean_spec_gstate,
                'mean_spec_estate':mean_spec_estate,
                'mean_spec_gstate_bg': mean_spec_gstate_bg,
                'mean_spec_estate_bg': mean_spec_estate_bg,
                'sig_gstate': sig_gstate,
                'sig_estate': sig_estate}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\02_coherent_readout',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)
    
        
    return bins_y, voltage, hist_g, hist_e, fid[index]      


def compare_ge_readout_FT_simultaneously_nosequence_weight(num_steps=5, numrec=2000, ifplot=1, ifplotSNR=1,ifsave=0):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    fg = 5.6185 - 5.604
    fe = 5.606 - 5.604
    
    gpeak = fg  # 10 dB
    epeak = fe    
    
    linewidth = 0.002 
        
        
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=num_steps, num_records_per_pattern=numrec))
    
    
    (spec, specbg, freq) = process_thermalreadout_raw_weight(rec_all_raw_ave, ifplot=0)
   
    
    
    
    spec_reshape = np.reshape(spec,[np.int(np.shape(spec)[0]/(2*num_steps)),2*num_steps,np.shape(spec)[1]])
    spec_gstate = spec_reshape[:,0:num_steps-0,:]
    spec_estate = spec_reshape[:,num_steps:2*num_steps-0,:]
    
    specbg_reshape = np.reshape(specbg,[np.int(np.shape(spec)[0]/(2*num_steps)),2*num_steps,np.shape(spec)[1]])
    specbg_gstate = specbg_reshape[:,0:num_steps-0,:]
    specbg_estate = specbg_reshape[:,num_steps:2*num_steps-0,:]
    
    mean_spec_gstate_subbg = np.mean(spec_gstate, axis = (0,1)) - np.mean(specbg_gstate, axis = (0,1))
    mean_spec_estate_subbg = np.mean(spec_estate, axis = (0,1)) - np.mean(specbg_estate, axis = (0,1))
    
    mean_spec_gstate = np.mean(spec_gstate, axis = (0,1)) 
    mean_spec_estate = np.mean(spec_estate, axis = (0,1)) 
    
    mean_spec_gstate_bg = np.mean(specbg_gstate, axis = (0,1)) 
    mean_spec_estate_bg = np.mean(specbg_estate, axis = (0,1)) 
    
    
    weight_spectrum = mean_spec_estate - mean_spec_gstate
    windex =  np.where(np.abs(weight_spectrum)<20)[0]
    weight_spectrum[windex] = 0
    
    findex =  np.where(np.logical_and(freq>-20,freq<20))[0]   
    weight = np.reshape(weight_spectrum,(1,np.size(weight_spectrum)))
    print(np.shape(weight))
    print(np.shape(spec))
    weight = np.tile(weight, (np.shape(spec)[0],1))
    print(np.shape(weight))
    sig = np.mean(spec[:,findex]*weight[:,findex],axis=1)
    
    sig_reshape = np.reshape(sig,[np.int(np.size(sig)/(2*num_steps)),2*num_steps])
    
#    sig_gstate = sig_reshape[:,0:num_steps-0].flatten()
#    sig_estate = sig_reshape[:,num_steps:2*num_steps-0].flatten()
    
    sig_gstate = sig_reshape[:,0:num_steps-0].flatten()
    sig_estate = sig_reshape[:,num_steps:2*num_steps-0].flatten()
    
    lim = np.max(sig_estate)*1.2   
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    
    
    
#    fg = 5.6185 - 5.6205
#    fe = 5.606 - 5.6205
    
    if ifplot:
        fig, ax = plt.subplots()
#        plt.plot(1e3*freq, mean_spec_gstate_subbg, 'b-')
#        plt.plot(1e3*freq, mean_spec_estate_subbg, 'r-')
        plt.plot(1e3*freq, weight_spectrum, 'k-')
#        plt.legend(['g state','e state'])
        plt.xlabel('Freq (MHz)')
        plt.ylabel('FT signal')
        plt.xlim([-20,20])
        plt.ylim([-1000, 3000])
#        autoscale_y(plt.gca())
        
        
        plt.plot(1e3*np.array([fg, fg]), plt.gca().get_ylim(), 'b--')
        plt.plot(1e3*np.array([-fg, -fg]), plt.gca().get_ylim(), 'b--')
        
        
        plt.plot(1e3*np.array([fe, fe]), plt.gca().get_ylim(), 'r--')
        plt.plot(1e3*np.array([-fe, -fe]), plt.gca().get_ylim(), 'r--')
        
        
            
        
        plt.show()
        
    if ifplot:
        fig, ax = plt.subplots()
        plt.plot(1e3*freq, mean_spec_gstate, 'b-')
        plt.plot(1e3*freq, mean_spec_estate, 'r-')
#        plt.plot(1e3*freq, mean_spec_gstate_bg, 'k--')
        
        
        plt.legend(['g state','e state'])
        
#        plt.plot(freq, spec_gstate[0,0,:], 'b--')
#        plt.plot(freq, spec_estate[0,0,:], 'r--')
        
        plt.xlabel('Freq (MHz)')
        plt.ylabel('FT signal')
        plt.xlim([-20,20])
        plt.ylim([-50, 10000])
        
        
        plt.plot([1e3*fg, 1e3*fg], plt.gca().get_ylim(), 'b--')
        plt.plot([-1e3*fg, -1e3*fg], plt.gca().get_ylim(), 'b--')
        
        
        plt.plot([1e3*fe, 1e3*fe], plt.gca().get_ylim(), 'r--')
        plt.plot([-1e3*fe, -1e3*fe], plt.gca().get_ylim(), 'r--')
#        autoscale_y(plt.gca())
        
        
        gpeak = gpeak * 1e3  # 10 dB
        epeak = epeak * 1e3    
        
        linewidth = linewidth * 1e3
        ylimp = plt.gca().get_ylim()
#        ax.add_patch(Rectangle((epeak-linewidth,ylimp[0]),2*linewidth,np.ptp(ylimp)
#                               ,alpha=0.1,color='r',edgecolor=None))
#        ax.add_patch(Rectangle((-epeak-linewidth,ylimp[0]),2*linewidth,np.ptp(ylimp)
#                               ,alpha=0.1,color='r',edgecolor=None))
#        ax.add_patch(Rectangle((gpeak-linewidth,ylimp[0]),2*linewidth,np.ptp(ylimp)
#                               ,alpha=0.1,color='b',edgecolor=None))
#        ax.add_patch(Rectangle((-gpeak-linewidth,ylimp[0]),2*linewidth,np.ptp(ylimp)
#                               ,alpha=0.1,color='b',edgecolor=None))
        
        
        plt.show()
        
        
        
    

    
    
    bins_y = np.linspace(-lim, lim, binnum)
    binwidth = bins_y[1] - bins_y[0]
    hist_g, bins = np.histogram(sig_gstate, bins=bins_y)
    hist_e, bins = np.histogram(sig_estate, bins=bins_y)
    intcount_g = np.zeros(np.size(hist_g))
    intcount_e = np.zeros(np.size(hist_g))
    
    for k in np.arange(np.size(hist_g)):
        intcount_g[k] = np.sum(hist_g[0:k+1])
        intcount_e[k] = np.sum(hist_e[0:k+1])
    
    if ifplot:
        plt.figure()    
        plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g), 'b-')
        plt.plot(bins_y[0:-1]+binwidth/2, intcount_e/np.max(intcount_e), 'r-')
        plt.xlabel('ch1 read')
        plt.ylabel('Integrated count')
        plt.legend(['g state', 'e state'])  
        plt.show()
    
    
    fid = intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e)
    voltage = bins_y[0:-1]+binwidth/2
    index = np.argmax(fid)
    if ifplot: 
        plt.figure()
        plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e), 'b-')
        plt.xlabel('ch1 read')
        plt.ylabel('readout fidelity')
        
        print('Cutoff voltage is ' + str(voltage[index]))
        print('Readout fidelity is ' + str(fid[index]))
        plt.plot([voltage[index], voltage[index]], plt.gca().get_ylim(), 'r')
        maxfid = np.max(fid)
        fid_g = intcount_g[index]/np.max(intcount_g)
        fid_e = 1-intcount_e[index]/np.max(intcount_e)
        plt.show()
    
    
    if ifplot:
        plt.figure()
        bins_y = np.linspace(-lim, lim, binnum)
        plt.hist(sig_gstate, bins=bins_y, color = 'b', histtype='step', orientation='vertical')
    
        plt.hist(sig_estate, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
        
        plt.xlabel('FT signal')
        plt.ylabel('Count')
        plt.legend(['g state', 'e state'],loc='upper left')    
        plt.title('Readout fidelity is ' + str(fid[index]))
        plt.show()
    
    
    
    xdata = voltage
    ydata_g = hist_g
    ydata_e = hist_e
    
    if ifplotSNR:
        startpoint = [np.max(ydata_g), 0, np.ptp(xdata)/4]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        startpoint = [np.max(ydata_e), 0, np.ptp(xdata)/4]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.figure()
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        plt.xlabel('FT signal')
        plt.ylabel('Count')
        
        plt.title('SNR is {}'.format(SNR))
        plt.show()
        
    if ifsave:
        data = {'freq':freq, 'mean_spec_gstate_subbg':mean_spec_gstate_subbg,
                'mean_spec_estate_subbg':mean_spec_gstate_subbg,
                'mean_spec_gstate':mean_spec_gstate,
                'mean_spec_estate':mean_spec_estate,
                'mean_spec_gstate_bg': mean_spec_gstate_bg,
                'mean_spec_estate_bg': mean_spec_estate_bg,
                'sig_gstate': sig_gstate,
                'sig_estate': sig_estate}
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\02_coherent_readout',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)
    
        
    return bins_y, voltage, hist_g, hist_e, fid[index] 

def compare_ge_readout_noFourierTransform_simultaneously_nosequence(num_steps=5, numrec=2000, ifplot=1, ifplotSNR=1):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=num_steps, num_records_per_pattern=numrec))
    
    
    sig = process_thermalreadout_raw_noFourierTransform(rec_all_raw_ave, ifplot=ifplot)
        
    sig_reshape = np.reshape(sig,[np.int(np.size(sig)/(2*num_steps)),2*num_steps])
    
    sig_gstate = sig_reshape[:,0:num_steps-0].flatten()
    sig_estate = sig_reshape[:,num_steps:2*num_steps-0].flatten()
        
#    mean_spec_gstate = np.mean(spec_gstate, axis = (0,1)) - np.mean(specbg_gstate, axis = (0,1))
#    mean_spec_estate = np.mean(spec_estate, axis = (0,1)) - np.mean(specbg_estate, axis = (0,1))
    
#    mean_spec_gstate = np.mean(spec_gstate, axis = (0,1)) 
#    mean_spec_estate = np.mean(spec_estate, axis = (0,1)) 

    lim = np.max(sig_estate)*1.2   
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    
    
    
#    if ifplot:
#        plt.figure()
#        plt.plot(freq, mean_spec_gstate, 'b-')
#        plt.plot(freq, mean_spec_estate, 'r-')
#        plt.legend(['g state','e state'])
#        plt.xlabel('Freq (GHz)')
#        plt.ylabel('FT signal')
#        plt.xlim([-0.02,0.02])
#        plt.ylim([-50, 10000])
##        autoscale_y(plt.gca())
#        plt.show()
    

    
    
    bins_y = np.linspace(-lim, lim, binnum)
    binwidth = bins_y[1] - bins_y[0]
    hist_g, bins = np.histogram(sig_gstate, bins=bins_y)
    hist_e, bins = np.histogram(sig_estate, bins=bins_y)
    intcount_g = np.zeros(np.size(hist_g))
    intcount_e = np.zeros(np.size(hist_g))
    
    for k in np.arange(np.size(hist_g)):
        intcount_g[k] = np.sum(hist_g[0:k+1])
        intcount_e[k] = np.sum(hist_e[0:k+1])
    
    if ifplot:
        plt.figure()    
        plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g), 'b-')
        plt.plot(bins_y[0:-1]+binwidth/2, intcount_e/np.max(intcount_e), 'r-')
        plt.xlabel('ch1 read')
        plt.ylabel('Integrated count')
        plt.legend(['g state', 'e state'])  
        plt.show()
    
    
    fid = intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e)
    voltage = bins_y[0:-1]+binwidth/2
    index = np.argmax(fid)
    if ifplot: 
        plt.figure()
        plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e), 'b-')
        plt.xlabel('ch1 read')
        plt.ylabel('readout fidelity')
        
        print('Cutoff voltage is ' + str(voltage[index]))
        print('Readout fidelity is ' + str(fid[index]))
        plt.plot([voltage[index], voltage[index]], plt.gca().get_ylim(), 'r')
        maxfid = np.max(fid)
        fid_g = intcount_g[index]/np.max(intcount_g)
        fid_e = 1-intcount_e[index]/np.max(intcount_e)
        plt.show()
    
    
    if ifplot:
        plt.figure()
        bins_y = np.linspace(-lim, lim, binnum)
        plt.hist(sig_gstate, bins=bins_y, color = 'b', histtype='step', orientation='vertical')
    
        plt.hist(sig_estate, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
        
        plt.xlabel('FT signal')
        plt.ylabel('Count')
        plt.legend(['g state', 'e state'],loc='upper left')    
        plt.title('Readout fidelity is ' + str(fid[index]))
        plt.show()
    
    
    
    xdata = voltage
    ydata_g = hist_g
    ydata_e = hist_e
    
    if ifplotSNR:
        startpoint = [np.max(ydata_g), 0, np.ptp(xdata)/4]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        startpoint = [np.max(ydata_e), 0, np.ptp(xdata)/4]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.figure()
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        
        
        plt.title('SNR is {}'.format(SNR))
        plt.show()
        
    return bins_y, voltage, hist_g, hist_e, fid[index]   
        
def compare_ge_readout_FT_notFFT(numrec=1000, ifplotSNR=0):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    tr.check_readout(state = 0, measure = 0, num_steps = 5)
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=5, num_records_per_pattern=numrec))
    
    
    sig_gstate = process_thermalreadout_raw_notFFT(rec_all_raw_ave, ifplot=0)
    
    
    tr.check_readout(state = 1, measure = 0, num_steps = 5)
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=5, num_records_per_pattern=numrec))
    
    
    sig_estate = process_thermalreadout_raw_notFFT(rec_all_raw_ave, ifplot=0)
    
    print(np.shape(sig_estate))
    print(np.shape(sig_gstate))
#    lim = np.max(np.concatenate((sig_estate,sig_gstate)))*1.2
    lim = np.max(sig_estate)*1.2   
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    

    
    plt.figure()
    bins_y = np.linspace(-lim, lim, binnum)
    binwidth = bins_y[1] - bins_y[0]
    hist_g, bins = np.histogram(sig_gstate, bins=bins_y)
    hist_e, bins = np.histogram(sig_estate, bins=bins_y)
    intcount_g = np.zeros(np.size(hist_g))
    intcount_e = np.zeros(np.size(hist_g))
    
    for k in np.arange(np.size(hist_g)):
        intcount_g[k] = np.sum(hist_g[0:k+1])
        intcount_e[k] = np.sum(hist_e[0:k+1])
        
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g), 'b-')
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_e/np.max(intcount_e), 'r-')
    plt.xlabel('ch1 read')
    plt.ylabel('Integrated count')
    plt.legend(['g state', 'e state'])  
    plt.show()
    
     
    plt.figure()
    plt.plot(bins_y[0:-1]+binwidth/2, intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e), 'b-')
    plt.xlabel('ch1 read')
    plt.ylabel('readout fidelity')
    fid = intcount_g/np.max(intcount_g)-intcount_e/np.max(intcount_e)
    voltage = bins_y[0:-1]+binwidth/2
    index = np.argmax(fid)
    print('Cutoff voltage is ' + str(voltage[index]))
    print('Readout fidelity is ' + str(fid[index]))
    plt.plot([voltage[index], voltage[index]], plt.gca().get_ylim(), 'r')
    maxfid = np.max(fid)
    fid_g = intcount_g[index]/np.max(intcount_g)
    fid_e = 1-intcount_e[index]/np.max(intcount_e)
    plt.show()
    
    
    
    plt.figure()
    bins_y = np.linspace(-lim, lim, binnum)
    plt.hist(sig_gstate, bins=bins_y, color = 'b', histtype='step', orientation='vertical')

    plt.hist(sig_estate, bins=bins_y, color = 'r', histtype='step', orientation='vertical')
    
    plt.xlabel('FT signal')
    plt.ylabel('Count')
    plt.legend(['g state', 'e state'],loc='upper left')    
    plt.title('Readout fidelity is ' + str(fid[index]))
    plt.show()
    
    
    
    xdata = voltage
    ydata_g = hist_g
    ydata_e = hist_e
    
    if ifplotSNR:
        startpoint = [np.max(ydata_g), 0, np.ptp(xdata)/4]
        fitfunc = Gaussian;
        para_g, pcov_g = opt.curve_fit(fitfunc, xdata, ydata_g, startpoint)
        perr_g = np.sqrt(np.diag(pcov_g))
        
        startpoint = [np.max(ydata_e), 0, np.ptp(xdata)/4]
        para_e, pcov_e = opt.curve_fit(fitfunc, xdata, ydata_e, startpoint)
        perr_e = np.sqrt(np.diag(pcov_e))
    
        SNR = np.abs((para_e[1]-para_g[1]))/((np.abs(para_g[2])+np.abs(para_e[2]))/2)
        plt.figure()
        plt.plot(xdata, ydata_g, 'b.')
        plt.plot(xdata, ydata_e, 'r.')
        plt.plot(xdata, fitfunc(xdata,*para_g), 'b-')
        plt.plot(xdata, fitfunc(xdata,*para_e), 'r-')
        
        
        plt.title('SNR is {}'.format(SNR))
        plt.show()
#    return maxfid, fid_g, fid_e

    
def readout_singleRun():
   
    
    
    
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=10, num_records_per_pattern=500))
    
    rec_rot,prob_e = process_readout_raw(rec_all_raw_ave)
    
    
    ch_a_g = rec_rot[0]
    ch_b_g = rec_rot[1]
    
    

    plt.figure()
    plt.plot(np.mean(rec_all_raw_ave[0][:,:],axis=0),'r')
    plt.plot(np.mean(rec_all_raw_ave[1][:,:],axis=0),'b')
    plt.show()
    
    plt.figure()
    plt.plot(rec_all_raw_ave[0][1000,:],'r')
    plt.plot(rec_all_raw_ave[1][1000,:],'b')
    plt.show()
    
    plt.figure()
    ax_scatter = plt.axes()
    ax_scatter.scatter(ch_a_g, ch_b_g, c='r',s=0.5, marker='.', alpha=0.2)
    
    limx = np.max(np.abs(ch_a_g))*1.2
    limy = np.max(np.abs(ch_b_g))*1.2
    lim = np.max([limx, limy])
    ax_scatter.set_xlim([-lim, lim])
    ax_scatter.set_ylim([-lim, lim])
    plt.gca().set_aspect('equal')
    plt.show()
    
    
    
    binnum = 50
    bins_x = np.linspace(-limx, limx, binnum)
    bins_y = np.linspace(-limy, limy, binnum)
    
    counts_x = np.histogram(ch_a_g)
    counts_y = np.histogram(ch_b_g)
    
#    counts_x = np.histogram(ch_a, bins_x)
#    counts_y = np.histogram(ch_b, bins_y)
#    plt.hist(bins[:-1], bins, weights=counts)

#    ax_histx.hist(ch_a, bins=bins_x, weights=counts_x)
#    ax_histy.hist(ch_b, bins=bins_y, weights=counts_y, orientation='horizontal')
    
#    ax_histx.hist(ch_a, bins=bins_x, weights=counts_x)
#    ax_histy.hist(bins_x[:-1], bins=bins_x, weights=counts_y, orientation='horizontal')
    
#    plt.hist(ch_a, bins=bins_x, histtype='step', orientation='vertical')
    plt.figure()
    plt.hist(ch_b_g, bins=bins_y, histtype='step', orientation='vertical')
    plt.show()
    
    print('Probability in e is ' + str(prob_e))


def readout_FT_singlerun(num_steps=5, numrec=2000, ifplot=1, ifsave = 0):
#    startindex = 500
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    
    daq_params, rec_all_raw_ave = (
            dp.run_daq_rawdata(num_patterns=num_steps, num_records_per_pattern=numrec))
    
    
    sig_subBg_gpeak, sig_subBg_epeak, spec, specbg, freq = process_thermalreadout_raw(rec_all_raw_ave, ifplot=0)
    sig = sig_subBg_epeak - sig_subBg_gpeak
    
    
#    xlim = 0.02
#    index = np.where(np.logical_and.reduce((freq>-xlim, freq<xlim)))
#    freq = freq[index[0]]
#    spec = spec[:,index[0]]
#    specbg = specbg[:,index[0]]
    
    sindex = np.argsort(freq)
    freq = freq[sindex]
    spec = spec[:,sindex]
    specbg = specbg[:,sindex]
    
      
    mean_spec_subbg = np.mean(spec, axis = 0) - np.mean(specbg, axis = 0)
    
    mean_spec = np.mean(spec, axis = 0) 

    lim = np.max(sig)*1.2   
    binnum = 50
    
    bins_x = np.linspace(-lim, lim, binnum)
    bins_y = np.linspace(-lim, lim, binnum)
    
    
    
    if ifplot:
        plt.figure()
        plt.plot(freq, mean_spec_subbg, 'b-')
#        plt.legend(['g state','e state'])
        plt.xlabel('Freq (GHz)')
        plt.ylabel('FT signal')
        plt.xlim([-0.05,0.05])
#        plt.ylim([-50, 2000])
        autoscale_y(plt.gca())
        plt.show()
        
    if ifplot:
        plt.figure()
        plt.plot(freq, mean_spec, 'k-')
        
        
#        plt.legend(['g state','e state'])
        
#        plt.plot(freq, spec[0,:], 'b--')
#        plt.plot(freq, spec_estate[0,0,:], 'r--')
        
        plt.xlabel('Freq (GHz)')
        plt.ylabel('FT signal')
        plt.xlim([-0.05,0.05])
#        plt.ylim([-50, 5000])
        autoscale_y(plt.gca())
        plt.show()
    
    
    data = {"freq":freq, "mean_spec_subbg":mean_spec_subbg, "mean_spec":mean_spec}
    
    if ifsave:
        root = Tk()
        root.wm_attributes('-topmost', 1)
        savename = asksaveasfile(initialdir='C:\\Data\\2020\\201021_thermal_readout\\02_coherent_readout',parent=root)
        root.destroy() 
        io.savemat(savename.name, data)
    

def rotate_iq(angle_deg=0., rec=[None, None]):
    rec_rotated = [None, None]
    
    ch_cmplx = rec[0] + 1j*rec[1]
    ch_cmplx_rot = abs(ch_cmplx)*np.exp(1j*np.angle(ch_cmplx))*np.exp(1j*np.pi*angle_deg/180)
    rec_rotated[0] = np.real(ch_cmplx_rot)
    rec_rotated[1] = np.imag(ch_cmplx_rot)
    rec_rotated = np.array(rec_rotated)
    return rec_rotated

def process_readout_raw(rec_all_raw_ave, ifplot=0):
    # This program can be used to determine the following parameters
   
   
    rec_int = np.mean(rec_all_raw_ave[:,:,startindex:endindex],axis=2) - np.mean(rec_all_raw_ave[:,:,bgstartindex:bgendindex],axis=2)
#    rec_int = np.mean(rec_all_raw_ave[:,:,startindex:endindex],axis=2) - 123
    rec_rot = rotate_iq(angle, rec_int) 
    
    ch_a = rec_rot[0]
    ch_b = rec_rot[1]
    
   
    prob_g = np.size(np.where(ch_b<=143))/np.size(ch_b)
    prob_e = 1-prob_g
    
    if ifplot:
        ax = plt.figure()
        plt.plot(np.mean(rec_all_raw_ave[0][:,:],axis=0),'r')
        plt.plot(np.mean(rec_all_raw_ave[1][:,:],axis=0),'b')
        plt.legend(['ch0','ch1'])
        plt.xlabel('Time (ns)')
        plt.ylabel('ALAZAR read')
        plt.show()
#        plt.gca.set_aspect('equal')
        
#        
#        plt.figure()
#        ax_scatter = plt.axes()
#        ax_scatter.scatter(ch_a, ch_b, s=0.1, c='r', marker='.', alpha=0.1)
#
#        lim = np.max(np.concatenate((np.abs(ch_a),np.abs(ch_b))))*1.2
#        ax_scatter.set_xlim([-lim, lim])
#        ax_scatter.set_ylim([-lim, lim])
#        
#        ax_scatter.set_aspect(aspect=1)
#        plt.xlabel('ch0 read')
#        plt.ylabel('ch1 read')
#        plt.show()
#        
#        plt.figure()
#        binwidth = 1
#        bins_y = np.arange(-lim, lim, binwidth)
#        hist, bins = np.histogram(ch_b, bins=bins_y)
#        count = plt.hist(ch_b, bins=bins_y, histtype='step', orientation='vertical')
#        plt.plot(bins_y[0:-1]+binwidth/2, count[0], 'b-')
#        intcount = np.zeros(np.size(count[0]))
#        for k in np.arange(np.size(count[0])):
##            intcount[k] = np.sum(count[0][0:k+1])
#            intcount[k] = np.sum(hist[0:k+1])
##        print(co1unt)
#        plt.xlabel('ch1 read')
#        plt.ylabel('Samples')
#        plt.show()
#        
#        plt.figure()
#        plt.plot(bins_y[0:-1]+binwidth/2, intcount/np.max(intcount), 'b-')
#        plt.xlabel('ch1 read')
#        plt.ylabel('Integrated histogram')
#        plt.show()
        
        
    return rec_rot, prob_e

def process_readout_pattern(rec_all_raw_ave, daq_params, ifplot=0):
    # This program can be used to process the data with a pattern

    num_patterns = daq_params.num_patterns
    num_records_per_pattern = daq_params.num_records_per_pattern
    
    
#    # Readout Integration Parameters 
#    startindex = 0
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 7000
#    bgendindex = bgstartindex + duration
#    
#    # Rotation angle in iq plane (unit: degree)
#    angle = 245
#    # cutoff value to determine qubit in g or e state (run compare_ge_readout to get this value)
#    cutoff = 5
    
#    rec_all_raw_ave = np.zeros((np.shape(rec_all_raw)[1],
#                               np.shape(rec_all_raw)[0]*np.shape(rec_all_raw)[2],np.shape(rec_all_raw)[3]))
#    for k in np.arange(np.shape(rec_all_raw)[0]):
#        rec_all_raw_ave[:,k*np.shape(rec_all_raw)[2]:(k+1)*np.shape(rec_all_raw)[2],:] = rec_all_raw[k]
    
    
#    print(np.shape(rec_all_raw_ave))
    rec_int = np.mean(rec_all_raw_ave[:,:,startindex:endindex],axis=2) - np.mean(rec_all_raw_ave[:,:,bgstartindex:bgendindex],axis=2)
#    print(np.shape(rec_int))
    rec_rot = rotate_iq(angle, rec_int) 
    rec_rot = rec_rot[:,0:num_patterns*num_records_per_pattern]
    
    
    ch_a = rec_rot[0]
    ch_b = rec_rot[1]
    
    ch_b = np.reshape(ch_b,(num_records_per_pattern,num_patterns))
    
    prob_e = np.zeros(num_patterns)
    
    for k in np.arange(num_patterns):
        prob_e[k] = np.size(np.where(ch_b[:,k]>cutoff))/num_records_per_pattern
    
    if ifplot:
          
                
        plt.figure()
        plt.plot(np.arange(num_patterns), prob_e, 'k.-')
        plt.xlabel('data #')
        plt.ylabel('P$_{e}$')
        plt.show()
        
        
    return rec_rot, prob_e    


def process_coherentReadout_pattern(rec_all_raw_ave, daq_params, ifplot=0):
    # This program can be used to process the data with a pattern

    num_patterns = daq_params.num_patterns
    num_records_per_pattern = daq_params.num_records_per_pattern
    
    
#    # Readout Integration Parameters 
#    startindex = 0
#    duration = 1000
#    endindex = startindex + duration
#    bgstartindex = 7000
#    bgendindex = bgstartindex + duration
#    
#    # Rotation angle in iq plane (unit: degree)
#    angle = 245
#    # cutoff value to determine qubit in g or e state (run compare_ge_readout to get this value)
#    cutoff = 5
    
#    rec_all_raw_ave = np.zeros((np.shape(rec_all_raw)[1],
#                               np.shape(rec_all_raw)[0]*np.shape(rec_all_raw)[2],np.shape(rec_all_raw)[3]))
#    for k in np.arange(np.shape(rec_all_raw)[0]):
#        rec_all_raw_ave[:,k*np.shape(rec_all_raw)[2]:(k+1)*np.shape(rec_all_raw)[2],:] = rec_all_raw[k]
    
    
#    print(np.shape(rec_all_raw_ave))
    rec_int = np.mean(rec_all_raw_ave[:,:,startindex:endindex],axis=2) - np.mean(rec_all_raw_ave[:,:,bgstartindex:bgendindex],axis=2)
#    print(np.shape(rec_int))
    rec_rot = rotate_iq(angle, rec_int) 
    rec_rot = rec_rot[:,0:num_patterns*num_records_per_pattern]
    
    
    
    ch_a = rec_rot[0]
    ch_b = rec_rot[1]
    
    
    print([num_records_per_pattern,num_patterns])
    
    ch_b = np.reshape(ch_b,(num_records_per_pattern,num_patterns))
    
    prob_e = np.zeros(num_patterns)
    for k in np.arange(num_patterns):
        prob_e[k] = np.mean(ch_b[:,k])
    
    if ifplot:
          
                
        plt.figure()
        plt.plot(np.arange(num_patterns), prob_e, 'k.-')
        plt.xlabel('data #')
        plt.ylabel('P$_{e}$')
        plt.show()
        
        
    return rec_rot, prob_e 
    
def process_thermalreadout_pattern(rec_all_raw_ave, daq_params, ifplot=0):
    # This program can be used to process the data with a pattern

    num_patterns = daq_params.num_patterns
    num_records_per_pattern = daq_params.num_records_per_pattern
    
    
#    # Readout Integration Parameters 
#    startindex = 300
#    duration = 2000
#    endindex = startindex + duration
#    bgstartindex = 4000
#    bgendindex = bgstartindex + duration
    
    sig =  np.zeros(num_records_per_pattern*num_patterns)
    bg =  np.zeros(num_records_per_pattern*num_patterns)
    
    
    
    
    
    
    start_time = time.time()
    
    
    xdata = np.arange(np.shape(rec_all_raw_ave)[2])        
    index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
            
            
    totspec = np.zeros((num_patterns,np.size(index)))
    totbg = np.zeros((num_patterns,np.size(index)))
    for l in np.arange(num_patterns):
        for m in np.arange(num_records_per_pattern):
            k = m * num_patterns + l
            xdata = np.arange(np.shape(rec_all_raw_ave)[2])
            signal_i = rec_all_raw_ave[0][k,:] - 1*np.mean(rec_all_raw_ave[0][k,bgstartindex:bgendindex])
            signal_q = rec_all_raw_ave[1][k,:] - 1*np.mean(rec_all_raw_ave[1][k,bgstartindex:bgendindex])      
            
            index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
            xdata = xdata[index[0]]
            signal = signal_i[index[0]] + 1j * signal_q[index[0]]
            fourier = np.fft.fft(signal)
            n = signal.size
            timestep = xdata[2]-xdata[1]
            freq = np.fft.fftfreq(n, d=timestep)
            spectrum = np.abs(fourier)
            totspec[l,:] = totspec[l,:] + spectrum
            
            xdata = np.arange(np.shape(rec_all_raw_ave)[2])
            signal_i = rec_all_raw_ave[0][k,:] - 1*np.mean(rec_all_raw_ave[0][k,bgstartindex:bgendindex])
            signal_q = rec_all_raw_ave[1][k,:] - 1*np.mean(rec_all_raw_ave[1][k,bgstartindex:bgendindex]) 
            
            index = np.where(np.logical_and(xdata>bgstartindex,xdata<bgendindex))
            xdata = xdata[index[0]]
            signal = signal_i[index[0]] + 1j * signal_q[index[0]]
            fourier = np.fft.fft(signal)
            n = signal.size
            timestep = xdata[2]-xdata[1]
            freq = np.fft.fftfreq(n, d=timestep)
            
            spectrum = np.abs(fourier)
            totbg[l,:] = totbg[l,:] + spectrum
    
#    spec = (totspec - totbg)/num_records_per_pattern
    spec = totspec/num_records_per_pattern
    
    sig_subBg_gstate = np.zeros(np.shape(spec)[0])
    
    
#    fe = 0.0105
    
    fg = 0.0139 
    fe = 0.0025
    
    w = 0.002
    # integration region
    freqstart = fe - w
    freqend = fe + w
    findex =  np.where(np.logical_and(freq>freqstart,freq<freqend))[0]
    sig_subBg_estate = np.mean(spec[:,findex],axis=1)
    
    
    freqstart = -fe - w
    freqend = -fe + w
    findex =  np.where(np.logical_and(freq>freqstart,freq<freqend))[0]
    sig_subBg_estate = sig_subBg_estate + np.mean(spec[:,findex],axis=1)
    
    freqstart = fg - w
    freqend = fg + w
    findex =  np.where(np.logical_and(freq>freqstart,freq<freqend))[0]
    sig_subBg_gstate = np.mean(spec[:,findex],axis=1)
    
    
    freqstart = -fg - w
    freqend = -fg + w
    findex =  np.where(np.logical_and(freq>freqstart,freq<freqend))[0]
    sig_subBg_gstate = sig_subBg_estate + np.mean(spec[:,findex],axis=1)
    
    print("--- Analysis time: %s seconds ---" % (time.time() - start_time))
    
#    sig_subBg = sig - bg
#        
#    sig_subBg = np.reshape(sig_subBg,(num_records_per_pattern,num_patterns))
#    sig_subBg = np.mean(sig_subBg[:,:], axis = 0)
    
     
    if ifplot:
#          
#                
        plt.figure()
        plt.plot(np.arange(num_patterns), sig_subBg_gstate, 'b.-')
        plt.plot(np.arange(num_patterns), sig_subBg_estate, 'r.-')
        plt.xlabel('data #')
        plt.ylabel('FT signal')
        plt.legend(['integrate on g peak','integrate on e peak'],loc='upper left')
        plt.show()
        
        
        plt.figure()
#        plt.plot(np.arange(num_patterns), sig_subBg_gstate, 'b.-')
        plt.plot(np.arange(num_patterns), sig_subBg_estate-sig_subBg_gstate, 'r.-')
        plt.xlabel('data #')
        plt.ylabel('FT signal')
        plt.show()
        
        plt.figure()
        
        index = np.where(freq>0)
        freq = freq[index]
        
        plt.figure()
        plt.plot(freq,np.squeeze(totspec[2,index]),'b')
        plt.plot(freq,np.squeeze(totbg[2,index]),'r')
        plt.xlim([0,0.05])
        
    return sig_subBg_gstate, sig_subBg_estate


def process_thermalreadout_raw(rec_all_raw_ave, gpeak= 0.0145, epeak= 0.002, linewidth=0.001,  ifplot=0):
    # This program can be used to process the data with a pattern

    
    start_time = time.time()
    
    xdata = np.arange(np.shape(rec_all_raw_ave)[2])        
    index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
    
    totrep = np.shape(rec_all_raw_ave)[1]

#    meansig_i = np.mean(rec_all_raw_ave[0][:,:],axis=0)- 1*np.mean(rec_all_raw_ave[0][:,bgstartindex:bgendindex],axis=(0,1))
#    meansig_q = np.mean(rec_all_raw_ave[1][:,:],axis=0)- 1*np.mean(rec_all_raw_ave[0][:,bgstartindex:bgendindex],axis=(0,1))
#    meansig = meansig_i + 1j * meansig_q
#    fourier = np.fft.fft(meansig)
#    n = meansig.size
#    timestep = xdata[2]-xdata[1]
#    freq = np.fft.fftfreq(n, d=timestep)
#    spectrum = np.abs(fourier)
#    plt.figure()
#    plt.plot(freq, spectrum)
#    plt.ylim([0,1000])
    
    
    for k in np.arange(totrep):
        xdata = np.arange(np.shape(rec_all_raw_ave)[2])
        signal_i = rec_all_raw_ave[0][k,:] - 1*np.mean(rec_all_raw_ave[0][k,bgstartindex:bgendindex])
        signal_q = rec_all_raw_ave[1][k,:] - 1*np.mean(rec_all_raw_ave[1][k,bgstartindex:bgendindex])     
        
        index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
        xdata = xdata[index[0]]
        signal = signal_i[index[0]] + 1j * signal_q[index[0]]
#        signal = signal_q[index[0]]
        
        
        
        fourier = np.fft.fft(signal)
        n = signal.size
        timestep = xdata[2]-xdata[1]
        freq = np.fft.fftfreq(n, d=timestep)
        spectrum = np.abs(fourier)
#        spectrum = fourier
#        indexf = np.where(freq>=0)
#        freq = freq[indexf[0]]
#        spectrum = spectrum[indexf[0]]
        
        if k == 0:
            totspec = np.zeros((totrep,np.size(spectrum)))
            totbg = np.zeros((totrep,np.size(spectrum)))
            
#            plt.figure()
#            plt.plot(signal_i, 'b-')
#            plt.plot(signal_q, 'r-')
#            plt.plot(np.sqrt(signal_i**2+signal_q**2), 'g-')
#            plt.xlim([0,3000])
#            plt.show()
        
        totspec[k,:] = spectrum
        
        xdata = np.arange(np.shape(rec_all_raw_ave)[2])
        signal_i = rec_all_raw_ave[0][k,:] - 1*np.mean(rec_all_raw_ave[0][k,bgstartindex:bgendindex])
        signal_q = rec_all_raw_ave[1][k,:] - 1*np.mean(rec_all_raw_ave[1][k,bgstartindex:bgendindex])
        
        index = np.where(np.logical_and(xdata>bgstartindex,xdata<bgendindex))
        xdata = xdata[index[0]]
        signal = signal_i[index[0]] + 1j * signal_q[index[0]]
#        signal = signal_q[index[0]]
        
        fourier = np.fft.fft(signal)
        n = signal.size
        timestep = xdata[2]-xdata[1]
        freq = np.fft.fftfreq(n, d=timestep)
        spectrum = np.abs(fourier)
#        spectrum = fourier
#        indexf = np.where(freq>=0)
#        freq = freq[indexf[0]]
#        spectrum = spectrum[indexf[0]]
        
        totbg[k,:] = spectrum
    
    
    spec = totspec
    specbg = totbg
    
    sindex = np.argsort(freq)
    freq = freq[sindex]
    spec = spec[:,sindex]
    specbg = specbg[:,sindex]
    
#    spec = totspec2 - totbg2
    
    if ifplot:
        plt.figure()
        plt.plot(freq, spec[0,:], '-')
        plt.plot(freq, np.mean(spec[:,:],axis=0), '--')
        plt.xlabel('Freq (GHz)')
        plt.ylabel('FT signal')
        plt.xlim([0,0.02])
        plt.show()
        
        plt.figure()
        plt.plot(freq, totspec[0,:], 'r-')
        plt.plot(freq, totbg[0,:], 'b-')
        plt.xlabel('Freq (GHz)')
        plt.ylabel('FT signal')
        plt.xlim([0.001,0.02])
        autoscale_y(plt.gca())
        plt.show()
    
    # integration region
#    gpeak = 0.011 # 0 dB
#    epeak = 0.004
#    linewidth = 0.001
    
    
#    gpeak = 0.015 # 10 dB
#    epeak = 0.002    
#    
#    linewidth = 0.001
    
#    gpeak = 0.014 # 10 dB
#    epeak = 0.0023    
#    
#    linewidth = 0.0002
    
    
#    gpeak = 0.0145 # off
#    epeak = 0.002
#    linewidth = 0.0002
    
#    sig_subBg_gpeak = np.zeros(np.shape(spec)[0])
    freqstart_g = gpeak-linewidth
    freqend_g = gpeak+linewidth
    findex =  np.where(np.logical_and(freq>freqstart_g,freq<freqend_g))[0]
    sig_subBg_gpeak = np.mean(spec[:,findex],axis=1)
    
    freqstart_g = -gpeak-linewidth
    freqend_g = -gpeak+linewidth
    findex =  np.where(np.logical_and(freq>freqstart_g,freq<freqend_g))[0]
    sig_subBg_gpeak = sig_subBg_gpeak + np.mean(spec[:,findex],axis=1)
    
    
    freqstart_e = epeak-linewidth
    freqend_e = epeak+linewidth
    findex =  np.where(np.logical_and(freq>freqstart_e,freq<freqend_e))[0]
    sig_subBg_epeak = np.mean(spec[:,findex],axis=1)
    
    freqstart_e = -epeak-linewidth
    freqend_e = -epeak+linewidth
    findex =  np.where(np.logical_and(freq>freqstart_e,freq<freqend_e))[0]
    sig_subBg_epeak = sig_subBg_epeak + np.mean(spec[:,findex],axis=1)
    
    
    
    print("--- Analysis time: %s seconds ---" % (time.time() - start_time))
        
    return sig_subBg_gpeak, sig_subBg_epeak, spec, specbg, freq


def process_thermalreadout_raw_weight(rec_all_raw_ave, ifplot=0):
    # This program can be used to process the data with a pattern

    
    start_time = time.time()
    
    xdata = np.arange(np.shape(rec_all_raw_ave)[2])        
    index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
    
    totrep = np.shape(rec_all_raw_ave)[1]

#    meansig_i = np.mean(rec_all_raw_ave[0][:,:],axis=0)- 1*np.mean(rec_all_raw_ave[0][:,bgstartindex:bgendindex],axis=(0,1))
#    meansig_q = np.mean(rec_all_raw_ave[1][:,:],axis=0)- 1*np.mean(rec_all_raw_ave[0][:,bgstartindex:bgendindex],axis=(0,1))
#    meansig = meansig_i + 1j * meansig_q
#    fourier = np.fft.fft(meansig)
#    n = meansig.size
#    timestep = xdata[2]-xdata[1]
#    freq = np.fft.fftfreq(n, d=timestep)
#    spectrum = np.abs(fourier)
#    plt.figure()
#    plt.plot(freq, spectrum)
#    plt.ylim([0,1000])
    
    
    for k in np.arange(totrep):
        xdata = np.arange(np.shape(rec_all_raw_ave)[2])
        signal_i = rec_all_raw_ave[0][k,:] - 1*np.mean(rec_all_raw_ave[0][k,bgstartindex:bgendindex])
        signal_q = rec_all_raw_ave[1][k,:] - 1*np.mean(rec_all_raw_ave[1][k,bgstartindex:bgendindex])     
        
        index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
        xdata = xdata[index[0]]
        signal = signal_i[index[0]] + 1j * signal_q[index[0]]
#        signal = signal_q[index[0]]
        
        
        
        fourier = np.fft.fft(signal)
        n = signal.size
        timestep = xdata[2]-xdata[1]
        freq = np.fft.fftfreq(n, d=timestep)
        spectrum = np.abs(fourier)
#        spectrum = fourier
#        indexf = np.where(freq>=0)
#        freq = freq[indexf[0]]
#        spectrum = spectrum[indexf[0]]
        
        if k == 0:
            totspec = np.zeros((totrep,np.size(spectrum)))
            totbg = np.zeros((totrep,np.size(spectrum)))
            
#            plt.figure()
#            plt.plot(signal_i, 'b-')
#            plt.plot(signal_q, 'r-')
#            plt.plot(np.sqrt(signal_i**2+signal_q**2), 'g-')
#            plt.xlim([0,3000])
#            plt.show()
        
        totspec[k,:] = spectrum
        
        xdata = np.arange(np.shape(rec_all_raw_ave)[2])
        signal_i = rec_all_raw_ave[0][k,:] - 1*np.mean(rec_all_raw_ave[0][k,bgstartindex:bgendindex])
        signal_q = rec_all_raw_ave[1][k,:] - 1*np.mean(rec_all_raw_ave[1][k,bgstartindex:bgendindex])
        
        index = np.where(np.logical_and(xdata>bgstartindex,xdata<bgendindex))
        xdata = xdata[index[0]]
        signal = signal_i[index[0]] + 1j * signal_q[index[0]]
#        signal = signal_q[index[0]]
        
        fourier = np.fft.fft(signal)
        n = signal.size
        timestep = xdata[2]-xdata[1]
        freq = np.fft.fftfreq(n, d=timestep)
        spectrum = np.abs(fourier)
#        spectrum = fourier
#        indexf = np.where(freq>=0)
#        freq = freq[indexf[0]]
#        spectrum = spectrum[indexf[0]]
        
        totbg[k,:] = spectrum
    
    
    spec = totspec
    specbg = totbg
    
    sindex = np.argsort(freq)
    freq = freq[sindex]
    spec = spec[:,sindex]
    specbg = specbg[:,sindex]
    
#    spec = totspec2 - totbg2
    
    if ifplot:
        plt.figure()
        plt.plot(freq, spec[0,:], '-')
        plt.plot(freq, np.mean(spec[:,:],axis=0), '--')
        plt.xlabel('Freq (GHz)')
        plt.ylabel('FT signal')
        plt.xlim([0,0.02])
        plt.show()
        
        plt.figure()
        plt.plot(freq, totspec[0,:], 'r-')
        plt.plot(freq, totbg[0,:], 'b-')
        plt.xlabel('Freq (GHz)')
        plt.ylabel('FT signal')
        plt.xlim([0.001,0.02])
        autoscale_y(plt.gca())
        plt.show()
    
        
    print("--- Analysis time: %s seconds ---" % (time.time() - start_time))
        
    return spec, specbg, freq

def process_thermalreadout_raw_notFFT(rec_all_raw_ave, ifplot=0):
    # This program can be used to process the data with a pattern

    
    start_time = time.time()
    
    xdata = np.arange(np.shape(rec_all_raw_ave)[2])        
    index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
    
    totrep = np.shape(rec_all_raw_ave)[1]

        
    totsig = np.zeros(totrep)
    totbg = np.zeros(totrep)
    
    totsig2 = np.zeros(totrep)
    totbg2 = np.zeros(totrep)
    
    omega = 0.004
    
    for k in np.arange(totrep):
        xdata = np.arange(np.shape(rec_all_raw_ave)[2])
        signal_i = rec_all_raw_ave[0][k,:] - np.mean(rec_all_raw_ave[0][k,bgstartindex:bgendindex])
        signal_q = rec_all_raw_ave[1][k,:] - np.mean(rec_all_raw_ave[1][k,bgstartindex:bgendindex])       
        
        index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
        xdata = xdata[index[0]]
#        signal = signal_i[index[0]] 
        signal = signal_i[index[0]] + 1j * signal_q[index[0]]
        
        FTsig = np.abs(np.sum(signal*np.exp(-2*np.pi*1j*omega*xdata)))
 
        totsig[k] = FTsig
        
#        xdata = np.arange(np.shape(rec_all_raw_ave)[2])
#        signal_i = rec_all_raw_ave[0][k,:]
#        signal_q = rec_all_raw_ave[1][k,:]
#        
#        index = np.where(np.logical_and(xdata>bgstartindex,xdata<bgendindex))
#        xdata = xdata[index[0]]
#        signal = signal_q[index[0]]
#        signal = signal_i[index[0]] + 1j * signal_q[index[0]]
#        
#        
#        FTbg = np.abs(np.sum(signal*np.exp(-2*np.pi*1j*omega*xdata)))
# 
#        totbg[k] = FTbg
        
#    for k in np.arange(totrep):
#        xdata = np.arange(np.shape(rec_all_raw_ave)[2])
#        signal_i = rec_all_raw_ave[0][k,:]
#        signal_q = rec_all_raw_ave[1][k,:]        
#        
#        index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
#        xdata = xdata[index[0]]
#        signal = signal_q[index[0]] 
#        
#        FTsig = np.abs(np.sum(signal*np.exp(-2*np.pi*1j*omega*xdata)))
# 
#        totsig2[k] = FTsig
#        
#        xdata = np.arange(np.shape(rec_all_raw_ave)[2])
#        signal_i = rec_all_raw_ave[0][k,:]
#        signal_q = rec_all_raw_ave[1][k,:]
#        
#        index = np.where(np.logical_and(xdata>bgstartindex,xdata<bgendindex))
#        xdata = xdata[index[0]]
#        signal = signal_q[index[0]]
#        
#        
#        FTbg = np.abs(np.sum(signal*np.exp(-2*np.pi*1j*omega*xdata)))
# 
#        totbg2[k] = FTbg
    
    sig = totsig

    
    print("--- Analysis time: %s seconds ---" % (time.time() - start_time))
        
    return sig


def process_thermalreadout_raw_noFourierTransform(rec_all_raw_ave, ifplot=0):
    # This program can be used to process the data with a pattern

    
    start_time = time.time()
    
    xdata = np.arange(np.shape(rec_all_raw_ave)[2])        
    index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
    
    totrep = np.shape(rec_all_raw_ave)[1]

        
    
    totsig = np.zeros(totrep)
    for k in np.arange(totrep):
        xdata = np.arange(np.shape(rec_all_raw_ave)[2])
        signal_i = rec_all_raw_ave[0][k,:] - 1*np.mean(rec_all_raw_ave[0][k,bgstartindex:bgendindex])
        signal_q = rec_all_raw_ave[1][k,:] - 1*np.mean(rec_all_raw_ave[1][k,bgstartindex:bgendindex])     
        
        index = np.where(np.logical_and(xdata>startindex,xdata<endindex))
        xdata = xdata[index[0]]
        signal = np.sum(np.sqrt(signal_i[index[0]]**2 + signal_q[index[0]]**2))
        totsig[k] = signal
        
    
    print("--- Analysis time: %s seconds ---" % (time.time() - start_time))
        
    return totsig

def program_stop_button():
    def start():
        """Enable scanning by setting the global flag to True."""
        global running
        running = True

    def stop():
        """Stop scanning by setting the global flag to False."""
        global running
        running = False
        
    def scanning():
        if running:  # Only do this if the Stop button has not been clicked
            print('start scan')
            readout_singleRun(ifWX = 0)
            root.after(10, scanning)
            
        # After 1 second, call scanning again (create a recursive loop)
        else:
            print('stop')
        
#    print("Initialize WX")
#    wx_programs.wx_initialize()
    global running
    running = True    
    root = Tk()
    root.wm_attributes('-topmost', 1)
    
    root.title("Title")
    root.geometry("500x500")
    
    app = tki.Frame(root)
    app.grid()
    
    start = tki.Button(app, text="Start Scan", command=scanning)
    stop = tki.Button(app, text="Stop", command=stop)
    
    
    start.grid()
    stop.grid()

    root.mainloop()

    
            
   
#    root.destroy() 

def Gaussian(x, A, c, sigma):
    return A * np.exp(-(x-c)**2/(2*sigma**2))

# Disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore print
def enablePrint():
    sys.stdout = sys.__stdout__   


def autoscale_y(ax,margin=0.05):
    """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""

    import numpy as np

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo,hi = ax.get_xlim()
        y_displayed = yd[((xd>=lo) & (xd<=hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed)-margin*h
        top = np.max(y_displayed)+margin*h
        return bot,top

    lines = ax.get_lines()
    bot,top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot: bot = new_bot
        if new_top > top: top = new_top

    ax.set_ylim(bot,top)
    
if __name__ == '__main__':
    pass
  #  pulsed_readout_ramsey_GE()
    #rabi_gaussian()
    #Ramsey()
    #T2_echo()
    #T1()
  #  ramsey_GE()
    #T1_constant()
   # rabi_constant()
    #ramsey_constant()
    #ramseyef_constant()
    #ramseygf_constant()
    #echogf_constant()
    #rabi()
    #ramsey_GE()