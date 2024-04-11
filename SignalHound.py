# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:57:16 2024

@author: crow104
"""

import numpy as np
import pyvisa 
rm = pyvisa.ResourceManager()

#Copied Code from BNC, and changed name of variables, now trying to find right address and parameters

def set_SignalHound_output(freq_GHz, power_dBm=13,\
                   SignalHound_addr='USB0::0x03EB::0xAFFF::141-216340000-0292::INSTR'):
    try:
        SignalHound_handle = rm.open_resource(SignalHound_addr)
        SignalHound_handle.write(f"sour:freq {freq_GHz*1E9}")
        SignalHound_handle.write(f"sour:pow {power_dBm}")
    finally:
        SignalHound_handle.close()
##END set_bnc_output
    
    
def get_SignalHound_freq_GHz(SignalHound_addr='USB0::0x03EB::0xAFFF::141-216340000-0292::INSTR'):
    
    try:
        SignalHound_handle = rm.open_resource(SignalHound_addr)
        freq = SignalHound_handle.query(f"sour:freq?")
        return np.float(freq)/1E9
    finally:
        SignalHound_handle.close()
##END get_bnc_freq
    

def configure_for_sweep(start_freq_GHz, stop_freq_GHz, power_dBm, \
                        dwell_sec=0.05, num_freq_points=201,\
                        SignalHound_addr='USB0::0x03EB::0xAFFF::421-4385A0002-0619::INSTR'):
    try:
        SignalHound_handle = rm.open_resource(SignalHound_addr)
        
        print("Configuring BNC for sweep")
        # 0. Rest instrument
        SignalHound_handle.write("*CLS")
        print("0, ", SignalHound_handle.query(":SYST:ERR?") )
        
        # 1. configure RF out
        SignalHound_handle.write(f"sour:freq {start_freq_GHz*1E9}")
        SignalHound_handle.write(f"sour:powe {power_dBm}")
        print("1, ", bnc_handle.query(":SYST:ERR?") )
        
        # 2. configure RF sweep
        SignalHound_handle.write(":SWE:DEL 0")
        SignalHound_handle.write(":SWE:DIR UP;  :SWE:SPAC LIN;  :SWE:COUN 1.0")
        SignalHound_handle.write(f":SWE:DWEL {dwell_sec};  :SWE:POIN {num_freq_points}")
        print("2, ", SignalHound_handle.query(":SYST:ERR?") )

        # 3. configure to trigger internally
        SignalHound_handle.write(":TRIG:TYPE NORM;  :TRIG:SLOP POS;  :TRIG:SOUR IMM;")
        SignalHound_handle.write(":TRIG:DEL  0.0;   :TRIG:ECO 1.0") # last command: "run (E)very Nth count"
        SignalHound_handle.write(":TRIG:OUTP:MODE POIN;   :TRIG:OUTP:POL NORM;")
        print("3, ", SignalHound_handle.query(":SYST:ERR?") )

        # 4. setup trigger to VNA (Low Frequency Output)
        SignalHound_handle.write(":LFO:STAT ON;  :LFO:SOUR TRIG;  :LFO:SHAP SQU;")
        SignalHound_handle.write(":LFO:FREQ 0.0;  :LFO:AMPL 1.0")
        print("4, ", SignalHound_handle.query(":SYST:ERR?") )
        
        # 5. setup frequency sweep
        SignalHound_handle.write("INIT:CONT OFF;  :FREQ:MODE FIX;  :FREQ:MODE SWE")
        SignalHound_handle.write(f":FREQ:STAR {start_freq_GHz*1E9};  :FREQ:STOP {stop_freq_GHz*1E9}")
        print("5, ", SignalHound_handle.query(":SYST:ERR?") )

        # 6. prepare to run
        SignalHound_handle.write(":OUTP ON")
        SignalHound_handle.write("INIT:CONT OFF") # do not repeat sweep
        SignalHound_handle.write(":INIT") # arm trigger
        print("6, ", SignalHound_handle.query(":SYST:ERR?") )
        
    finally:
        SignalHound_handle.close()
##END configure_for_Sweep
#
#    bnc_handle.write("*CLS;")
#    print("1, ", bnc_handle.query(":SYST:ERR?") )
#    bnc_handle.write(":FREQ 4.000000 GHZ;:POW -15.000000;")
#    print("2, ", bnc_handle.query(":SYST:ERR?") )
#    bnc_handle.write(":SWE:DEL  0.000000;:SWE:DWEL 0.050000;:SWE:DIR UP;:SWE:SPAC LIN;:SWE:POIN 251.000000;:SWE:COUN 1.000000")
#    print("3, ", bnc_handle.query(":SYST:ERR?") )
#    bnc_handle.write(":TRIG:TYPE NORM;:TRIG:SLOP POS;:TRIG:SOUR IMM;:TRIG:DEL 0.000000;:TRIG:ECO 1.000000")
#    print("4, ", bnc_handle.query(":SYST:ERR?") )
#    bnc_handle.write(":TRIG:OUTP:MODE POIN;:TRIG:OUTP:POL NORM;")
#    print("5, ", bnc_handle.query(":SYST:ERR?") )
#    bnc_handle.write(":LFO:STAT ON;:LFO:SOUR TRIG;:LFO:SHAP SQU;:LFO:FREQ 0.000000;:LFO:AMPL 1.000000")
#    print("6, ", bnc_handle.query(":SYST:ERR?") )
#    bnc_handle.write(":OUTP ON;")
#    print("7, ", bnc_handle.query(":SYST:ERR?") )
#    bnc_handle.write("INIT:CONT OFF;:FREQ:MODE FIX;:FREQ:STAR 4.000000 GHZ;:FREQ:STOP 4.500000 GHZ;:FREQ:MODE SWE")
#    print("8, ", bnc_handle.query(":SYST:ERR?") )
#    bnc_handle.write("INIT:CONT OFF;:INIT")
#    print("9, ", bnc_handle.query(":SYST:ERR?") )