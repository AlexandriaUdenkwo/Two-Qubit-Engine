# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 18:52:34 2020

@author: crow104
"""

import sys, os

def blockPrint():
    sys.stdout = open(os.devnull, 'w')
    
def enablePrint():
    sys.stdout = sys.__stdout__
    
    
    
print('ss')
blockPrint()
print('2')
enablePrint()
print('12')
