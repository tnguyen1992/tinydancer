# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:11:46 2023

@author: tnguyen
"""

 from mne.fiff import read_evoked
    evoked = read_evoked('ctf-ave.fif')
    print(type(evoked))