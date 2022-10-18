#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 09:46:24 2022

@author: james
"""

import requests
import numpy as np
import os

class PTF(object):
    def queryPTF(RA,Dec,rad='1'):
        # you can use astropy IRSA to do this, but I found that it timed out more often for some reason. 
        # WISE currently uses this, but I might change in future to a manual search
        # https://astroquery.readthedocs.io/en/latest/ipac/irsa/irsa.html
        
        RA=RA.split(":")
        Dec=Dec.split(":")
        
        
        if len(RA[0])==1:
            RA[0]="0"+RA[0]
        if len(RA[1])==1:
            RA[1]="0"+RA[1]
        CheckSec=RA[2].split(".")
        if len(CheckSec[0])==1:
            RA[2]="0"+RA[2]
        
        if len(Dec[0])==1:
            Dec[0]="0"+Dec[0]
        if len(Dec[1])==1:
            Dec[1]="0"+Dec[1]
        CheckSec=Dec[2].split(".")
        if len(CheckSec[0])==1:
            Dec[2]="0"+Dec[2]
        
        
        requests_string = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog=ptf_lightcurves&spatial=cone&radius='+rad+'&radunits=arcsec&objstr='+RA[0]+'h+'+RA[1]+'m+'+RA[2]+'s+'+Dec[0]+'d+'+Dec[1]+'m+'+Dec[2]+'s&outfmt=1'
        res = requests.get(requests_string)#selcols=ra,dec,w1mpro,w1sigmpro,w1snr,w2mpro,w2sigmpro,w2snr,w3mpro,w3sigmpro,w3snr,w4mpro,w4sigmpro,w4snr')
        
        with open("PTF.txt", "w") as f:
            f.write(res.text)
        
        
    def splitRandG_PTF():
        if "PTF.txt" in os.listdir(os.getcwd()):
            b = np.loadtxt("PTF.txt", unpack=True, skiprows=101, dtype=str)
            try:
                qual = b[-4]=="1"
                obsmjd=b[0].astype(np.float64)[qual]
                mag_autocorr=b[1].astype(np.float64)[qual]
                magerr_auto=b[2].astype(np.float64)[qual]
                
                r=b[8]=="2"
                g=b[8]=="1"
                if len(obsmjd[r])>0:
                    np.savetxt("rPTF.csv", np.array([obsmjd[r], mag_autocorr[r], magerr_auto[r]]).T)
                if len(obsmjd[g])>0:
                    np.savetxt("gPTF.csv", np.array([obsmjd[g], mag_autocorr[g], magerr_auto[g]]).T)
            
            except: None


