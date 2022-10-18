#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 29 12:56:02 2022

@author: james
"""


import numpy as np
import matplotlib.pyplot as plt
from getZTF import ZTF
from getPanstarrs import getPanstarrs
from getWISE import WISE
from getASASSNweb import getASASSNweb
from astropy.stats import sigma_clip

class plotEverything(object):
    def plot():
        plt.clf()
        try: getPanstarrs.plot()
        except: None
        
        try: ZTF.plot_photometry("", filename="ZTF.csv", saveit=False)
        except: None
        
        try: #this is mjd, but it doesn't matter for plotting purposes. i should convert Panstars.csv times to BJD and then save rPTF as BJD in the future
            obsmjd, mag_autocorr, magerr_auto = np.loadtxt("rPTF.csv", unpack=True)
            plt.scatter(obsmjd, mag_autocorr)
            plt.errorbar(obsmjd, mag_autocorr, yerr=magerr_auto, ls=" ", label="rPTF")
        except: None
        
        try:
            BJD_V, mag_V, mag_err_V = np.loadtxt("ASASSNv_lc.dat", unpack=True)
            plt.scatter(BJD_V,mag_V)
            plt.errorbar(BJD_V,mag_V, yerr=mag_err_V, ls=" ", label="vASASSN")
        except: None
        
        try:
            BJD_G, mag_G, mag_err_G = np.loadtxt("ASASSNg_lc.dat", unpack=True)
            plt.scatter(BJD_G,mag_G)
            plt.errorbar(BJD_G,mag_G, yerr=mag_err_G, ls=" ", label="gASASSN")
        except: None
        
        try: 
            BJD, Mag, Magerr= np.genfromtxt('CatalinaDataBJD.csv', unpack=True)
            plt.scatter(BJD, Mag)
            plt.errorbar(BJD,Mag,yerr=Magerr,ls=' ', label="vCat")
        except: None
        
        try:
            mjdFileo, magFileo, dmagFileo, uJyFileo, duJyFileo = np.loadtxt("ATLAS_filtO_flux_and_err.dat", unpack=True)
            plt.scatter(mjdFileo, magFileo)
            plt.errorbar(mjdFileo, magFileo, yerr=dmagFileo, ls=" ", label="ATLASo")
        except: None
        
        try:
            mjdFilec, magFilec, dmagFilec, uJyFilec, duJyFilec = np.loadtxt("ATLAS_filtC_flux_and_err.dat", unpack=True)
            plt.scatter(mjdFilec, magFilec)
            plt.errorbar(mjdFilec, magFilec, yerr=dmagFilec, ls=" ", label="ATLASc")
        except: None
        
        try: 
            BJD, mag, magerr, filtASASSN = getASASSNweb.processDownloadedWebFile(0,0)
            if len(BJD) > 20:
                filtered_data = sigma_clip(mag, cenfunc='median', sigma_lower=5, sigma_upper=5, maxiters=2, masked=True)
                mask=mag==filtered_data
                BJD=BJD[mask]
                mag=mag[mask]
                magerr=magerr[mask]
                filtASASSN=filtASASSN[mask]
                
                for y in np.unique(filtASASSN):
                    plt.scatter(BJD[filtASASSN==y], mag[filtASASSN==y])
                    plt.errorbar(BJD[filtASASSN==y], mag[filtASASSN==y], yerr=magerr[filtASASSN==y], ls=" ", label=str(y))
                
        except Exception as e: print(e)
            
        
        WISE.Plot()
        
        plt.gca().invert_yaxis()
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        plt.xlabel("BJD")
        plt.ylabel("Magnitude")
        plt.autoscale()
        plt.savefig("AllPhot.png")
        plt.clf()
        plt.close()
        



