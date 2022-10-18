#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 15 16:46:40 2022

@author: james
"""

from astropy.io.votable import parse
import urllib.request 
import numpy as np
import matplotlib.pyplot as plt


class photometricSED():
    def get_url_CDS(RADec, rad):
        RADec_split = RADec.split(":")
        RADec_new_split=RADec_split[2].split(" ")
        url="http://vizier.u-strasbg.fr/viz-bin/sed?-c="
        url+=RADec_split[0]
        url+="%20"
        url+=RADec_split[1]
        url+="%20"
        url+=RADec_new_split[0]
        if "+" in RADec:
            url+="%20%2B"
            url+=RADec_new_split[1][0:]
        elif "-" in RADec:
            url+="%20"
            url+=RADec_new_split[1]
        url+="%20"
        url+=RADec_split[3]
        url+="%20"
        url+=RADec_split[4]
        url+="&-c.rs="
        url+=str(rad)
        return url
    
    
    
    
    def plot_SED(url):
        urllib.request.urlretrieve(url, "photSED.vot")
        votab = parse("photSED.vot")
        
        
        for table in votab.iter_tables():
            data = table.array
            #print(data)
            ra=(data["_RAJ2000"])
            dec=(data["_DEJ2000"])
            sedfreq=data["sed_freq"]
            sedflux=data["sed_flux"]
            sedfluxe=data["sed_eflux"]
            sedfilter=data["sed_filter"]
            
            #log_sedfluxe=1/np.log(10) * sedfluxe/sedflux
        
        
        sed_wl = 299792458/sedfreq  # c/f to get wl
        
        mask=sed_wl<1400
        ra=ra[mask]
        dec=dec[mask]
        sedfreq=sedfreq[mask]
        sedflux=sedflux[mask]
        sedfluxe=sedfluxe[mask]
        sedfilter=sedfilter[mask]
        sed_wl=sed_wl[mask]
        
        unique_filt = np.unique(sedfilter)
        
        colours=[]
        for i in unique_filt:
            colours.append(np.random.rand(3,))
        
        plt.clf()
        for count, filt in enumerate(sedfilter):
            a=np.argwhere(unique_filt==filt)[0]
            col=colours[np.argwhere(unique_filt==filt)[0][0]]
            plt.errorbar(sed_wl[count], (sedflux[count]), yerr=sedfluxe[count], fmt='.', c=col,label=sedfilter[count])
            
            
        plt.grid()
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("log Flux [Jy]")
        plt.yscale('log')
        plt.savefig("PhotSED.pdf")
        #plt.legend()
        plt.clf()
        
        np.savetxt("SED_CDS.csv", np.array([sed_wl, sedflux, sedfluxe, sedfilter ]).T, fmt="%s")
    




