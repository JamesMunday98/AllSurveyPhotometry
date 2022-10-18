#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 29 17:46:57 2022

@author: james
"""


import numpy as np
import matplotlib.pyplot as plt
from astroquery.ipac.irsa import Irsa
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation
from miscAstro import * # analysis:ignore

#WISE 5σ photometric sensitivity is estimated to be 0.068, 0.098, 0.86 and 5.4 mJy 
#(16.6, 15.6, 11.3, 8.0 Vega mag) at 3.4, 4.6, 12 and 22 μm in unconfused regions on the 
#ecliptic plane


# PLEASE READ THIS if you wonder why I change the MJD https://wise2.ipac.caltech.edu/docs/release/neowise/expsup/sec1_2.html#exptime_offset
# new data release every spring https://irsa.ipac.caltech.edu/Missions/wise.html
class NEOWISE(object):
    def queryWise(RAdeg,Decdeg):
        Irsa.ROW_LIMIT = 5000   # 5000 is the new value for row limit here.
        Irsa.TIMEOUT = 200
        t = Irsa.query_region(SkyCoord(RAdeg, Decdeg, unit=(u.deg,u.deg)), spatial="Cone",
                                                            radius=1.5 * u.arcsec, catalog="neowiser_p1bs_psd", 
                                                            selcols="ra,dec,mjd,w1mpro,w1sigmpro,w2mpro,w2sigmpro") # these selcols are ignored only here?
        
        t.write('NEOWISE.csv', overwrite=True)


    def saveFilters(RAdeg,Decdeg,gmag):
        table=Table.read('NEOWISE.csv')
        
        #IMPORTANT NOTE: WISE WAS A SPACECRAFT AND THE STORED TIME WAS MJD. I DO NOT KNOW THE EXACT LOCATION FOR BJD
        
        if gmag<16.6: #not optimal, but if the gmag is 16.6 and i look for white dwarfs, will be brighter in blue
            try:
                mask1=((table['w1mpro'] > 0)      &     (table['w1mpro'] > gmag-2.25))
                mjd=table['mjd'][mask1]       + 0.57    # this is correct for WISE and NEOWISE W1 and W2
                
                for NEOcount, mjd_indiv in enumerate(mjd):
                    if mjd_indiv < 57196.48882: # THIS IS ONLY FOR NEOWISE - "Leap Second Error"
                        mjd[NEOcount]+=1/86400
                
                bjd=miscAstro.jd_corr(mjd,RAdeg,Decdeg,loc=EarthLocation.of_site("lapalma")).value
                w1mag=table['w1mpro'][mask1]
                w1mage=table['w1sigmpro'][mask1]
                np.savetxt("NEOWISE1_340nm.csv", np.array([bjd, w1mag, w1mage]).T)
            except Exception as e: print("exception here Neowise"); print(e)
        
        if gmag<15.6:
            try:
                mask2=((table['w2mpro'] > 0)      &     (table['w2mpro'] > gmag-2))
                mjd=table['mjd'][mask2]       + 0.57    # this is correct for WISE and NEOWISE W1 and W2
                
                for NEOcount, mjd_indiv in enumerate(mjd):
                    if mjd_indiv < 57196.48882: # THIS IS ONLY FOR NEOWISE - "Leap Second Error" in 2012
                        mjd[NEOcount]+=1/86400
                
                bjd=miscAstro.jd_corr(mjd,RAdeg,Decdeg,loc=EarthLocation.of_site("lapalma")).value
                w2mag=table['w2mpro'][mask2]
                w2mage=table['w2sigmpro'][mask2]
                np.savetxt("NEOWISE2_460nm.csv", np.array([bjd, w2mag, w2mage]).T)
            except Exception as e: print("exception here2 Neowise"); print(e)
        

    def Plot():
        try:
            mjd,mag,mage=np.loadtxt("NEOWISE1_340nm.csv",unpack=True)
            plt.scatter(mjd,mag, c='b', label="WISE1")
            plt.errorbar(mjd,mag,yerr=mage,ls=" ", c='b')
        except: None
        
        try:
            mjd,mag,mage=np.loadtxt("NEOWISE2_460nm.csv",unpack=True)
            plt.scatter(mjd,mag, c='g', label="WISE2")
            plt.errorbar(mjd,mag,yerr=mage,ls=" ", c='g')
        except: None
        



