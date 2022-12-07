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
import dateutil
import astropy
#WISE 5σ photometric sensitivity is estimated to be 0.068, 0.098, 0.86 and 5.4 mJy 
#(16.6, 15.6, 11.3, 8.0 Vega mag) at 3.4, 4.6, 12 and 22 μm in unconfused regions on the 
#ecliptic plane



# PLEASE READ THIS if you wander why I change the MJD https://wise2.ipac.caltech.edu/docs/release/neowise/expsup/sec1_2.html#exptime_offset
class WISE(object):
    def queryWise(RAdeg,Decdeg, ref_epoch=2020,pmra=0,pmdec=0):
        ############## correct for proper motion (rough, but better than nothing)
        dt = dateutil.parser.parse(str(int(ref_epoch))+'.01.01')
        dt2=dateutil.parser.parse('2010.06.01') #WISE was 2010, use a guess that the mean location was that at 2010.5
        time = astropy.time.Time(dt)
        time2= astropy.time.Time(dt2)
        deltamjd_yr=(time.jd - time2.jd)/365
        
        
        pos_ra_change = ((pmra/3600) /1000)  *(-1*deltamjd_yr)
        RAdeg+=pos_ra_change
        
        pos_dec_change = ((pmdec/3600) /1000)  *(-1*deltamjd_yr)
        Decdeg+=pos_dec_change
        ##############
        
        Irsa.ROW_LIMIT = 5000   # 1000 is the new value for row limit here.
        t = Irsa.query_region(SkyCoord(RAdeg, Decdeg, unit=(u.deg,u.deg)), spatial="Cone",
                                                            radius=3 * u.arcsec, catalog="allwise_p3as_mep", 
                                                            selcols="ra,dec,mjd,w1mpro_ep,w1sigmpro_ep,w2mpro_ep,w2sigmpro_ep,w3mpro_ep,w3sigmpro_ep,w4mpro_ep,w4sigmpro_ep")
        
        t.write('WISE.csv', overwrite=True)


    def saveFilters(RAdeg,Decdeg,gmag):
        table=Table.read('WISE.csv')
        
        #IMPORTANT NOTE: WISE WAS A SPACECRAFT AND THE STORED TIME WAS MJD. I DO NOT KNOW THE EXACT LOCATION FOR BJD
        
        if gmag<16.6: #not optimal, but if the gmag is 16.6 and i look for white dwarfs, will be brighter in blue
            try:
                mask1=((table['w1mpro_ep'] > 0)      &     (table['w1mpro_ep'] > gmag-2.25))
                mjd=table['mjd'][mask1]     + 0.57    # this is correct for WISE and NEOWISE W1 and W2
                bjd=miscAstro.jd_corr(mjd,RAdeg,Decdeg,loc=EarthLocation.of_site("lapalma")).value
                w1mag=table['w1mpro_ep'][mask1]
                w1mage=table['w1sigmpro_ep'][mask1]
                np.savetxt("WISE1_340nm.csv", np.array([bjd, w1mag, w1mage]).T)
            except: None
        
        if gmag<15.6:
            try:
                mask2=((table['w2mpro_ep'] > 0)      &     (table['w2mpro_ep'] > gmag-2))
                mjd=table['mjd'][mask2]      + 0.57    # this is correct for WISE and NEOWISE W1 and W2
                bjd=miscAstro.jd_corr(mjd,RAdeg,Decdeg,loc=EarthLocation.of_site("lapalma")).value
                w2mag=table['w2mpro_ep'][mask2]
                w2mage=table['w2sigmpro_ep'][mask2]
                np.savetxt("WISE2_460nm.csv", np.array([bjd, w2mag, w2mage]).T)
            except: None
        
        if gmag<10.5:
            try:
                mask3=((table['w3mpro_ep'] > 0)      &     (table['w3mpro_ep'] > gmag-2))
                mjd=table['mjd'][mask3]
                bjd=miscAstro.jd_corr(mjd,RAdeg,Decdeg,loc=EarthLocation.of_site("lapalma")).value
                w3mag=table['w3mpro_ep'][mask3]
                w3mage=table['w3sigmpro_ep'][mask3]
                np.savetxt("WISE3_1200nm.csv", np.array([bjd, w3mag, w3mage]).T)
            except: None
        
        if gmag<7:
            try:
                mask4=((table['w4mpro_ep'] > 0)      &     (table['w4mpro_ep'] > gmag-2))
                mjd=table['mjd'][mask4]
                bjd=miscAstro.jd_corr(mjd,RAdeg,Decdeg,loc=EarthLocation.of_site("lapalma")).value
                w4mag=table['w4mpro_ep'][mask4]
                w4mage=table['w4sigmpro_ep'][mask4]
                np.savetxt("WISE4_2200nm.csv", np.array([bjd, w4mag, w4mage]).T)
            except: None

    def Plot():
        try:
            mjd,mag,mage=np.loadtxt("WISE1_340nm.csv",unpack=True)
            plt.scatter(mjd,mag, c='b', label="WISE1")
            plt.errorbar(mjd,mag,yerr=mage,ls=" ", c='b')
        except: None
        
        try:
            mjd,mag,mage=np.loadtxt("WISE2_460nm.csv",unpack=True)
            plt.scatter(mjd,mag, c='g', label="WISE2")
            plt.errorbar(mjd,mag,yerr=mage,ls=" ", c='g')
        except: None
        
        try:
            mjd,mag,mage=np.loadtxt("WISE3_1200nm.csv",unpack=True)
            plt.scatter(mjd,mag, c='orange', label="WISE3")
            plt.errorbar(mjd,mag,yerr=mage,ls=" ", c='orange')
        except: None
        
        try:
            mjd,mag,mage=np.loadtxt("WISE4_2200nm.csv",unpack=True)
            plt.scatter(mjd,mag, c='brown', label="WISE4")
            plt.errorbar(mjd,mag,yerr=mage,ls=" ", c='brown')
        except: None
        
        



