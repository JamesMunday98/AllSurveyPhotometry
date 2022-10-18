#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 11:16:28 2022

@author: james
"""

import numpy as np
import matplotlib.pyplot as plt
from miscAstro import *
from astropy.coordinates import EarthLocation
from astropy.timeseries import LombScargle
import astropy.units as u


#a = np.genfromtxt("forcedphotometry_req00177348_lc.txt",dtype=str)
#a = np.genfromtxt("forcedphotometry_req00015197_lc.txt",dtype=str)
#a = np.genfromtxt("forcedphotometry_req00188520_lc.txt",dtype=str) # abel 70
a = np.genfromtxt("forcedphotometry_v407_lc.txt",dtype=str) # abel 70

b=a[1:]

bmasknull=((b.T[20] != "null") & (b.T[24] !="null") & (b.T[25] != "null"))
b=b[bmasknull]


filt=b.T[4]
zpdiff = b.T[20].astype(float)
jd=b.T[22].astype(float)
forcediffimflux=b.T[24].astype(float)
forcediffimfluxunc=b.T[25].astype(float)
SNT=3
SNU=5



mask= (forcediffimflux / forcediffimfluxunc) > SNT
jd=jd[mask]
zpdiff=zpdiff[mask]
forcediffimflux=forcediffimflux[mask]
forcediffimfluxunc=forcediffimfluxunc[mask]
filt=filt[mask]



mag = zpdiff - 2.5*np.log10(forcediffimflux)
umag = 1.0857* forcediffimfluxunc / forcediffimflux



#ra=121.59563339553992
#dec=15.458611722079173
ra=307.8883201214792
dec=-7.088376336794571
ZTFloc=EarthLocation.from_geodetic(lat=33.3563, lon=-116.8650, height=1712)
bjd=miscAstro.jd_corr(jd-2400000.5, ra, dec, ZTFloc, jd_type='bjd').value



for i in ["ZTF_i", "ZTF_r", "ZTF_g"]:
    try:
        mask=filt==i
        if i=="ZTF_r":
            col='r'
        elif i=="ZTF_g":
            col='g'
        elif i=="ZTF_i":
            col='orange'
        #plt.errorbar(jd[mask],mag[mask],yerr=umag[mask],ls=' ', c=col)
        plt.scatter(bjd[mask], mag[mask], label=i, c=col)
        plt.errorbar(bjd[mask], mag[mask], yerr=umag[mask], ls=' ', c=col)
        
    except: None
    
plt.gca().invert_yaxis()
plt.legend()
plt.xlabel("BJD")
plt.ylabel("Mag")
plt.show()


plt.clf()

mask=filt=="ZTF_g" # tweak this
maskagain=((mag[mask] > 0) & (mag[mask] <27))



# =============================================================================
# fdot = 3.558E-16
# tminust0insec = (bjd-np.amin(jd-2400000.5))*86400
# 
# timediff=fdot*(tminust0insec**2)/2 * (1/0.00311)/86400
# 
# 
# bjd -= timediff
# =============================================================================







#frequency = np.linspace(2E-3, 4E-3, 1000000)*u.Hz
#frequency = np.linspace(0.4/86400, 0.8/86400, 1000000)*u.Hz
frequency = 1/np.linspace(200/86400, 10000/86400, 10000)*u.Hz
ls=LombScargle(bjd[mask][maskagain]*u.day, 
                               mag[mask][maskagain]*u.mag, umag[mask][maskagain]*u.mag)
power = ls.power(frequency)

plt.plot(1/frequency, power)

probabilities = [0.1, 0.05, 0.0001]

tenpercentFivepercentOnepercent=ls.false_alarm_level(probabilities)
for k in tenpercentFivepercentOnepercent:

    plt.axhline(k, c='r')

#plt.xlim(1/4E-3,1/2E-3)

#plt.axvline(177811, c='g', lw=2)
plt.xlabel("Period (s)")
plt.ylabel("Power")
plt.show()


print(len(frequency), frequency.min(), frequency.max()  )


iii = ls.false_alarm_probability(power)


plt.clf()
plt.plot(1/frequency,-np.log(iii))



for l in ls.false_alarm_probability(tenpercentFivepercentOnepercent):
    plt.axhline(-np.log(l), c='r')



#plt.xlim(1/4E-3,1/2E-3)
#plt.axvline(177811, c='g', lw=2)
plt.xlabel("Period (s)")
plt.ylabel("FAP")
plt.show()








