#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 17:16:54 2022

@author: james
"""

import numpy as np
import csv, requests, os, datetime, jdcal
import matplotlib.pyplot as plt


class ZTF(object):
    # create request to ZTF UNFINISHED
    def create_url(ID, CIRCLE,BAND=False,BANDNAME=False,MAG=False,NUM_OBS=False,TIME=False,BAD_CATFLAGS_MASK=False,COLLECTION=False,FORMAT=False): #POS,
        # unfinished business here
        url="https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?"
        if ID:
            url+="ID="
            url+=str(ID)
            url+="&"
        
        if CIRCLE:
            RA, Dec, radius=CIRCLE[0], CIRCLE[1], CIRCLE[2] # radius in deg
            url+="POS=CIRCLE "
            url+=str(RA)
            url+=" "
            url+=str(Dec)
            url+=" "
            url+=str(radius)
            url+="&"
        
        
        #if POS:
        #    POS_ID=str(POS)
        
        
        if BAND: #in nanometres
            minwl,maxwl=BAND[0]/100,BAND[1]/100
            url+=str(minwl)
            url+="e-7 "
            url+=str(maxwl)
            url+="e-7&"
        
        
        if BANDNAME:
            url+="BANDNAME="
            url+=str(BANDNAME)
            url+="&"
        
        
        if MAG:
            minmag,maxmag=MAG[0],MAG[1]
            url+="MAG="
            url+=str(minmag)
            url+=" "
            url+=str(maxmag)
            url+="&"
        
        
        if NUM_OBS:
            url+="NUM_OBS="
            url+=str(NUM_OBS)
            url+="&"
        
        if TIME:
            "enter TIME   as  [minT maxT] or as [0 2y]"
            minT=TIME[0]
            maxT=TIME[1]
            if "last" in TIME:
                if "y" in TIME:
                    yplace=TIME.find("y")
                    lastplace=TIME.find("last")
                    numyears=TIME[3:yplace]
            
                    dt=datetime.datetime.now()
                    mjd=sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)) -  2400000.5
                    
                    MJDmax=mjd
                    MJDmin=mjd-float(numyears)*365
                else:
                    raise ValueError
                    
            else:
                MJDmax=minT
                MJDmin=maxT
            url+="TIME="
            url+=str(MJDmin)
            url+=" "
            url+=str(MJDmax)
            url+="&"
            
            
        
        
        if BAD_CATFLAGS_MASK:
            url+="BAD_CATFLAGS_MASK=32768&"
        
        if COLLECTION:
            None
        
        return url+"FORMAT=CSV"
    
    
    
        
    #grab the data
    
    def get_data(url, auth):
        response = requests.get(url,auth=auth)
        return response
    
    
    #save the data to new directory
    
    def save_data(count, raw_data):
        
        filename="ZTF.csv"
        with open(filename, 'w') as f:
            writer = csv.writer(f)
            for line in raw_data.iter_lines():
                writer.writerow(line.decode('utf-8').split(','))
        
        
        
        ZTF.plot_photometry(os.getcwd(), filename)
    
    
    def plot_photometry(cwd, filename, saveit=True):
        data=np.genfromtxt(filename, skip_header=1, unpack=True, encoding = None)#dtype=None)
        
    
        with open("ZTF.csv",'r') as dest_f:
            data_iter = csv.reader(dest_f,
                                   quotechar = '"')
            next(data_iter)
            data = [data for data in data_iter]
        
        oid, expid, hjd, mjd, mag, magerr, catflags, filtercode, ra, dec, chi, sharp, filefracday, field, ccdid, qid, limitmag, magzp, magzprms, clrcoeff, clrcounc, exptime, airmass, programid= np.asarray(data).T
        
        if saveit==True:
            plt.clf()
        
        for filt in np.unique(filtercode):
            
            mask=((filtercode==filt) & (((catflags=="0")) | (catflags=="2") | (catflags=="16") | (airmass.astype(np.float64) < 2.5) | (magzprms.astype(np.float64)<0.04)))
            
            # remove any flagged data (page 79) https://web.ipac.caltech.edu/staff/fmasci/ztf/ztf_pipelines_deliverables.pdf 
            
            plt.scatter(mjd.astype(float)[mask], mag.astype(float)[mask], label=filt)
            plt.errorbar(mjd.astype(float)[mask], mag.astype(float)[mask], yerr=magerr.astype(float)[mask], ls=' ')
            if saveit==True:
                np.savetxt("ZTF_"+str(filt)+".csv", np.array([mjd.astype(float)[mask], mag.astype(float)[mask], magerr.astype(float)[mask] ]).T, fmt="%s")
        if saveit==True:
            plt.legend(loc="upper right")
            plt.xlabel("MJD")
            plt.ylabel("mag")
            plt.gca().invert_yaxis()
            plt.autoscale()
            plt.savefig("ZTFphot.pdf")
            plt.close()





