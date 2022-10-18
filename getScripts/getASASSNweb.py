#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 16:09:43 2022

@author: james
"""

import numpy as np
import matplotlib.pyplot as plt
import requests
from bs4 import BeautifulSoup
import json
import csv
from astropy.coordinates import EarthLocation
from miscAstro import *
import os
import pandas as pd
# this is the structure. on the site it says "good to 200s"!!!

class getASASSNweb(object):
    def get_urls(RAdeg,Decdeg,rad=4/60):
        base_url="https://asas-sn.osu.edu/photometry?utf8=%E2%9C%93&ra="+str(RAdeg)+"&dec="+str(Decdeg)+"&radius="+str(rad)+"&vmag_min=&vmag_max=&epochs_min=&epochs_max=&rms_min=&rms_max=&sort_by=raj2000"
        
        r = requests.get(base_url)
        soup = BeautifulSoup(r.text,'html.parser')
        for a in soup.find_all("a", href=True):
            #print(a)
            if "photometry.json?" in a['href']:
                jsonlink=a['href']
            if "photometry/" in a['href']:
                link_to_lightcurve=a['href']
        try:
            if jsonlink:
                None
        except: jsonlink=False
        try:
            if link_to_lightcurve:
                None
        except: link_to_lightcurve=False
        
        
        return "https://asas-sn.osu.edu/"+jsonlink, "https://asas-sn.osu.edu/"+link_to_lightcurve
        
    def get_json_info(json_link, lightcurve_link):
        req = requests.get(json_link).text
        
        d = json.loads(req)
        
        if len(d['results'])==1:
            res1 = d['results'][0]
            
            if res1['blend']==False:
                epochs=res1['epochs']
                ra=res1['raj2000']
                dec=res1['dej2000']
                np.savetxt("ASASSNweb_info.txt", np.array([epochs,ra,dec, lightcurve_link]).T, header="num epochs, raASASSN, decASASSN", fmt="%s")
        else: np.savetxt("ASSASNweb_empty.txt")
    def get_web_csv(link_to_lightcurve):
        
        req_new = requests.get(link_to_lightcurve+".csv")
        
        with open('ASASSNweb.csv', 'w') as f:
            writer = csv.writer(f)
            for line in req_new.iter_lines():
                writer.writerow(line.decode('utf-8').split(','))
        
        
        hjd, camera, mag, magerr, flux, fluxerr = np.loadtxt("ASASSNweb.csv", unpack=True)
        plt.scatter(hjd, mag)
        plt.errorbar(hjd,mag,yerr=magerr,ls=' ')
        plt.savefig("ASASNweb_notLS.png")
    
    
    
    def ASASSN_apphot_url(RAdeg, Decdeg):
        url="https://asas-sn.osu.edu/?dec="
        url+=str(Decdeg)
        url+="&ra="
        url+=str(RAdeg)
        return url

    
    def createShortcut(url, destination):
        text = '[InternetShortcut]\nURL={}\nIconIndex=0'.format(url)
        with open(destination + 'ASASSN_apphot.url', 'w') as fw:
            fw.write(text)
        return
    
    
    def processDownloadedWebFile(RAdeg, Decdeg):
        #https://iopscience.iop.org/article/10.1088/1538-3873/aa80d9/pdf
        Hawaii_Haleakala_Observatory = EarthLocation.from_geodetic(lat=20.7082,lon=-156.2568,height=3052)
        McDonald = EarthLocation.from_geodetic(lat=30.6797,lon=-104.0247,height=2077)
        SAAO=EarthLocation.from_geodetic(lat=-33.9345,lon=18.4771,height=2077)
        Cerro_Tololo_International_Observatory=EarthLocation.from_geodetic(lat=-30.169661,lon=-70.806525,height=2207)
        #ba, bb, bc, bd, = Hawaii_Haleakala_Observatory
        #bi,bj,bk,bl = McDonald
        #bm,bn,bo,bp=SAAO
        #be, bf, bg, bh, bq,br,bs,bt= Cerro_Tololo_International_Observatory
        
        
        for file in os.listdir(os.getcwd()):
            if file.startswith("light_curve_"):
        
                #HJD, UTdate, camera, FWHM, Limit, mag, mag_err, flux, flux_err, filt = np.genfromtxt(file, unpack=True, skip_header=1, delimiter=",", dtype=None)
# =============================================================================
#                 HJD=HJD.astype(np.float64)
#                 FWHM=FWHM.astype(np.float64)
#                 Limit=Limit.astype(np.float64)
#                 mag=mag.astype(np.float64)
#                 mag_err=mag_err.astype(np.float64)
#                 flux=flux.astype(np.float64)
#                 flux_err=flux_err.astype(np.float64)
# =============================================================================
                df = pd.read_csv(file)
                HJD=df['HJD'].to_numpy()
                camera=df['Camera'].to_numpy()
                mag=df['mag'].to_numpy()
                mag_err=df['mag_err'].to_numpy()
                filt=df['Filter'].to_numpy()
                
                
                
                m = mag_err==99.99
                HJD=HJD[~m]
                camera=camera[~m]
                mag=mag[~m].astype(np.float64)
                mag_err=mag_err[~m]
                filt=filt[~m]
                
                for count,i in enumerate(filt):
                    filt[count]="ASweb"+filt[count]
                
                
                tel=[]
                for c in camera:
                    if c.lower()=="ba" or c.lower()=="bb" or c.lower()=="bc" or c.lower()=="bd":
                        tel.append([Hawaii_Haleakala_Observatory.lon.value,Hawaii_Haleakala_Observatory.lat.value,Hawaii_Haleakala_Observatory.height.value])
                    elif c.lower()=="bi" or c.lower()=="bj" or c.lower()=="bk" or c.lower()=="bl":
                        tel.append([McDonald.lon.value,McDonald.lat.value,McDonald.height.value])
                    elif c.lower()=="bm" or c.lower()=="bn" or c.lower()=="bo" or c.lower()=="bp" :
                        tel.append([SAAO.lon.value,SAAO.lat.value,SAAO.height.value])
                    elif c.lower()=="be" or c.lower()=="bf" or c.lower()=="bg" or c.lower()=="bh" or c.lower()=="bq" or c.lower()=="br" or c.lower()=="bs" or c.lower()=="bt":
                        tel.append([Cerro_Tololo_International_Observatory.lon.value,Cerro_Tololo_International_Observatory.lat.value,Cerro_Tololo_International_Observatory.height.value])
                    else: print(c)
                
                
                
                tel=np.asarray(tel).T
                #flux=np.asarray(flux); flux_err=np.asarray(flux_err)
                #FWHM=np.asarray(FWHM); IMAGE=np.asarray(IMAGE)
                if RAdeg!=0 and Decdeg!=0: # I use this as a caveat to get plotting work
                    BJD=miscAstro.hjd_to_bjd(RAdeg,Decdeg,HJD,tel)
                else:
                    BJD=HJD
                
                return BJD-2400000.5, np.asarray(mag), np.asarray(mag_err), filt
                break
    
    
    
