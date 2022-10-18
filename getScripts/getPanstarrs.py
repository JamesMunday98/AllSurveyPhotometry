#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 28 16:14:23 2022

@author: james
"""

import urllib.request 
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import EarthLocation
from miscAstro import * # analysis:ignore
import os

# saturates at magnitude 12-14
# 30 to 60s exposure times
class getPanstarrs(object):
    def getAllData(RAdeg,Decdeg, rad=1/3600):
         # this is in degrees
        
        url="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/detection.csv?ra="+str(RAdeg)+"&dec="+str(Decdeg)+"&radius="+str(rad)+"&nDetections.gte=10"
        
        #https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/detection.votable?ra=210.802429&dec=54.348750&radius=0.0083333
        urllib.request.urlretrieve(url, "Panstars.csv")


        urlmean="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/mean.csv?ra="+str(RAdeg)+"&dec="+str(Decdeg)+"&radius="+str(rad)+"&nDetections.gte=10"
        urllib.request.urlretrieve(urlmean, "Panstarsmean.csv")
        
    def getPanstarrsMeanMags():
        s=Table.read('Panstarsmean.csv')
        if len(s)>0:
            gMeanMag=s['gMeanPSFMag']
            rMeanMag=s['rMeanPSFMag']
            iMeanMag=s['iMeanPSFMag']
            zMeanMag=s['zMeanPSFMag']
            yMeanMag=s['yMeanPSFMag']
            
            gMeanMagErr=s['gMeanPSFMagErr']
            rMeanMagErr=s['rMeanPSFMagErr']
            iMeanMagErr=s['iMeanPSFMagErr']
            zMeanMagErr=s['zMeanPSFMagErr']
            yMeanMagErr=s['yMeanPSFMagErr']
            
            
            mags=[gMeanMag.data[0],rMeanMag.data[0],iMeanMag.data[0],zMeanMag.data[0],yMeanMag.data[0]]
            magerrs=[gMeanMagErr.data[0],rMeanMagErr.data[0],iMeanMagErr.data[0],zMeanMagErr.data[0],yMeanMagErr.data[0]]
            
            np.savetxt("PanstarrsMeanMagsGRIZY.dat", np.array([mags,magerrs]).T)
        else:
            os.remove('Panstarsmean.csv')


    def getPanstarrsLCs(RAdeg,Decdeg):
        t = Table.read('Panstars.csv')
        if len(t)>0:
            t=t[((t['psfFlux']>0))]
            t["obsTime"]-=37/86400 # put this into UTC
            
            Hawaii_Haleakala_Observatory= EarthLocation.from_geodetic(lat=20.7082,lon=-156.2568,height=3052)
            BJD=miscAstro.jd_corr(t["obsTime"], RAdeg, Decdeg, Hawaii_Haleakala_Observatory).value
            
            
            
            
            AB=2.5*(23-np.log10(t['psfFlux']))-48.6
            dAB=np.abs(-2.5/np.log(10) * 1/t['psfFlux'] * t['psfFluxErr'])
            
            
            for count, i in enumerate(AB):
                if int(t['filterID'][count])==1:
                    colour="g"
                elif int(t['filterID'][count])==2:
                    colour="r"
                elif t['filterID'][count]==3:
                    colour="orange"
                elif t['filterID'][count]==4:
                    colour="k"
                elif t['filterID'][count]==5:
                    colour="grey"
                
                plt.scatter(BJD[count],i,c=colour)
                plt.errorbar(BJD[count],i,yerr=dAB[count],c=colour)
            
            plt.gca().invert_yaxis()
            plt.xlabel("MJD")
            plt.ylabel("Magnitude")
            plt.savefig("Panstars.pdf")
            plt.clf()
            
            
            # What are the brightest and faintest stars for which the data are reliable?
            # The answer is it depends on band and FWHM. A very conservative estimate for the bright limit is as follows:
            # g r i z y
            # 14.5	15	15	14	13
            
            mask1=((t['filterID'] == 1) & (AB>14.5))
            mask2=((t['filterID'] == 2) & (AB>15))
            mask3=((t['filterID'] == 3) & (AB>15))
            mask4=((t['filterID'] == 4) & (AB>14))
            mask5=((t['filterID'] == 5) & (AB>13))
            
            if len(AB[mask1])>0:
                np.savetxt("PANSTARRS_g.dat", np.array([BJD[mask1], AB[mask1], dAB[mask1]]).T)
            if len(AB[mask2])>0:
                np.savetxt("PANSTARRS_r.dat", np.array([BJD[mask2], AB[mask2], dAB[mask2]]).T)
            if len(AB[mask3])>0:
                np.savetxt("PANSTARRS_i.dat", np.array([BJD[mask3], AB[mask3], dAB[mask3]]).T)
            if len(AB[mask4])>0:
                np.savetxt("PANSTARRS_z.dat", np.array([BJD[mask4], AB[mask4], dAB[mask4]]).T)
            if len(AB[mask5])>0:
                np.savetxt("PANSTARRS_y.dat", np.array([BJD[mask5], AB[mask5], dAB[mask5]]).T)
            
        else:
            os.remove('Panstars.csv')



    def plot():
        BJD, Panstarrsmag, PanstarrsmagErr = np.loadtxt("PANSTARRS_g.dat", unpack=True)
        
        plt.scatter(BJD,Panstarrsmag)
        plt.errorbar(BJD,Panstarrsmag,yerr=PanstarrsmagErr,ls=" ", label="gPSTRS")
        
        
        
        
        BJD, Panstarrsmag, PanstarrsmagErr = np.loadtxt("PANSTARRS_r.dat", unpack=True)
        
        plt.scatter(BJD,Panstarrsmag)
        plt.errorbar(BJD,Panstarrsmag,yerr=PanstarrsmagErr,ls=" ", label="rPSTRS")
        
        
        
        BJD, Panstarrsmag, PanstarrsmagErr = np.loadtxt("PANSTARRS_i.dat", unpack=True)
        
        plt.scatter(BJD,Panstarrsmag)
        plt.errorbar(BJD,Panstarrsmag,yerr=PanstarrsmagErr,ls=" ", label="iPSTRS")
        
        
        
        BJD, Panstarrsmag, PanstarrsmagErr = np.loadtxt("PANSTARRS_z.dat", unpack=True)
        
        plt.scatter(BJD,Panstarrsmag)
        plt.errorbar(BJD,Panstarrsmag,yerr=PanstarrsmagErr,ls=" ", label="zPSTRS")
        
        
        
        BJD, Panstarrsmag, PanstarrsmagErr = np.loadtxt("PANSTARRS_y.dat", unpack=True)
        
        plt.scatter(BJD,Panstarrsmag)
        plt.errorbar(BJD,Panstarrsmag,yerr=PanstarrsmagErr,ls=" ", label="yPSTRS")
        
        
        
