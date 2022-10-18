#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 10:55:09 2022

@author: james
"""

import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import EarthLocation
from miscAstro import * # analysis:ignore
import matplotlib.pyplot as plt


# ASASSN can do 8-18 mag objects

class ASASSN(object):
    def compressASASSN():
        t = Table.read('/Users/james/Desktop/ASASSN/file_list/asassn_catalog_full.csv')
        
        print(t.colnames)
        
        idtab=t['id']
        source_idtab=t['source_id']
        asassn_name=t['asassn_name']
        edr3_id=t['edr3_source_id']
        galex_id=t['galex_id']
        mean_vmag=t['mean_vmag']
        apass_dr9_id=t['apass_dr9_id']
        
        ra=t['raj2000']
        dec=t['dej2000']
        
        
        tbhdu = fits.BinTableHDU.from_columns(
             [fits.Column(name='id', format='37A', array=idtab),
              fits.Column(name='source_id', format='7A', array=source_idtab),
              fits.Column(name='asassn_name', format='29A', array=asassn_name),
              fits.Column(name='edr3_source_id', format='25A', array=edr3_id),
              fits.Column(name='galex_id', format='23A', array=galex_id),
              fits.Column(name='mean_vmag', format='E', array=mean_vmag),
              fits.Column(name='apass_dr9_id', format='E', array=apass_dr9_id),
              fits.Column(name='raj2000', format='E', array=ra),
              fits.Column(name='dej2000', format='E', array=dec)])
        tbhdu.writeto('compressed_asassn_catalog_full.fits', overwrite=True)
        # sometimes edr3_source_id and galex_id are N/A
    
    
    
    def isItInASASSN(GaiaEDR3id, RAdeg, Decdeg, rad):#rad in deg
        t = Table.read('/Users/james/Desktop/ASASSN/file_list/compressed_asassn_catalog_full.fits')
        
        #GaiaEDR3id="a" # set this to a random string if you want to check if you want to check what happens when there is no gaia crossmatch
        #read gaia source id here
        arg=np.argwhere(t['edr3_source_id']==GaiaEDR3id)
        if len(arg)>0:
            ASASSN_name=t['asassn_name'][arg][0]
            ASASSN_name=np.char.decode(ASASSN_name)
        else:
            ra=t['raj2000'].value # in degrees!
            dec=t['dej2000'].value # in degrees!
            if Decdeg>=0:
                mask = ((ra < RAdeg+rad) & (ra>RAdeg-rad) & (dec>Decdeg-rad) & (dec<Decdeg+rad))
            else:
                mask = ((ra < RAdeg+rad) & (ra>RAdeg-rad) & (np.abs(dec)>abs(Decdeg)-rad) & (np.abs(dec)<abs(Decdeg)+rad) & (dec < 0))
            ASASSN_name=t['asassn_name'][mask]
        
        
        if len(ASASSN_name)==1:
            return str(ASASSN_name[0])
        else:
            try:
                ratar=ra[mask]
                dectar=dec[mask]
                absval=np.sqrt((ratar-RAdeg)**2 + (dectar-Decdeg)**2)
                
                closestASASSN_name=ASASSN_name[np.argmin(absval)]
                
                print("NO FOUND GAIAEDR3-ASASSN CROSS-MATCH AND MULTIPLE STARS IN SEARCH RADIUS: ASSUMING THAT THE CLOSEST STAR IS THE DESIRED (WITHIN THE SEARCH RADIUS", np.round(rad,3)," DEGREES)")
                np.savetxt("ASASSN_warning.dat", np.array([np.full((len(ratar),),GaiaEDR3id), np.full((len(ratar),),RAdeg), np.full((len(ratar),),Decdeg),ratar,dectar,np.full((len(ratar),),closestASASSN_name),np.full((len(ratar),),rad)]).T, fmt="%s", header="attempted_GaiaEDR3id,  input_RA,  input_Dec,  RAs_tar_ASASSN,  Decs_tar_ASASSN,  closest_ASASSN, search_radius")
                return closestASASSN_name
            except Exception as e: # can't do arithmetic on a value that doesn't exist, will error if ratar is empty
                print(e)    
                return "NO ASASSN"
    
    
    
    def getASASSN(filename):
        # ASASSN-V_J161010.64+053305.0.dat       g_band_lcs format
        # ASASSN-VJ235959.33+532742.6.dat        vardb_files format
        
        # first get g band
        filename12=filename.split(" ")
        filenameG=filename12[0]+"_"+filename12[1]+".dat"
        
        # then get v band
        filenameV=filename12[0]+filename12[1]+".dat"
        
        return filenameG, filenameV
     
     
    
    def getGband(RAdeg, Decdeg,filename):
        #https://iopscience.iop.org/article/10.1088/1538-3873/aa80d9/pdf
        Hawaii_Haleakala_Observatory = EarthLocation.from_geodetic(lat=20.7082,lon=-156.2568,height=3052)
        McDonald = EarthLocation.from_geodetic(lat=30.6797,lon=-104.0247,height=2077)
        SAAO=EarthLocation.from_geodetic(lat=-33.9345,lon=18.4771,height=2077)
        Cerro_Tololo_International_Observatory=EarthLocation.from_geodetic(lat=-30.169661,lon=-70.806525,height=2207)
        #ba, bb, bc, bd, = Hawaii_Haleakala_Observatory
        #bi,bj,bk,bl = McDonald
        #bm,bn,bo,bp=SAAO
        #be, bf, bg, bh, bq,br,bs,bt= Cerro_Tololo_International_Observatory
         
        
        HJD, camera, mag, mag_err=[],[],[],[]#, flux, flux_err, FWHM, IMAGE=[],[],[],[],[],[],[],[]
        tel=[]
        with open("/Users/james/Desktop/ASASSN/g_band_lcs/"+str(filename)) as file:
            for cnt, line in enumerate(file):
                if cnt!=0:
                    # the float below are to catch  if there is any dodgy entry. sometimes if there is a dodgy measurement, it looks like this:
                    # 2458393.50283	bH	>18.133	99.99000	-0.07700	0.04100	1.94000	coadd_bH001836
                    # ValueError: could not convert string to float: '>18.13'
                    # I ignore these completely
                    # a good example:
                    # ASASSN/g_band_lcs/ASASSN-V_J161008.57-561533.6.dat
                    try:
                        intermediate_camera=line[14:16]
                        
                        float(line[:13])
                        float(line[17:23])
                        float(line[24:31])
                        #float(line[32:39])
                        #float(line[40:47])
                        
                        HJD.append(float(line[:13]) )
                        camera.append(line[14:16])
                        mag.append(float(line[17:23]))
                        mag_err.append(float(line[24:31]))
                        # check these again if you need to use: doesn't account for any negative flux sign etc
                        #flux.append(float(line[32:39]))
                        #flux_err.append(float(line[40:47]))
                        #FWHM.append(float(line[48:55]))
                        #IMAGE.append(line[56:70])
                        
                        if intermediate_camera.lower()=="ba" or intermediate_camera.lower()=="bb" or intermediate_camera.lower()=="bc" or intermediate_camera.lower()=="bd":
                            tel.append([Hawaii_Haleakala_Observatory.lon.value,Hawaii_Haleakala_Observatory.lat.value,Hawaii_Haleakala_Observatory.height.value])
                        elif intermediate_camera.lower()=="bi" or intermediate_camera.lower()=="bj" or intermediate_camera.lower()=="bk" or intermediate_camera.lower()=="bl":
                            tel.append([McDonald.lon.value,McDonald.lat.value,McDonald.height.value])
                        elif intermediate_camera.lower()=="bm" or intermediate_camera.lower()=="bn" or intermediate_camera.lower()=="bo" or intermediate_camera.lower()=="bp" :
                            tel.append([SAAO.lon.value,SAAO.lat.value,SAAO.height.value])
                        elif intermediate_camera.lower()=="be" or intermediate_camera.lower()=="bf" or intermediate_camera.lower()=="bg" or intermediate_camera.lower()=="bh" or intermediate_camera.lower()=="bq" or intermediate_camera.lower()=="br" or intermediate_camera.lower()=="bs" or intermediate_camera.lower()=="bt":
                            tel.append([Cerro_Tololo_International_Observatory.lon.value,Cerro_Tololo_International_Observatory.lat.value,Cerro_Tololo_International_Observatory.height.value])
                        else: print(intermediate_camera)
                        
                    except: None
        
        tel=np.asarray(tel).T
        #flux=np.asarray(flux); flux_err=np.asarray(flux_err)
        #FWHM=np.asarray(FWHM); IMAGE=np.asarray(IMAGE)
        
        BJD=miscAstro.hjd_to_bjd(RAdeg,Decdeg,HJD,tel)
        
        return BJD-2400000.5, np.asarray(mag), np.asarray(mag_err)
    
    
    
    
    def getVband(RAdeg, Decdeg, filenameV):
        HJD, mag, mag_err, flux, flux_err= np.loadtxt("/Users/james/Desktop/ASASSN/vardb_files/"+str(filenameV), unpack=True, skiprows=2)
        
        Hawaii_Haleakala_Observatory = EarthLocation.from_geodetic(lat=20.7082,lon=-156.2568,height=3052)
        
        # this isn't really necessary, but included incase you want to separate site with extra knowledge OR according to declination (if it's -90 dec, safe to say it should be chile)
        tel=np.tile([Hawaii_Haleakala_Observatory.lon.value,Hawaii_Haleakala_Observatory.lat.value,Hawaii_Haleakala_Observatory.height.value],(len(HJD),1)).T
        
        BJD=miscAstro.hjd_to_bjd(RAdeg,Decdeg,HJD,tel)
    
        return BJD-2400000.5, mag, mag_err
    
    
    
    def plotVband(RAdeg,Decdeg,filenameV):
        BJD_V, mag_V, mag_err_V = ASASSN.getVband(RAdeg,Decdeg, filenameV)
        plt.clf()
        plt.title("Vband + " + filenameV)
        plt.scatter(BJD_V,mag_V)
        plt.errorbar(BJD_V,mag_V, yerr=mag_err_V,ls=" ")
        plt.gca().invert_yaxis()
        plt.savefig("ASASSNv.png")
        plt.clf()
        plt.close()
        np.savetxt("ASASSNv_lc.dat", np.array([BJD_V, mag_V, mag_err_V]).T)
    
    def plotGband(RAdeg,Decdeg,filenameG):
        BJD_G, mag_G, mag_err_G = ASASSN.getGband(RAdeg,Decdeg,filenameG)
        
        plt.clf()
        plt.title("Gband " + filenameG)
        plt.scatter(BJD_G,mag_G)
        plt.errorbar(BJD_G,mag_G, yerr=mag_err_G,ls=" ")
        plt.gca().invert_yaxis()
        plt.savefig("ASASSNg.png")
        plt.clf()
        plt.close()
        np.savetxt("ASASSNg_lc.dat", np.array([BJD_G, mag_G, mag_err_G]).T)
    
    

    def __init__():
        print("ASASSN EXPOSURES ARE 90s")
