#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 14 01:08:29 2022

@author: james
"""

from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
import requests
from bs4 import BeautifulSoup
from PIL import Image
#https://www.sdss.org/dr12/tutorials/quicklook/

class SDSSclass(object):
    def search_SDSS_phot(query):
        #note: can add colour/mag limits to this query
        res = SDSS.query_sql(query)
        return res
    
    
    def plot_SDSS_spec(name,spec):
        hdulist=spec[0]
        c0 = hdulist[0].header['coeff0']
        c1 = hdulist[0].header['coeff1']
        npix = hdulist[1].header['naxis2']
        wave = 10.**(c0 + c1 * np.arange(npix))
        
        bunit = hdulist[0].header['bunit']
        data = hdulist[1].data
        flux=data["flux"]
        ivar = data["ivar"]
        data2=hdulist[3].data
        line=data2["LINENAME"]
        linewave=data2["LINEWAVE"]
        
        SPECra=hdulist[0].header['PLUG_RA']
        SPECdec=hdulist[0].header['PLUG_DEC']
        
        justH=True
        plt.plot(wave, (flux) , 'k')
        for count, i in enumerate(line):
            if linewave[count] < np.amax(wave) and linewave[count] > np.amin(wave):
                if justH==True:
                    if "H_" in i:
                        plt.axvline(linewave[count], label="H", c='r', alpha=0.5)
                else:
                    if "H_" in i:
                        plt.axvline(linewave[count], label="H", c='r', alpha=0.5)
                    elif "He_" in i:
                        plt.axvline(linewave[count], label="He", c='b', alpha=0.5)
                    elif "O_" in i:
                        plt.axvline(linewave[count], label="O", c='g', alpha=0.5)
                    elif "N_" in i:
                        plt.axvline(linewave[count], label="N", c='y', alpha=0.5)
                    elif "S_" in i:
                        plt.axvline(linewave[count], label="S", c='orange', alpha=0.5)
                    elif "Ar_" in i:
                        plt.axvline(linewave[count], label="Ar", c='brown', alpha=0.5)
                    elif "Ne_" in i:
                        plt.axvline(linewave[count], label="Ne", c='brown', alpha=0.5)
                    else:
                        plt.axvline(linewave[count], label=str(i), c='red', alpha=0.5)
        
        
        handles, labels = plt.gca().get_legend_handles_labels()
        newLabels, newHandles = [], []
        for handle, label in zip(handles, labels):
          if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)
        plt.legend(newHandles, newLabels)
        
        plt.xlabel("Wavelength (AA)")
        try:
            plt.ylabel(bunit)
        except:
            plt.ylabel("Flux")
        np.savetxt("full_linelist.txt", np.asarray(line), fmt="%s")
        
        plt.xlim(np.amin(wave), np.amax(wave))
        plt.title("RA " + str(SPECra) + ", " + "Dec " + str(SPECdec))
        plt.savefig(str(name)+".png")
        hdulist.close()
        
    
    def get_SDSSmagsUGRIZ(name, RA,Dec):
        co = coords.SkyCoord(ra=RA*u.degree, dec=Dec*u.degree,frame="icrs")
        #result = SDSS.query_crossid(co, photoobj_fields=['modelMag_g', 'modelMag_i'])
        
        result = SDSS.query_crossid(co)
        
        magU=result["psfMag_u"].value[0]
        magUe=result["psfMagerr_u"].value[0]
        magG=result["psfMag_g"].value[0]
        magGe=result["psfMagerr_g"].value[0]
        magR=result["psfMag_r"].value[0]
        magRe=result["psfMagerr_r"].value[0]
        magI=result["psfMag_i"].value[0]
        magIe=result["psfMagerr_i"].value[0]
        magZ=result["psfMag_z"].value[0]
        magZe=result["psfMagerr_z"].value[0]
        
        
        mags=[magU,magG,magR,magI,magZ]
        magsE=[magUe,magGe,magRe,magIe,magZe]
        np.savetxt(str(name)+ ".csv", np.array([mags,magsE]).T)
        
        return magU,magG,magR,magI,magZ
    
    
    def search_SDSS_spectrum(query,RA, Dec):
        co = coords.SkyCoord(ra=RA*u.degree, dec=Dec*u.degree)
        result = SDSS.query_region(co, spectro=True)
        if result:
            spec = SDSS.get_spectra(matches=result, radius=1*u.arcsec)
        else:
            spec="NoSPEC"
        return spec
    
    def get_SDSSquery(RA,Dec,rad):
        #http://skyserver.sdss.org/dr12/en/tools/search/form/searchform.aspx
        query="select top 10 p.objid, p.ra, p.dec, p.u, p.g, p.r, p.i, p.z      from star p, dbo.fgetNearByObjEq("
        query+=str(RA)
        query+=","
        query+=str(Dec)
        query+=","
        query+=str(rad)
        query+=") n         where p.objid=n.objid"
        return query
        
    
    def get_findingchart_url(RA,Dec,scale):
        url="skyserver.sdss.org/dr16/en/tools/chart/image.aspx?ra="
        url+=str(RA)
        url+="&dec="
        url+=str(Dec)
        url+="&scale="
        url+=str(scale)
        url+="&height=512&width=512&opt=GO"
        
        return url
    
    
    def createShortcut(url, destination):
        text = '[InternetShortcut]\nURL=https://{}\nIconIndex=0'.format(url)
        with open(destination + 'SDSS.url', 'w') as fw:
            fw.write(text)
        return


    def scrapeImageFromShortcut(base_url):
        r = requests.get("http://"+base_url)
        soup = BeautifulSoup(r.text,'html.parser')
        for image_src in soup.find_all("img"):
            if "http://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/getjpeg?TaskName=Skyserver.Chart" in image_src['src']:
                #print(image_src['src'])
                try:
                    img = Image.open(requests.get(image_src['src'], stream = True).raw)
                    img.save('SDSS.png')
                except: None




        




