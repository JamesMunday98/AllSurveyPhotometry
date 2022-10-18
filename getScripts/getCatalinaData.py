#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 09:36:39 2022

@author: james
"""

import requests
import urllib
import numpy as np
import matplotlib.pyplot as plt
from miscAstro import *
from astropy.coordinates import EarthLocation
import astropy.time
import dateutil.parser


# "In cases where objects are within a 2-3 arcseconds, they will always be blended in CSS images. 
# Such cases can often be seen by comparison with higher resolution images from SDSS, or other 
# surveys."

class getCatalinaData(object):
    # IMPORTANT NOTE: if the search radius captures two stars, both stars are included in the CatalinaDataset.csv file
    # I will deal with this issue, but if you just take CatalinaDataset.csv then it will have contaminents
    
    def getData(RAdeg,Decdeg,ref_epoch=2020, pmra=0,pmdec=0, rad=3/60): # rad is in arcminutes! 
        
        ############## correct for proper motion (rough, but better than nothing)
        dt = dateutil.parser.parse(str(int(ref_epoch))+'.01.01')
        dt2=dateutil.parser.parse('2007.06.01') #Catalina was 2003-2012, use a guess that the mean location was that at 2007.5
        time = astropy.time.Time(dt)
        time2= astropy.time.Time(dt2)
        deltamjd_yr=(time.jd - time2.jd)/365
        
        
        # please ensure that dt > dt2 for this to work
        pos_ra_change = ((pmra/3600) /1000)  *(-1*deltamjd_yr)
        RAdeg+=pos_ra_change
        
        pos_dec_change = ((pmdec/3600) /1000)  *(-1*deltamjd_yr)
        Decdeg+=pos_dec_change
        ##############
    
    
    
    
        # if you make rad large, this script will not work. Instead, the html returned says about how
        # it is a large request and may take some time. doing a 2arcmin request is useless for me, 
        # so i don't include this case. if it is needed for you, feel free to develop this simple code
        # and let me know so I can also include it!
    
        url='http://nunuku.caltech.edu/cgi-bin/getcssconedb_release_img.cgi/post'
        # please note: to get the cgi file, I installed software called wireshark to track http activity
        # for readbility, you can add a display filter to show html only (e.g. html.request)
        # perform the action of manually getting the photometry from the webform at (url)
        # then I clicked the stop capturing packets button
        # now, go file -> export objects -> html
        # save these somewhere
        # find the file you want from the saved items
        # pass this below as the input format to the server at (url)
        
        files={'files': open('../../getScripts/getCatalinaCGI.cgi','rb')}
        # the keywords are simply the keywords inside the .cgi file, e.g. name="OUT"
        values={'upload_file' : 'getCatalinaCGI.cgi' , 'DB':'photcat' , 'OUT':'csv' , 'SHORT':'short', 'RA':str(RAdeg), 'Dec':str(Decdeg), 'Rad':str(rad)}
        r=requests.post(url,files=files,data=values)
    
        #print(r.text)
        a=r.text.split("web_file")
        
        
        try:
            link=a[0][-52:]+"web_file"+a[2][:-40]; #print(link)
            
            with urllib.request.urlopen(link) as testfile, open('CatalinaDataset.csv', 'w') as f:
                f.write(testfile.read().decode())
        except : None
            #np.savetxt("CATALINA_NOT_FOUND_ON_QUERY.txt", np.array(["this is weird because it was found in a vizier search"]), fmt="%s")
            
        
        try: files.close()
        except: None
        
        
        
    def plot(savefig=True):
        fi= np.genfromtxt('CatalinaDataset.csv', delimiter=',', unpack=True, skip_header=1)
        
        # mask to avoid any blended items
        MasterID, Mag, Magerr, RA, Dec, MJD, Blend =fi.T[np.array(fi[6]==0)].T
        
        

        for i in np.unique(MasterID):
            plt.title("if this is multicoloured, you are seeing multiple sources!!!")
            plt.scatter(MJD[MasterID==i],Mag[MasterID==i])
            plt.errorbar(MJD[MasterID==i],Mag[MasterID==i], yerr=Magerr[MasterID==i], ls=" ")
        
        if len(np.unique(MasterID)) > 1:
            np.savetxt("CATALINA_WARNING.txt", np.array(["two masterid entries are found here"]))
        
        
    
        plt.gca().invert_yaxis()
        plt.xlabel("BJD")
        plt.ylabel("Magnitude")
        plt.savefig("Catalina.pdf")
        plt.close()
        plt.clf()
            
    
    
    def handleData(RAdeg,Decdeg):
        fi= np.genfromtxt('CatalinaDataset.csv', delimiter=',', unpack=True, skip_header=1)
        
        # mask to avoid any blended items
        MasterID, Mag, Magerr, RA, Dec, MJD, Blend =fi.T[np.array(fi[6]==0)].T
        
        
        CatalinaLoc=EarthLocation.from_geodetic(lat=32.417, lon=-110.733, height=2510.028)
        
        
        BJD = miscAstro.jd_corr(MJD, RAdeg, Decdeg, CatalinaLoc, jd_type='bjd').value
        np.savetxt('CatalinaDataBJD.csv', np.array([BJD,Mag,Magerr]).T)
        
        



# do this in the future: http://nesssi.cacr.caltech.edu/DataRelease/FAQ.html#improve
# to get better mags. this sufficices for now though with difference mag being the only interest



