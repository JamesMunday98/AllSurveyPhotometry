#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 15 23:40:25 2022

@author: james
"""

# general libraries
import os, sys, sqlite3, matplotlib
sys.path.insert(0, os.getcwd()+"/getScripts")
import matplotlib.pyplot as plt
#from termcolor import colored
import numpy as np
from pathlib import Path
from astropy.table import Table
from multiprocessing import Process
from datetime import datetime

# classes used (mine)
from getK2 import K2
from getTESSnew import TESS
from getZTF import ZTF
from getATLASforcedPHOT import getATLASforcedPHOT
from miscAstro import * # analysis:ignore
from get_SDSS import SDSSclass
from get_photometric_SED_CDS import photometricSED
from OverplotGaiaLocation import OverplotGaia
from getPTF import PTF
from getPanstarrs import getPanstarrs
from plotEverything import plotEverything
from getASASSN import ASASSN
from getWISE import WISE
from getNEOWISE import NEOWISE
from getCatalinaData import getCatalinaData
from getCDS import CDS
from checkLocalStars import checkLocalStars
from getpwds import getpwds
from getASASSNweb import getASASSNweb
from distutils.dir_util import copy_tree
from getGaiaDatalink import GaiaDatalink


matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rc('font', size=14)
#plt.xticks(fontsize=12)
#plt.yticks(fontsize=12)

import warnings
warnings.filterwarnings("ignore")

#address warnings:
#import warnings
#warnings.filterwarnings('error')



def Final_SDSS(RAdeg, Decdeg):
    # get_SDSS
    if wantSDSS==True: # this is globally set at the start
        rad=1/60 #arcmin
        scale=0.1
        try:
            query=SDSSclass.get_SDSSquery(RAdeg,Decdeg,rad)
            photSDSS = SDSSclass.search_SDSS_phot(query)
            
            if photSDSS:
                spec=SDSSclass.search_SDSS_spectrum(query, RAdeg, Decdeg)
                if not spec=="NoSPEC":
                    SDSSclass.plot_SDSS_spec("spec",spec)
                
                
                uv,g,r,i,z=SDSSclass.get_SDSSmagsUGRIZ("ugriz",RAdeg, Decdeg)
                
                url=SDSSclass.get_findingchart_url(RAdeg,Decdeg,scale)
                SDSSclass.createShortcut(url, str(cwd)+"/")
                SDSSclass.scrapeImageFromShortcut(url)
                return uv,g,r,i,z
            else: 
                return "a", "a", "a", "a", "a"
        except: None


def Final_phot_SED_CDS(RADec):
    # get_photometric_SED_CDS
    search_radius_SED = 3 # arcsec
    try:
        if wantSED: # this is globally set at the start
            url=photometricSED.get_url_CDS(RADec, rad=search_radius_SED)
            photometricSED.plot_SED(url)
    except: None


def Final_K2(RADec):
    # getK2 
    search_radius_K2 = 4 # arcsec
    if wantK2==True: # this is globally set at the start
        split=RADec.split(":")
        split2=split[2].split(" ")
        # kepler coverage: RA=19h22m40s and Dec=+44˚30'00"(J2000)
        
        if int(split[0]) >=18 and int(split[0]) <=20 and int(split2[1]) >= 40 and int(split2[1])<=50: # this is where K2 imaged
            for time in exptimes:
                K2.get_K2(RADec, exptime=time, radius=4, ignore_any_dodgyness=False) # radius in arcsec
        else: None #print("Out of range of K2", RADec)



def Final_TESS(RADec, g):
    ########### getTESSnew  NOTE: made a limiting mag for tess of 18
    radius_TESS = 4 # arcseconds
    if wantTess==True:
        try: # break if g is "a", don't pass if g<18
            if g < 18: # this is globally set at the start
                for time in exptimes:
                    TESS.get_tess(RADec, time=time, radius=radius_TESS, ignore_any_dodgyness=False) # ignore_any_dodgyness= True # this command will not process any lightcurve that includes nans at any point
                # radius in arcsec
        except:# g is "a"
            for time in exptimes:
                TESS.get_tess(RADec, time=time, radius=radius_TESS, ignore_any_dodgyness=False) # ignore_any_dodgyness= True # this command will not process any lightcurve that includes nans at any point



def Final_ZTF(RAdeg, Decdeg, RA, Dec):
    #getPTF
    radius_ZTF = 2.5/3600   #  arcseconds in degrees
    if Decdeg > -34:
        if wantPTF==True: # this is globally set at the start
            PTF.queryPTF(RA,Dec)
            PTF.splitRandG_PTF()
        # getZTF
        if wantZTF==True: # this is globally set at the start
            urlZTF=ZTF.create_url(None, [RAdeg,Decdeg,radius_ZTF], BAD_CATFLAGS_MASK=True) # radius in deg
            try: ZTF.save_data(RADec,     ZTF.get_data(urlZTF, (getpwds.ZTF()[0], getpwds.ZTF()[1])))
            except: None

        

def Final_ATLAS_forced(RAdeg, Decdeg, RA, Dec, reference_epoch, pmra, pmdec): #come back to this and also don't recompute stuff I already have
    # get ATLAS forced photometry
    # you can change how you want the photometry to be computed - difference imaging or placing an aperture. check the online docs and the script that handles ATLAS things
    # https://fallingstar-data.com/forcedphot/apiguide/
    try:
        if wantATLASforced==True: # this is globally set at the start
        
            if "data.txt" in os.listdir(os.getcwd()) and not "dataATLASalreadyprocessed.txt" in os.listdir(os.getcwd()): #if I manually had to get the file from ATLAS because of a break
                a = dfresult = pd.read_csv("data.txt", delim_whitespace=True)
                if minimumMJD==50000:
                    getATLASforcedPHOT.plot_and_save_data(a,"FirstTime")
                else:
                    getATLASforcedPHOT.plot_and_save_data(a,"AddToOriginal")
                np.savetxt("dataATLASalreadyprocessed.txt", np.array(["nope"]))
        
            else:
                if "ATLAS_filtC_flux_and_err.dat" in os.listdir(os.getcwd()) or "ATLAS_filtO_flux_and_err.dat" in os.listdir(os.getcwd()):
                    try:
                        MJD_c,uJy_c, duJy_c, RA_c, Dec_c = np.loadtxt("ATLAS_filtC_flux_and_err.dat", unpack=True)
                        maxMJD_c=np.amax(MJD_c)+0.1 # +0.1 just to ignore the first measurement
                    except: None
                    
                    try:
                        MJD_o,uJy_o, duJy_o, RA_o, Dec_o = np.loadtxt("ATLAS_filtO_flux_and_err.dat", unpack=True)
                        maxMJD_o=np.amax(MJD_o)+0.1 # +0.1 just to ignore the first measurement
                    except: None
                    
                    # first get the largest MJD value between both files if they exist
                    try: minimumMJD=np.amax(np.array([maxMJD_c,maxMJD_o]))
                    except:
                        try: #otherwise get the largest of the _c file
                            minimumMJD=maxMJD_c*-1
                        except: #otherwise get the largest of the _o file
                            minimumMJD=maxMJD_o*-1        
                else: minimumMJD=50000 # else there is no prior entry and we ask for all the data
                
                dt = datetime.today().strftime('%Y-%m-%d')
                today=(list(sqlite3.connect(":memory:").execute("select julianday('" + dt + "')"))[0][0] -2400000.5)
                
                if abs(today-minimumMJD)>182.5:
                    a=getATLASforcedPHOT.ATLAS(getpwds.ATLAS()[0], getpwds.ATLAS()[1], RAdeg, Decdeg, 
                                               reference_epoch, pmra, pmdec, minMJD=minimumMJD)#minimumMJD)
                    
                    if minimumMJD==50000:
                        getATLASforcedPHOT.plot_and_save_data(a,"FirstTime")
                    else:
                        getATLASforcedPHOT.plot_and_save_data(a,"AddToOriginal")
                
    except:
        with open("../../list_of/bad_atlas.txt", "a") as atlasfile:
            atlasfile.write(str(RA) + " " + str(Dec)+ "\n")



def Final_Catalina(RAdeg,Decdeg,ref_epoch,pmra,pmdec):
    # get Catalina
    radius_catalina = 2.5/60  # passed as arcsec, so this is 3 arc seconds
    if wantCatalina==True and not "Catalina.csv" in os.listdir(os.getcwd()): # this is globally set at the start
        try:
            getCatalinaData.getData(RAdeg,Decdeg,ref_epoch,pmra,pmdec,radius_catalina)
            getCatalinaData.plot()
            getCatalinaData.handleData(RAdeg,Decdeg)
        except: None
    

def FinalGAIA(RADec, BP_RP, Abs_g, gmag):
    if wantGaiaHR ==True:
        OverplotGaia.plotGaia(RADec, BP_RP, Abs_g, gmag)

def FinalPanstarrs(RAdeg,Decdeg):
    radius_Panstarrs = 3/3600
    if wantPanstarrs==True :
        getPanstarrs.getAllData(RAdeg,Decdeg, rad=radius_Panstarrs)
        getPanstarrs.getPanstarrsLCs(RAdeg,Decdeg)
        getPanstarrs.getPanstarrsMeanMags()

def FinalASASSN(RAdeg, Decdeg, eDR3name="a"):
    radius_ASASSN = 5/3600
    if wantASASSN==True and not ("ASASSNv_lc.dat" in os.listdir(os.getcwd()) or "ASASSNg_lc.dat" in os.listdir(os.getcwd())):
        
        ASASSN_name = ASASSN.isItInASASSN(eDR3name, RAdeg, Decdeg, radius_ASASSN)
        filenameG, filenameV = ASASSN.getASASSN(ASASSN_name)
        
        try: ASASSN.plotVband(RAdeg,Decdeg,filenameV)
        except: 
            try: 
                try: json_link, lightcurve_link = getASASSNweb.get_urls(RAdeg, Decdeg)
                except: np.savetxt("ASSASNweb_empty.txt", np.array([1]))
                getASASSNweb.get_json_info(json_link, lightcurve_link)
                ASASSNweb.get_web_csv(lightcurve_link)
            except: None
        
        try: ASASSN.plotGband(RAdeg,Decdeg,filenameG)
        except: None # NOTE for the future : I don't know if the catalogue has any G only entries
        
        apphoturl=getASASSNweb.ASASSN_apphot_url(RAdeg,Decdeg)
        getASASSNweb.createShortcut(apphoturl, str(cwd)+"/")
        

def FinalWISE(RAdeg, Decdeg, gmag):
    if wantWISE == True and not "WISE.csv" in os.listdir(os.getcwd()) and gmag<16.6:
        try:
            WISE.queryWise(RAdeg,Decdeg,ref_epoch=2020,pmra=0,pmdec=0)
            WISE.saveFilters(RAdeg,Decdeg, gmag)
        except: None

def FinalNEOWISE(RAdeg, Decdeg, gmag):
    if wantNEOWISE == True and gmag<16.6:#and not "NEOWISE.csv" in os.listdir(os.getcwd()):
        try:
            NEOWISE.queryWise(RAdeg,Decdeg)
            NEOWISE.saveFilters(RAdeg,Decdeg, gmag)
        except Exception as e: print("error neowise"); print(e)

def FinalCDS(RAdeg, Decdeg):
    if wantCDS==True:
        url=CDS.CDSurl(RAdeg,Decdeg)
        CDS.createShortcut(url, str(cwd)+"/")






if __name__ == '__main__':
    # don't touch these two
    joinTESS=True
    joinK2=True
    
    # What data do you want?
    wantSDSS=True
    wantK2=True
    wantTess=True
    wantZTF=True
    wantATLASforced=True # this might take a couple of minutes as a request is queued to their server
    wantCatalina=True
    wantSED=True
    wantPTF=True
    wantPanstarrs=True
    wantGaiaHR=True
    wantWISE=True
    wantNEOWISE=True # this can be long to query
    wantASASSN=True
    wantCDS=True
    wantGaiaDatalink=True
    
    remove_old_dir=False  # are you sure? put twice to make sure you are positive. this will junk the contents of the created folders
    remove_old_dir=False # are you sure? put twice to make sure you are positive
    
    # TESS exposure times to search for? not all are always available
    exptimes=['fast', 'short', 'long']
    t = Table.read("example.fits")


    NameOfNewDir = "Objects"
    try: os.mkdir(NameOfNewDir)
    except: None

    # this can be named anything, but this place has to be a separate folder to dump files into
    os.chdir(NameOfNewDir)
    homedirectory=os.getcwd()
    
    
    for count, (RAdeg, Decdeg) in enumerate(zip(t['ra'].value,t['dec'].value)):
        # Many of these things rely off of Gaia metrics. These are a) to handle proper motion in search queries b) plot on the Gaia HR diagram
        # You can set these to different values if you just care about getting the data... but the options are included as they are typical things to plot
        # Generic values could be:
        # Gaia_Gmag = XX  # you need this number! I apply saturation limits and minimum magnitudes for e.g. TESS, WISE
        # ref_epoch = 2020
        # propermRA = 0
        # propermDec = 0
        # probWD = 0
        # BPRP = 0
        # GaiaSourceID = YY # you need this if you want to get Gaia datalink data
        # GaiaABS_G = 0
        # Teff = 0
        
        Gaia_Gmag = t['phot_g_mean_mag_corrected'].value[count]   #   the Gaia magnitude in the Gband
        ref_epoch = t['ref_epoch'].value[count]   #   the Gaia reference epoch (e.g. 2016)
        propermRA = t['pmra'].value[count]     # the Gaia proper motion in RA
        propermDec = t['pmdec'].value[count]   # the Gaia proper motion in Dec
        probWD = t['Pwd'].value[count]    # This is specific for white dwarfs and is just for the title of a graph. set this equal to 0 if you do not care about WDs. I use it as the probability of a white dwarf
        BPRP = t['BP_RP'].value[count]    # Gaia BP-RP
        GaiaSourceID = t['source_id'].value[count]    # Gaia source ID. Working as of DR3
        GaiaABS_G = t['M_G'].value[count]   # Absolute magnitude in Gaia G
        Teff = t['teff_H'].value[count]     # Teff of the object
        
        

        RA,Dec = miscAstro.ra_dec_deg_to_hr(RAdeg,Decdeg)
        
        if Decdeg<0:    RADec=str(RA) + " " + str(Dec)
        else:  RADec=str(RA) + " +" + str(Dec)
        
        folder = os.getcwd()+"/"+str(RADec)
        
        try:
            if remove_old_dir==True:     miscAstro.remDir(folder)
        except:None
        
        # make new dir and go to dir. if this breaks, dir already exists and we enter current dir
        try: os.mkdir(str(RADec))
        except: None
        os.chdir(folder)
        
        try: os.mkdir("files")
        except: None
        
        cwd=os.getcwd()
        
        # bring back all files I neatened
        if remove_old_dir==False:
            for filename2 in os.listdir(cwd+"/files/"):
                Path(cwd+"/files/"+filename2).rename(cwd+"/"+filename2)
        
        
        np.savetxt('TargetRADecDegrees.dat', np.array([RAdeg,Decdeg]).T)
        
        
        #print("We have ", str(len(t['ra'].value) - count), "out of ", str(len(t['ra'].value)), "remaining (progress = ", str(np.round(100*count/len(t['ra'].value),2)), "% )")
        #print(colored((RA,Dec),'cyan')); print(RAdeg, Decdeg)
        
        try: uv,g,r,i,z=Final_SDSS(RAdeg,Decdeg)
        except: None
        
        
        
        p0 = Process(target = FinalNEOWISE(RAdeg, Decdeg, gmag=Gaia_Gmag))
        p0.start()
        
        p1 = Process(target = Final_phot_SED_CDS(RADec))
        p1.start()
        try:
            if wantTess == True or wantK2 == True or wantATLASforced==True:
                obj, star_mag = checkLocalStars.find_star_in_gaia_edr3(RAdeg,Decdeg) #added in case there is a formatting mismatch, otherwise you could use Nicola's catalogue
        except: None
        
        if wantTess==True:
            returnClause = checkLocalStars.localTESS(obj,star_mag)
            if returnClause == "Good":
                p2 = Process(target = Final_TESS(RADec,Gaia_Gmag))
                p2.start()
                joinTESS=True
        
        if wantK2==True:
            returnClause = checkLocalStars.localK2(obj,star_mag)
            if returnClause == "Good":
                p3 = Process(target = Final_K2(RADec))
                p3.start()
                joinK2=True
        
        if Gaia_Gmag >12.5: # saturation limit
            # note that I put my own quality control cut on the ztf data by airmass and zeropoint rms
            p4 = Process(target = Final_ZTF(RAdeg,Decdeg, RA, Dec))
            p4.start()
        try:
            if Gaia_Gmag >12.5 and Decdeg>=-45: # saturation limit
                if wantATLASforced==True:
                    returnClause = checkLocalStars.localATLAS(obj,star_mag)
                    if returnClause == "Good":
                        p5 = Process(target = Final_ATLAS_forced(RAdeg,Decdeg,RA,Dec,reference_epoch=ref_epoch, pmra=propermRA, pmdec=propermDec))
                        #p5 = Process(target = Final_ATLAS_forced(RAdeg,Decdeg,RA,Dec,reference_epoch=2016, pmra=-146.303, pmdec=-155.864))
                        p5.start()
        except: None
        
        if Gaia_Gmag >= 13: # saturation limit
            p6 = Process(target = Final_Catalina(RAdeg,Decdeg, ref_epoch=ref_epoch, pmra=propermRA, pmdec=propermDec)) # bit larger since not as good astrometric solution
            p6.start()
        
        if Gaia_Gmag >= 12: # saturation limit
            p7 = Process(target = FinalPanstarrs(RAdeg, Decdeg))
            p7.start()
        
        #if Gaia_Gmag >= 11: # saturation limit
        #    p8 = Process(target = FinalASASSN(RAdeg, Decdeg, eDR3name="EDR3 "+str(GaiaSourceID)))
        #    p8.start()
        
        p9 = Process(target = FinalWISE(RAdeg, Decdeg, gmag=Gaia_Gmag))
        p9.start()
        
        # plot gaia hr
        try:
            p10 = Process(target = FinalGAIA(RADec, BPRP, GaiaABS_G, gmag=Gaia_Gmag))
            p10.start()
        except: None
        
        # get CDS shortcut
        p11 = Process(target = FinalCDS(RAdeg,Decdeg))
        p11.start()
        
        
        p0.join()
        p1.join() # these make it so that the bunch terminates when the final process pX does
        if wantTess==True and joinTESS==True:  
            try: p2.join()
            except: None
        if wantK2==True  and joinK2==True:  
            try: p3.join()
            except: None
        if Gaia_Gmag >12.5: # saturation limit
            try: p4.join()
            except: None
        if wantATLASforced==True:  
            try: p5.join()
            except: None
        if Gaia_Gmag >= 13: # saturation limit
            try: p6.join()
            except: None
        if Gaia_Gmag >= 12: # saturation limit
            try: p7.join()
            except: None
        #if Gaia_Gmag >= 11: # saturation limit
        #    try: p8.join()
        #    except: None
        p9.join()
        
        
        
        pX = Process(target=plotEverything.plot());       pX.start()
        
        if wantGaiaDatalink==True:
            try:
                source_id_input=str(GaiaSourceID)
                GaiaDatalink.getData(source_id_input)
            except: None
        
        
        
        
        # here is an example period search following the bare basics. you need to import the file you want.
        try:
            from astropy.timeseries import LombScargle
            from astropy import units as u
            MJD, mag, mage = np.loadtxt("ZTF_zg.csv", unpack=True)
            # remember that this is MJD!!! you should always convert your system to the barycentric reference frame.
            # see tdb Barycentric Dynamical Time (TDB)   under https://docs.astropy.org/en/stable/time/index.html
            # you can do it all with astropy routines
            frequency, power = LombScargle(MJD*u.d, mag, dy=mage).autopower()
            plt.plot(frequency, power)
            plt.xlabel("Frequency " +str(frequency.unit))
            plt.ylabel("Power")
            plt.savefig("ZTFperiodogram.png")
            plt.clf()
            
            
            # fold data at the highest peak... this may not be the true frequency of the system
            best_frequency = frequency[power==np.amax(power)]
            period=1/best_frequency.value[0]
            plt.errorbar((MJD%period)/period, mag, yerr=mage, fmt='.k')
            plt.savefig("ZTFphasefold.png")
            plt.clf()
            
            # and if you want to get into periodograms and period searching more seriously, PLEASE read these
            # https://ui.adsabs.harvard.edu/abs/2015ApJ...812...18V/abstract
            # https://iopscience.iop.org/article/10.3847/1538-4365/aab766     <- particularly this one. it's fantastic
            
            # lastly I recommend investigating the BLS search, typical for eclipsing systems/exoplanets
        
        except: None
        

        # neaten file list
        exceptions=["files", "GaiaLoc.png", "PhotSED.pdf", "SDSS.png", "ZTFphot.pdf",
                    "AllPhot.png", "CDS_objectofinterest.url"]
    
        for filename in os.listdir(cwd):
            if filename not in exceptions:
                Path(cwd+"/"+filename).rename(cwd+"/files/"+filename)
                
        

        # delete all memory of variables so I never get confused with next loop. variables defined locally in a function are never saved globally
        try: del respAll
        except: None
        try: del respZTF
        except: None
        
        del RA; del Dec
        del RAdeg; del Decdeg
        del folder; del RADec
        
        try: del uv; del g; del r; del i; del z
        except: None
        
        try: del p1
        except: None
        try: del p2
        except: None
        try: del p3
        except: None
        try: del p4
        except: None
        try: del p5
        except: None
        try: del p6
        except: None
        try: del p7
        except: None
        try: del p8
        except: None
        try: del p9
        except: None
        try: del p10
        except: None
        try: del p11
        except: None
        
    
        # back to the start
        os.chdir(homedirectory)



# observing biases:
# asassn - 30s
# goto - 60s, maybe changed in data out of comissioning
# gaia - 4s
# ztf - 
# crts - 
# asassn - 

# todo:
# query cds references, e.g.: http://simbad.u-strasbg.fr/simbad/sim-id-refs?submit=sort+references&Ident=SDSS%20J053332.05%2B020911.5

