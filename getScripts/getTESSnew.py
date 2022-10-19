#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 23:21:30 2022

@author: james
"""

import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk
#from astropy.timeseries import LombScargle

plt.ioff()

# https://docs.lightkurve.org/tutorials/2-creating-light-curves/2-1-combining-multiple-quarters.html
# https://heasarc.gsfc.nasa.gov/docs/tess/Aperture-Photometry-Tutorial.html


#transit_mask = corrected_lc.create_transit_mask(period=[3.2277, 7.3745],
                                                #duration=[0.25, 0.25],
                                                #transit_time=[2385.6635, 2389.9635])
#https://docs.lightkurve.org/tutorials/2-creating-light-curves/2-3-k2-pldcorrector.html?highlight=to_corrector


class TESS(object):
    def get_tess(RADec, time, radius=5, ignore_any_dodgyness=False):
        print(RADec)
        search_result_lc=lk.search_lightcurve(RADec, radius=radius,author="TESS", exptime=time) # https://docs.lightkurve.org/reference/api/lightkurve.search_lightcurve.html
        lcfound=False
        
        # inspect tpf: https://docs.lightkurve.org/tutorials/1-getting-started/interactively-inspecting-data.html
        if len(search_result_lc)==0:
            search_result = lk.search_targetpixelfile(RADec, mission="TESS", exptime=time)
            if search_result:
                #print(search_result)
                try:
                    all_tpf = search_result.download_all("hard")
                except Exception as e:
                    print("YIKES");  print(e)
                    all_tpf=search_result[0].download("hard")    # THIS IS BAD. I think it only comes in when there's a ram issue
                all_lcs=[]
                
                if all_tpf:  # if this is not empty
                    for tpf in all_tpf:
# =============================================================================
#                         print("###################################")
#                         print(tpf)
#                         print(tpf.flux.value)
#                         
#                         masked=tpf.flux_err.value > 0
#                         print(masked)
#                         print(len(tpf))
#                         #tpf=tpf[tpf.flux_err.value > 0]
# =============================================================================
                        #try:
                        bkg = tpf.get_bkg_lightcurve()
                        for z in range(3):
                            medianbkg=np.median(bkg.flux).value
                            sigmabkg=np.std(bkg.flux).value
                            mask=((bkg.flux.value > (medianbkg-3.5*sigmabkg)) & (bkg.flux.value < (medianbkg+2*sigmabkg)))
                            bkg=bkg[mask]
                            tpf=tpf[mask]
                            
                        
                        
                        aper = tpf.create_threshold_mask()
                        raw_lc = tpf.to_lightcurve(aperture_mask=aper)
                        raw_lc = raw_lc.remove_nans()
                        try:
                            lc=raw_lc.to_corrector(method="sff").correct() #cbv #pld
                        except:
                            try:
                                lc=raw_lc.to_corrector(method="cbv").correct()
                            except:
                                None
                        
                        try:
                            lc=lc.remove_outliers(sigma=5)
                            lc=lc.normalize()
                            
                            all_lcs.append(lc)
                        except: None
                        
                    try:
                        lc_obj=lk.LightCurveCollection(all_lcs).stitch()
                        ax=lc_obj.scatter(title=time)
                        ax.figure.savefig('AllLC'+str(time)+'.png')
                        plt.close()
                        lcfound=True
                    except: None

            elif len(search_result_lc) > 1:
                lc_obj=search_result_lc.stitch().remove_outliers().remove_nans()   # the flux here is meant to automatically be PDCSAP
                lcfound=True
                
                
                
                
                
                
                
            
            if lcfound==True:
                np.savetxt("TESSdata_"+str(time)+".csv", np.array([lc_obj.time + 2457000- 2400000.5, lc_obj.flux, lc_obj.flux_err]).T, fmt="%s")
                # time is output as BJD 
            
            
            
            
            
            
            
#old code here to do periodograms with lightkurve         
                #try: # this will only break if lc_obj does not exist
                    ##ts=lc_obj.time
                    ##fluxes=lc_obj.flux
                    ##fluxes_e=lc_obj.flux_err
                    
                    ##times_to_sample=np.linspace(1/86400,1/250,1000)*u.Hz
                    ##ls = LombScargle(ts,fluxes,fluxes_e)
                    ##power=ls.power(times_to_sample)
                    
                    
                    ##plt.plot((1/times_to_sample)/86400, power)
                    ##plt.show()
                
                
                
                
                    #periodogram=lc_obj.to_periodogram("bls", minimum_period=250/86400, maximum_period=2)
                    #plt.clf()
                    #print("I GOT HERE")
                    #periodogram.plot().figure.savefig('Periodogram_'+str(time)+'.png')
                    #period=periodogram.period_at_max_power
                    #folded_lc=lc_obj.fold(period)
                    #folded_lc.scatter(title=str(period) + "   " + str(time)).savefig('folded_lc_'+str(time)+'.png')
                    #plt.close()
                    ##folded_lc.plot_river(method='sigma')
                #except Exception as exc:
                #    print(exc)
            
            
            
            
            
            # for kepler data, https://docs.lightkurve.org/tutorials/2-creating-light-curves/2-3-k2-sffcorrector.html
            
            
            # limitation: might struggle with crowded fields... you can use difference imaging with lightkurve https://docs.lightkurve.org/tutorials/3-science-examples/periodograms-verifying-the-location-of-a-signal.html
        
    


#TESS.get_tess("23:53:00.9 -38:51:47.67", "short")
            