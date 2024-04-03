#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 23:21:30 2022

@author: james
"""

import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk
from astropy.timeseries import LombScargle

# https://docs.lightkurve.org/tutorials/2-creating-light-curves/2-1-combining-multiple-quarters.html
# https://heasarc.gsfc.nasa.gov/docs/tess/Aperture-Photometry-Tutorial.html


#transit_mask = corrected_lc.create_transit_mask(period=[3.2277, 7.3745],
                                                #duration=[0.25, 0.25],
                                                #transit_time=[2385.6635, 2389.9635])
#https://docs.lightkurve.org/tutorials/2-creating-light-curves/2-3-k2-pldcorrector.html?highlight=to_corrector



class K2(object):
    def get_K2(RADec, exptime, radius=3, ignore_any_dodgyness=False):
        search_result_lc=lk.search_lightcurve(RADec, radius=radius,author="Kepler", exptime=exptime) # https://docs.lightkurve.org/reference/api/lightkurve.search_lightcurve.html
        
        
        # inspect tpf: https://docs.lightkurve.org/tutorials/1-getting-started/interactively-inspecting-data.html
        if len(search_result_lc)==0:
            search_result = lk.search_targetpixelfile(RADec, mission="Kepler", exptime=exptime)
            if search_result:
                all_tpf = search_result.download_all()
                all_lcs=[]
                
                
                if all_tpf:  # if this is not empty
                    for tpf in all_tpf:
                        try:
                            lc=tpf.to_lightcurve(method="pld") #cbv #pld
                            #print(lc.info())
                            lc=lc.remove_nans()
                            
                            all_lcs.append(lc)
                            
                                
                            
                        except Exception as e:
                            print(e)
                            print("trying to deal with nans")
                            lc=tpf.to_lightcurve().remove_nans()
                            mask=~(np.isnan(lc["flux_err"]) | np.isnan(lc["flux"]))
                            lc=lc[mask].remove_outliers()
                            ax=lc.scatter()
                            ax.set_title(str(tpf))
                            if ignore_any_dodgyness == True:
                                ax.figure.savefig('K2_maybeDodgy_IgnoredFromPeriodogram.png')
                            else:
                                ax.figure.savefig('K2_maybeDodgy_butIncludedInPeriodogram.png')
                                all_lcs.append(lc)
                        
                        
                    lc_obj=lk.LightCurveCollection(all_lcs).stitch().remove_outliers()
                    ax=lc_obj.scatter(title=exptime)
                    ax.figure.savefig('K2_AllLC'+str(exptime)+'.png')
                        
            elif len(search_result_lc) > 1:
                lc_obj=search_result_lc.stitch().remove_outliers().remove_nans()   # the flux here is meant to automatically be PDCSAP
            plt.close()
            
            try:     np.savetxt("K2data_"+str(exptime)+".csv", np.array([lc_obj.time, lc_obj.flux, lc_obj.flux_err]).T, fmt="%s")
            except:  print("No K2 data found.\n")
            
            try: # this will only break if lc_obj does not exist
                #ts=lc_obj.time
                #fluxes=lc_obj.flux
                #fluxes_e=lc_obj.flux_err
                
                #times_to_sample=np.linspace(1/86400,1/250,1000)*u.Hz
                #ls = LombScargle(ts,fluxes,fluxes_e)
                #power=ls.power(times_to_sample)
                
                
                #plt.plot((1/times_to_sample)/86400, power)
                #plt.show()
            
            
            
            
                periodogram=lc_obj.to_periodogram("bls", minimum_period=8*60/86400, maximum_period=0.5)
                plt.clf()
                
                #plt_periodogram=periodogram.plot()
                #plt_periodogram.axhline(ls.false_alarm_level(1-0.9987), c="r")
                
                periodogram.plot().figure.savefig('K2_Periodogram_'+str(exptime)+'.png')
                
                period=periodogram.period_at_max_power
                folded_lc=lc_obj.fold(period)
                folded_lc.scatter(title=str(period) + "   " + str(exptime))
                folded_lc.scatter(title=str(period) + "   " + str(exptime)).figure.savefig('K2_folded_lc_'+str(exptime)+'.png')
                
                #folded_lc.plot_river(method='sigma')
                plt.close()
            except Exception as exc:
                None#print(exc)
            
            
            
            
            
            #lightkurve.LightCurve.create_transit_mask
            
            # for kepler data, https://docs.lightkurve.org/tutorials/2-creating-light-curves/2-3-k2-sffcorrector.html
            
            
            # limitation: might struggle with crowded fields... one can use difference imaging with lightkurve https://docs.lightkurve.org/tutorials/3-science-examples/periodograms-verifying-the-location-of-a-signal.html
        

    

    
            
