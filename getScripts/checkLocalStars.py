#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 10:01:23 2022

@author: james
"""

from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
import numpy as np

class checkLocalStars(object):
    def find_star_in_gaia_edr3(RAdeg,Decdeg):
        v = Vizier(columns=['Gmag'])
        obj=coord.SkyCoord(ra=RAdeg, dec=Decdeg, unit=(u.deg, u.deg), frame='icrs')
        
        # grab star information
        result = v.query_region(obj, radius=0.2 * u.arcsec, catalog='I/350/gaiaedr3')[0] # gaia edr3   https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/350&-out.max=50&-out.form=HTML%20Table&-out.add=_r&-out.add=_RAJ,_DEJ&-sort=_r&-oc.form=sexa
        return [obj, result['Gmag']]
        
        
    def localTESS(obj, star_mag):
        v = Vizier(columns=['Gmag'], column_filters={'Gmag': '<19.5'})
        
        # I use massive to make sure there are no nearby super bright stars that could cause saturation/diffraction spikes
        result_TESS = v.query_region(obj, radius=22*2 * u.arcsec, catalog='I/350/gaiaedr3')[0] # tess 21 arcsec per pixel, +2 for proper motion
        result_TESS_MASSIVE = v.query_region(obj, radius=22*7 * u.arcsec, catalog='I/350/gaiaedr3')[0] # tess 21 arcsec per pixel, +2 for proper motion
        
        for count, row in enumerate(result_TESS): 
            if row['Gmag'] == star_mag:  result_TESS.remove_row(count)
        
        for count, row in enumerate(result_TESS_MASSIVE): 
            if row['Gmag'] == star_mag:  result_TESS_MASSIVE.remove_row(count)
        
        
        if result_TESS_MASSIVE:
            if np.all(star_mag < result_TESS['Gmag']-0.5) and star_mag<=18 and np.all(result_TESS_MASSIVE['Gmag'] > 10):  # if it is easily the brightest star there
                print("passed checks with other stars");
                return "Good"
            elif np.all(star_mag < result_TESS['Gmag']+0.5) and star_mag<=18 and np.all(result_TESS_MASSIVE['Gmag'] > 10): # elif it is half a mag brighter than the others in the local area
                print("maybe, take care, not included in the combined period search")
                return "Maybe"
            else: # don't try, too much faff if dimmer than 18 for little gain
                print("nope")
                return "No"
        else: 
            #print("only star in field, good to go")
            return "Good"


    def localK2(obj,star_mag):
        v = Vizier(columns=['Gmag'], column_filters={'Gmag': '<19.5'})
        # I use massive to make sure there are no nearby super bright stars that could cause saturation/diffraction spikes
        result_K2 = v.query_region(obj, radius=4.98*2 * u.arcsec, catalog='I/350/gaiaedr3',column_filters={'Gmag': '<19.5'})[0] # Kepler 3.98 arcsec per pixel, +2 for proper motion
        result_K2_MASSIVE = v.query_region(obj, radius=4.48*7 * u.arcsec, catalog='I/350/gaiaedr3',column_filters={'Gmag': '<19.5'})[0] # K2 3.98 arcsec per pixel
        
        for count, row in enumerate(result_K2): 
            if row['Gmag'] == star_mag:
                result_K2.remove_row(count)
        
        for count, row in enumerate(result_K2_MASSIVE): 
            if row['Gmag'] == star_mag:
                result_K2_MASSIVE.remove_row(count)
        
        
        if result_K2_MASSIVE:
            if np.all(star_mag < result_K2['Gmag']-0.5) and star_mag<=18 and np.all(result_K2_MASSIVE['Gmag'] > 10):  # if it is easily the brightest star there
                #print("passed checks with other stars")
                return "Good"
            elif np.all(star_mag < result_K2['Gmag']+0.5) and star_mag<=18 and np.all(result_K2_MASSIVE['Gmag'] > 10): # elif it is half a mag brighter than the others in the local area
                #print("maybe, take care, not included in the combined period search")
                return "Maybe"
            else: # don't try, too much faff
                #print("nope")
                return "No"
        else:
            #print("only star, fine to go"); 
            return "Good"
    
    
    
    def localATLAS(obj,star_mag):
        v = Vizier(columns=['Gmag'], column_filters={'Gmag': '<19.5'})
        # I use massive to make sure there are no nearby super bright stars that could cause saturation/diffraction spikes
        result_ATLAS = v.query_region(obj, radius=1.86*5 * u.arcsec, catalog='I/350/gaiaedr3',column_filters={'Gmag': '<19.5'})[0] # Kepler 3.98 arcsec per pixel, +2 for proper motion
        result_ATLAS_MASSIVE = v.query_region(obj, radius=1.86*14 * u.arcsec, catalog='I/350/gaiaedr3',column_filters={'Gmag': '<19.5'})[0] # K2 3.98 arcsec per pixel
        
        for count, row in enumerate(result_ATLAS): 
            if row['Gmag'] == star_mag:
                result_ATLAS.remove_row(count)
        
        for count, row in enumerate(result_ATLAS_MASSIVE): 
            if row['Gmag'] == star_mag:
                result_ATLAS_MASSIVE.remove_row(count)
        
        
        if result_ATLAS_MASSIVE:
            if np.all(star_mag < result_ATLAS['Gmag']-1) and np.all(result_ATLAS_MASSIVE['Gmag'] > 13):  # if it is easily the brightest star there
                #print("passed checks with other stars")
                return "Good"
            elif np.all(star_mag < result_ATLAS['Gmag']+0.5) and np.all(result_ATLAS_MASSIVE['Gmag'] > 13): # elif it is half a mag brighter than the others in the local area
                #print("maybe, take care, not included in the combined period search")
                print("stopped atlas maybe")
                return "Maybe"
            else: # don't try, too much faff
                print("stopped atlas bad")
                #print("nope")
                return "No"
        else:
            #print("only star, fine to go"); 
            return "Good"
        
        
# We have  2347 out of  2784 remaining (progress =  15.7 % )
#('20:05:51.08095021', '-29:35:00.88146411')
#301.46283729253634 -29.583578184473748
#attempt to get argmin of an empty sequence


