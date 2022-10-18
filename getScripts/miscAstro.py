#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 21:01:18 2022

@author: james
"""

from astropy import units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.time import Time
import os
from PIL import Image
import shutil

class miscAstro(object):
    def ra_dec_deg_to_hr(RAdeg,Decdeg):
        RA=Angle(RAdeg/15*u.deg).to_string(unit=u.degree, sep=':')
        Dec = Angle(Decdeg*u.deg).to_string(unit=u.degree, sep=':')
        return RA, Dec
        
    
    def ra_dec_hr_to_deg(ra,dec):
        # string formats "00:00:00 +00:00:00"
        c = SkyCoord(ra,dec, frame='icrs', unit=(u.hourangle, u.deg))
        return c.ra.deg, c.dec.deg
    
    
    def jd_corr(mjd, ra, dec, loc, jd_type='bjd'):
        # in this script, I work in BJD UTC. Never is it BJD TDB. I want to change this
        target=SkyCoord(ra, dec,unit=(u.deg, u.deg), frame='icrs')
        jd=Time(mjd,format='mjd',scale='utc',location=loc)
        if jd_type == 'bjd':
          corr=jd.light_travel_time(target,kind='barycentric')
        elif jd_type == 'hjd':
          corr=jd.light_travel_time(target,kind='heliocentric')
        new_jd = jd + corr
        return new_jd
  
    
    # assisted by format of https://github.com/WarwickAstro/time-conversions/blob/master/convert_times.py
    def getLightTravelTimes(ra, dec, time_to_correct):
        target = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
        ltt_bary = time_to_correct.light_travel_time(target)
        ltt_helio = time_to_correct.light_travel_time(target, 'heliocentric')
        return ltt_bary, ltt_helio


    def hjd_to_bjd(RAdeg,Decdeg,HJD,tel):
        # input times, subtract off the HJD correction frmo the JD time
        all_HJD = Time(HJD, format='jd', scale='utc', location=tel)
        _, ltt_helio = miscAstro.getLightTravelTimes(RAdeg, Decdeg, all_HJD)
        all_JD = Time(all_HJD.utc - ltt_helio, format='jd', scale='utc', location=tel)
        
        # add the ltt correction to get barycentric time, save results
        ltt_bary, _ = miscAstro.getLightTravelTimes(RAdeg, Decdeg, all_JD)
        #BJD = (all_JD.tdb + ltt_bary).value      #tdb
        BJD = (all_JD + ltt_bary).value            #Whatever ASASSN has given the units in, whether tdb or utc, this maintains those units    -    I assume what is given by ASASSN is HJD UTC (I should check this to see if it improves or worsens a period search, the difference is >30s time offset so will impact things a lot)
    
        #print((BJD-HJD)*86400)   # looks right
        return BJD

      
        
    def MergeImages(list_im, filename):
        # robbed from stackoverflow with my own comments added: 
        # https://stackoverflow.com/questions/30227466/combine-several-images-horizontally-with-python
        
        imgs    = [ Image.open(i) for i in list_im ]
        # pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
        min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
        imgs_comb = np.hstack( [(np.asarray( i.resize(min_shape) ) for i in imgs ) ])
        
        ## save that beautiful picture
        #imgs_comb = Image.fromarray( imgs_comb)
        #imgs_comb.save( 'GaiaLoc_Periodogram.png' )
        
        # for a vertical stacking it is simple: use vstack
        imgs_comb = np.vstack([ (np.asarray( i.resize(min_shape) ) for i in imgs ) ])
        imgs_comb = Image.fromarray( imgs_comb)
        imgs_comb.save( filename )
            
        
    
    def append_images(images, direction='horizontal',
                      bg_color=(255,255,255), aligment='center'):
        """
        Appends images in horizontal/vertical direction.
    
        Args:
            images: List of PIL images
            direction: direction of concatenation, 'horizontal' or 'vertical'
            bg_color: Background color (default: white)
            aligment: alignment mode if images need padding;
               'left', 'right', 'top', 'bottom', or 'center'
    
        Returns:
            Concatenated image as a new PIL image object.
        """
        widths, heights = zip(*(i.size for i in images))
    
        if direction=='horizontal':
            new_width = sum(widths)
            new_height = max(heights)
        else:
            new_width = max(widths)
            new_height = sum(heights)
    
        new_im = Image.new('RGB', (new_width, new_height), color=bg_color)
    
    
        offset = 0
        for im in images:
            if direction=='horizontal':
                y = 0
                if aligment == 'center':
                    y = int((new_height - im.size[1])/2)
                elif aligment == 'bottom':
                    y = new_height - im.size[1]
                new_im.paste(im, (offset, y))
                offset += im.size[0]
            else:
                x = 0
                if aligment == 'center':
                    x = int((new_width - im.size[0])/2)
                elif aligment == 'right':
                    x = new_width - im.size[0]
                new_im.paste(im, (x, offset))
                offset += im.size[1]
    
        return new_im
    
# =============================================================================
#     def MergeIms_Folded_periodsZTF():
#         images=[Image.open('best1.png'),Image.open('best2.png'),Image.open('best3.png')]
#         images2=[Image.open('best4.png'),Image.open('best5.png'),Image.open('best6.png')]
#         images3=[Image.open('best7.png'),Image.open('best8.png'),Image.open('best9.png')]
#         combo_1 = miscAstro.append_images(images, direction='horizontal')
#         combo_2 = miscAstro.append_images(images2, direction='horizontal')#, aligment='top',  bg_color=(220, 140, 60))
#     
#         combo_3 = miscAstro.append_images(images3, direction='horizontal')
#         combo_4 = miscAstro.append_images([combo_1, combo_2,combo_3], direction='vertical')
#         combo_4.save('bestPeriods.png')
#         os.remove('best1.png'); os.remove('best2.png')
#         os.remove('best3.png'); os.remove('best4.png')
#         os.remove('best5.png'); os.remove('best6.png')
#         os.remove('best7.png'); os.remove('best8.png')
#         os.remove('best9.png')
# =============================================================================
    
    def MergeIms_Periodograms(filein1, filein2, fileout):
        images=[Image.open(filein1),Image.open(filein2)]
        combo_1 = miscAstro.append_images(images, direction='horizontal')
        combo_1.save(fileout)
        os.remove(filein1)
        os.remove(filein2)
        
    
        
    def MergeIms_Folded_periodsZTF(output,stringtoadd=''):
        images=[Image.open('best1'+stringtoadd+'.png'),Image.open('best2'+stringtoadd+'.png'),Image.open('best3'+stringtoadd+'.png')]
        images2=[Image.open('best4'+stringtoadd+'.png'),Image.open('best5'+stringtoadd+'.png'),Image.open('best6'+stringtoadd+'.png')]
        images3=[Image.open('best7'+stringtoadd+'.png'),Image.open('best8'+stringtoadd+'.png'),Image.open('best9'+stringtoadd+'.png')]
        combo_1 = miscAstro.append_images(images, direction='horizontal')
        combo_2 = miscAstro.append_images(images2, direction='horizontal')#, aligment='top',  bg_color=(220, 140, 60))
    
        combo_3 = miscAstro.append_images(images3, direction='horizontal')
        combo_4 = miscAstro.append_images([combo_1, combo_2,combo_3], direction='vertical')
        combo_4.save('bestPeriods'+output+'.png')
        os.remove('best1'+stringtoadd+'.png'); os.remove('best2'+stringtoadd+'.png')
        os.remove('best3'+stringtoadd+'.png'); os.remove('best4'+stringtoadd+'.png')
        os.remove('best5'+stringtoadd+'.png'); os.remove('best6'+stringtoadd+'.png')
        os.remove('best7'+stringtoadd+'.png'); os.remove('best8'+stringtoadd+'.png')
        os.remove('best9'+stringtoadd+'.png')
    
        
    def MergeIms_Folded_periodsTESS():
        images=[Image.open('TESSbest1.png'),Image.open('TESSbest2.png'),Image.open('TESSbest3.png')]
        images2=[Image.open('TESSbest4.png'),Image.open('TESSbest5.png'),Image.open('TESSbest6.png')]
        images3=[Image.open('TESSbest7.png'),Image.open('TESSbest8.png'),Image.open('TESSbest9.png')]
        combo_1 = miscAstro.append_images(images, direction='horizontal')
        combo_2 = miscAstro.append_images(images2, direction='horizontal')#, aligment='top',  bg_color=(220, 140, 60))
    
        combo_3 = miscAstro.append_images(images3, direction='horizontal')
        combo_4 = miscAstro.append_images([combo_1, combo_2,combo_3], direction='vertical')
        combo_4.save('bestPeriodsTESS.png')
        os.remove('TESSbest1.png'); os.remove('TESSbest2.png')
        os.remove('TESSbest3.png'); os.remove('TESSbest4.png')
        os.remove('TESSbest5.png'); os.remove('TESSbest6.png')
        os.remove('TESSbest7.png'); os.remove('TESSbest8.png')
        os.remove('TESSbest9.png')
        
    def remDir(folder):
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except: None

        # yeet the old stuff
        os.rmdir(folder)

    



