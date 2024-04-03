#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 23:00:00 2022

@author: james
"""

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.insert(0, "/home/astro/phrsjf/PeriodSearcher2/getScripts/getScripts")

class OverplotGaia(object):
    def plotGaia(RADec, BPRP, Abs_g, gmag):
        img = mpimg.imread('../../getScripts/Screenshot.png')
        
        RADecsplit=str(RADec).split(":")
        
        RADecTitle = RADecsplit[0] + ":" + RADecsplit[1] + ":" + RADecsplit[2][:4] + " " + RADecsplit[2][-3:]+":"+RADecsplit[3]+":"+    RADecsplit[4][:4]     
        plt.clf()
        imgplot = plt.imshow(img, extent=[-0.75,2.2,17,5], aspect=0.16)
        plt.scatter(BPRP,Abs_g, marker="*", s=70, c='r')
        
        # rough ZZ ceti from the diagram of https://arxiv.org/pdf/2107.08458.pdf
        plt.plot([0.1,2.5514E-3], [11, 14.36], c='r', alpha=0.6)
        plt.plot([-0.08, 0.032], [14.1, 10.9], c='blue', alpha=0.6)
        
        
        #plt.xlim(-0.75,2.2)
        #plt.ylim(5,17)
        try:
            plt.title(str(RADecTitle) + ", gmag = " + str(gmag), fontsize=10)
        except:
            try:
                plt.title(str(RADecTitle) + ", gmag = " + str(gmag), fontsize=10)
            except: None
        
        
        plt.ylabel("Absolute G [mag]")
        plt.xlabel("BP - RP [mag]")
        plt.savefig("GaiaLoc.png", dpi=200)
        plt.clf()
        plt.close()
