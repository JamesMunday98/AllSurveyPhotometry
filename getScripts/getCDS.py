#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 13:57:06 2022

@author: james
"""

class CDS(object):
    def CDSurl(RAdeg, Decdeg):
        url="http://cdsportal.u-strasbg.fr/?target="
        url+=str(RAdeg)
        url+="%20"
        
        if Decdeg < 0:
            url+=str(Decdeg)
        else:
            url+="%2B"
            url+=str(Decdeg)
        
        return url

    
    def createShortcut(url, destination):
        text = '[InternetShortcut]\nURL={}\nIconIndex=0'.format(url)
        with open(destination + 'CDS_objectofinterest.url', 'w') as fw:
            fw.write(text)
        return


