#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:08:18 2022

@author: james
"""

from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table

#https://www.cosmos.esa.int/web/gaia-users/archive/datalink-products#datalink_jntb_get_all_prods

class GaiaDatalink(object):
    def extract_dl_ind(datalink_dict, key, RVS, figsize = [15,5], fontsize = 12, linewidth = 2, show_legend = True, show_grid = True):
        ""
        "Extract individual DataLink products and export them to an Astropy Table"
        ""
        dl_out  = datalink_dict[key][0].to_table()
        if 'time' in dl_out.keys():
            GaiaDatalink.plot_e_phot(dl_out, colours  = ['green', 'red', 'blue'], title = 'Epoch photometry', fontsize = fontsize, show_legend = show_legend, show_grid = show_grid, figsize = figsize)
        if 'wavelength' in dl_out.keys():
            if len(dl_out) == 343:  title = 'XP Sampled'
            if len(dl_out) == 2401: title = 'RVS'
            GaiaDatalink.plot_sampled_spec(dl_out, RVS,color = 'blue', title = title, fontsize = fontsize, show_legend = False, show_grid = show_grid, linewidth = linewidth, legend = '', figsize = figsize)
        return dl_out
    
    
    def plot_e_phot(inp_table, colours  = ['green', 'red', 'blue'], title = 'Epoch photometry', fontsize = 12, show_legend = True, show_grid = True, figsize = [15,5]):
        ""
        "Epoch photometry plotter. 'inp_table' MUST be an Astropy-table object."
        ""
        fig      = plt.figure(figsize=figsize)
        xlabel   = f'JD date [{inp_table["time"].unit}]'
        ylabel   = f'magnitude [{inp_table["mag"].unit}]'
        gbands   = ['G', 'RP', 'BP']
        colours  = iter(colours)
    
        plt.gca().invert_yaxis()
        for band in gbands:
            phot_set = inp_table[inp_table['band'] == band]
            plt.plot(phot_set['time'] + 2455197.5- 2400000.5, phot_set['mag'], 'o', label = band, color = next(colours))
        GaiaDatalink.make_canvas(title = title, xlabel = xlabel, ylabel = ylabel, fontsize= fontsize, show_legend=show_legend, show_grid = show_grid)
        plt.savefig("gaiaEpochPhot.png")
        plt.clf()
        #plt.show()
    
    
    def plot_sampled_spec(inp_table, RVS,color = 'blue', title = '', fontsize = 14, show_legend = True, show_grid = True, linewidth = 2, legend = '', figsize = [12,4], show_plot = True):
        ""
        "RVS & XP sampled spectrum plotter. 'inp_table' MUST be an Astropy-table object."
        ""
        if show_plot:
            fig      = plt.figure(figsize=figsize)
        xlabel   = f'Wavelength [{inp_table["wavelength"].unit}]'
        ylabel   = f'Flux [{inp_table["flux"].unit}]'
        plt.plot(inp_table['wavelength'], inp_table['flux'], '-', linewidth = linewidth, label = legend)
        GaiaDatalink.make_canvas(title = title, xlabel = xlabel, ylabel = ylabel, fontsize= fontsize, show_legend=show_legend, show_grid = show_grid)
        if show_plot:
            if RVS==True: plt.savefig("gaiaRVS.png")
            else: plt.savefig("gaiaSpectrum.png")
            plt.clf()
            #plt.show()
    
    
    def make_canvas(title = '', xlabel = '', ylabel = '', show_grid = False, show_legend = False, fontsize = 12):
        ""
        "Create generic canvas for plots"
        ""
        plt.title(title,    fontsize = fontsize)
        plt.xlabel(xlabel,  fontsize = fontsize)
        plt.ylabel(ylabel , fontsize = fontsize)
        plt.xticks(fontsize = fontsize)
        plt.yticks(fontsize = fontsize)
        if show_grid:
            plt.grid()
        if show_legend:
            plt.legend(fontsize = fontsize*0.75)
        
    
    
    def getData(source_id_input):
        query = f"SELECT source_id, ra, dec, pmra, pmdec, parallax \
        FROM gaiadr3.gaia_source \
        WHERE source_id = "+str(source_id_input)
        
        
        job     = Gaia.launch_job_async(query)
        results = job.get_results()
        
        
            
        retrieval_type = 'ALL'          # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
        data_structure = 'INDIVIDUAL'   # Options are: 'INDIVIDUAL', 'COMBINED', 'RAW'
        data_release   = 'Gaia DR3'
        
        
        datalink = Gaia.load_data(ids=results['source_id'], data_release = data_release, retrieval_type=retrieval_type, data_structure = data_structure, verbose = False, output_file = None)
        dl_keys  = [inp for inp in datalink.keys()]
        dl_keys.sort()
        
        
        dl_key=['XP_CONTINUOUS-Gaia DR3 '+ source_id_input+'.xml', 'EPOCH_PHOTOMETRY-Gaia DR3 '+ source_id_input+'.xml','MCMC_MSC-Gaia DR3 '+ source_id_input+'.xml','MCMC_GSPPHOT-Gaia DR3 '+ source_id_input+'.xml','RVS-Gaia DR3 '+ source_id_input+'.xml','XP_SAMPLED-Gaia DR3 '+ source_id_input+'.xml']
        for i in dl_key:
            if "RVS" in i:
                RVS=True
            else:
                RVS=False
            try:
                dl_out  = GaiaDatalink.extract_dl_ind(datalink, i,RVS, figsize=[20,7])
                ascii.write(dl_out,i[:-4]+".dat", overwrite=True)
            except:
                None

# =============================================================================
# 
# source_id_input="4911590910260264960"
# 
# GaiaDatalink.getData(source_id_input)
# 
# import os
# 
# 
# for file in os.listdir(os.getcwd()):
#     if "EPOCH_PHOTOMETRY-Gaia" in file and not ".xml" in file:
#         print(file)
#         t=ascii.read(file)
# print(t.colnames)
# 
# t=t[(t['mag']>0) & (t['rejected_by_variability']=="False") & (t['rejected_by_photometry']=="False")]
# 
# time=t['time'] + 2455197.5- 2400000.5
# band = t['band']
# mag=t['mag']
# 
# import numpy as np
# magrand=np.log10(t['flux'])
# 
# 
# plt.clf()
# mage=np.asarray(np.abs(-2.5/np.log(10) * t['flux_error']/t['flux']))
#     
# 
# 
# =============================================================================





