#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 18:32:03 2022

@author: james
"""
import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk

RADec = "18:36:51.92 +51:20:24.438"
radius = 5
def TessRunOnceWithCBV():
    for time in ["fast", "short", "long"]:
    
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
                        lc=raw_lc.to_corrector(method="cbv").correct() #cbv #pld
                        
                        
                        lc=lc.remove_outliers(sigma=5)
                        lc=lc.normalize()
                        
                        all_lcs.append(lc)
                        
                        
                    lc_obj=lk.LightCurveCollection(all_lcs).stitch()
                    ax=lc_obj.scatter(title=time)
                    ax.figure.savefig('AllLC'+str(time)+'.png')
                    plt.close()
                    lcfound=True
        
            elif len(search_result_lc) > 1:
                lc_obj=search_result_lc.stitch().remove_outliers().remove_nans()   # the flux here is meant to automatically be PDCSAP
                lcfound=True
                
                
            if lcfound==True:
                lc_obj.plot()
                
                plt.clf()
                np.savetxt("TESSdata_"+str(time)+".csv", np.array([lc_obj.time + 2457000- 2400000.5, lc_obj.flux, lc_obj.flux_err]).T, fmt="%s")
                # time is output as BJD 


t, y, dy = np.loadtxt("TESSdata_short.csv", unpack=True)
plt.errorbar(t,y,yerr=dy, ls=' ', fmt='.k')
plt.title("All TESS")
plt.show()

from astropy import units as u
from astropy.timeseries import LombScargle
freq = 1/np.linspace(5.75*60/86400, 6.5*60/86400, 40000)
ls = LombScargle(t*u.day, y, dy)
power = ls.power(freq*(1/u.day))
plt.plot((1/freq)*24*60, power)

print(np.amax(power))
v=np.argwhere(power==np.amax(power))
P = (1/freq[np.argwhere(power==np.amax(power))[0][0]])*24*60

plt.xlim(5,7)
plt.axvline(P, c='r')
plt.show()


from PyAstronomy.pyasl import binningx0dt
phase = (t%P)/P

r1, dt1 = binningx0dt(phase, y, yerr=dy, nbins=100, x0=np.amin(phase), useBinCenter=True)


binnedBJDfold = r1.T[0]
binnedflux = r1.T[1]
binnedfluxerrs = r1.T[2]



plt.errorbar(binnedBJDfold, binnedflux, yerr=binnedfluxerrs, ls=' ', fmt='.k')

plt.axhline(1, c='grey', ls='--')

plt.title(str(P)+" minutes")
plt.show()



from scipy.optimize import curve_fit
import pylab as plt

guess_freq = 1
guess_amplitude = 0
guess_phase = 0
guess_offset = 1

p0=[guess_freq, guess_amplitude,
    guess_phase, guess_offset]

# create the function we want to fit
def my_sin(x, freq, amplitude, sinphase, offset):
    return np.sin(x * freq + sinphase) * amplitude + offset

# now do the fit
fit = curve_fit(my_sin, binnedBJDfold, binnedflux,sigma=binnedfluxerrs, p0=p0)


# recreate the fitted curve using the optimized parameters
data_fit = my_sin(binnedBJDfold, *fit[0])

print(fit[0][2]%1)

plt.errorbar(binnedBJDfold, binnedflux, yerr=binnedfluxerrs, ls=' ', fmt='.k')
plt.plot(binnedBJDfold, data_fit, c='k')
plt.title(str(P)+" minutes")
plt.axhline(1, c='grey')
plt.ylabel("Flux")
plt.xlabel("Phase")
plt.show()








# short = 2min
# long = 30/10min depending on sector
# fast = 20s



        