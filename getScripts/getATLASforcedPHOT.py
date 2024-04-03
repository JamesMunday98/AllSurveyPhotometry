#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 11:16:23 2022

@author: james
"""

#!/usr/bin/env python3
import os, sys, re, time, requests
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clip

class getATLASforcedPHOT(object):
    def ATLAS(username, password, RA, Dec, reference_epoch, pmra, pmdec, minMJD=50000): #pm in mas/yr
        BASEURL = "https://fallingstar-data.com/forcedphot"
        # BASEURL = "http://127.0.0.1:8000"

        if os.environ.get('ATLASFORCED_SECRET_KEY'):
            token = os.environ.get('ATLASFORCED_SECRET_KEY')
            #print('Using stored token')
        else:
            data = {'username':username, 'password': password}

            resp = requests.post(url=f"{BASEURL}/api-token-auth/", data=data)

            if resp.status_code == 200:
                token = resp.json()['token']
                #print(f'Your token is {token}')
                #print('Store this by running/adding to your .zshrc file:')
                #print(f'export ATLASFORCED_SECRET_KEY="{token}"')
            else:
                print(f'ERROR {resp.status_code}')
                print(resp.text)
                sys.exit()


        headers = {'Authorization': f'Token {token}', 'Accept': 'application/json'}

        task_url = None
        while not task_url:
            with requests.Session() as s:
                # alternative to token auth
                # s.auth = ('USERNAME', 'PASSWORD')
                resp = s.post(f"{BASEURL}/queue/", headers=headers, data={
                    'ra': RA, 'dec': Dec, 'mjd_min': minMJD, "propermotion_ra": pmra,
                    "propermotion_dec": pmdec,  "radec_epoch_year": reference_epoch,
                    'send_email': False, "use_reduced": True})

                if resp.status_code == 201:  # successfully queued
                    task_url = resp.json()['url']
                    #print(f'The task URL is {task_url}')
                elif resp.status_code == 429:  # throttled
                    message = resp.json()["detail"]
                    #print(f'{resp.status_code} {message}')
                    t_sec = re.findall(r'available in (\d+) seconds', message)
                    t_min = re.findall(r'available in (\d+) minutes', message)
                    if t_sec:
                        waittime = int(t_sec[0])
                    elif t_min:
                        waittime = int(t_min[0]) * 60
                    else:
                        waittime = 10
                    #print(f'Waiting {waittime} seconds')
                    time.sleep(waittime)
                else:
                    print(f'ERROR {resp.status_code}')
                    print(resp.text)
                    sys.exit()


        result_url = None
        taskstarted_printed = False
        found_val=0
        while not result_url: #and found_val<2:
            with requests.Session() as s:
                resp = s.get(task_url, headers=headers,timeout=7200)

                if resp.status_code == 200:  # HTTP OK
                    if resp.json()['finishtimestamp']:
                        result_url = resp.json()['result_url']
                        #print(f"Task is complete with results available at {result_url}")
                        #found_val+=2
                    elif resp.json()['starttimestamp']:
                        if not taskstarted_printed:
                            #print(f"Task is running (started at {resp.json()['starttimestamp']})")
                            taskstarted_printed = True
                        time.sleep(2)
                    else:
                        #print(f"Waiting for job to start (queued at {resp.json()['timestamp']})")
                        time.sleep(4)
                else:
                    print(f'ERROR {resp.status_code}')
                    print(resp.text)
                    sys.exit()


        with requests.Session() as s:
            textdata = s.get(result_url, headers=headers, timeout=7200).text

            # if we'll be making a lot of requests, keep the web queue from being
            # cluttered (and reduce server storage usage) by sending a delete operation
            #s.delete(task_url, headers=headers).json()

        dfresult = pd.read_csv(StringIO(textdata.replace("###", "")), delim_whitespace=True)

        return dfresult

    def plot_and_save_data(a, string):
        # flux is the usable thing, and all errors are valid as it is difference imaging
        # 30s exposure
        # c for cyan, o for orange
        MJD=a["MJD"];   filt=a["F"];    uJy=a["uJy"];    duJy=a["duJy"]
        RA=a["RA"];   Dec=a["Dec"];    mag=a["m"];    dmag=a['dm'];    sigma5=a['mag5sig']

        mask5sigma=((mag<=sigma5) & (dmag<0.35) & (mag>8) & (mag<19.5))

        MJD=MJD[mask5sigma];  filt=filt[mask5sigma];  uJy=uJy[mask5sigma];   duJy=duJy[mask5sigma];   RA=RA[mask5sigma]
        Dec=Dec[mask5sigma];   mag=mag[mask5sigma];   dmag=dmag[mask5sigma];   sigma5=sigma5[mask5sigma]

        o=filt=="o";  c=filt=="c"

        MJDo=MJD[o];  uJyo=uJy[o];   duJyo=duJy[o];   mago=mag[o];   dmago=dmag[o]
        MJDc=MJD[c];   uJyc=uJy[c];   duJyc=duJy[c];   magc=mag[c];  dmagc=dmag[c]


        if string=="AddToOriginal":
            try:
                #mjdFileo, magFileo, dmagFileo, uJyFileo, duJyFileo = np.loadtxt("ATLAS_filtO_flux_and_err.dat", unpack=True)
                mjdFileo, magFileo, dmagFileo = np.loadtxt("ATLAS_o.dat", unpack=True)

                MJDo=np.concatenate((MJDo, mjdFileo))
                mago=np.concatenate((mago,magFileo))
                dmago=np.concatenate((dmago,dmagFileo))
                #uJyo=np.concatenate((uJyo,uJyFileo))
                #duJyo=np.concatenate((duJyo,duJyFileo))
            except: None

            try:
                #mjdFilec, magFilec, dmagFilec, uJyFilec, duJyFilec = np.loadtxt("ATLAS_filtC_flux_and_err.dat", unpack=True)
                mjdFilec, magFilec, dmagFilec = np.loadtxt("ATLAS_c.dat", unpack=True)

                MJDc=np.concatenate((MJDc, mjdFilec))
                magc=np.concatenate((magc,magFilec))
                dmagc=np.concatenate((dmagc,dmagFilec))
                #uJyc=np.concatenate((uJyc,uJyFilec))
                #duJyc=np.concatenate((duJyc,duJyFilec))
            except: None

            maxit=1
        else:  maxit=2


        #masked_o = sigma_clip(mago,sigma_lower=3.2, sigma_upper=3.5, cenfunc='median', maxiters=maxit).mask
        #masked_c = sigma_clip(magc,sigma_lower=3.2, sigma_upper=3.5, cenfunc='median', maxiters=maxit).mask

        ##plt.errorbar(MJDc[~masked_c],magc[~masked_c], yerr=dmagc[~masked_c], fmt=".g", label="c")
        ##plt.errorbar(MJDo[~masked_o],mago[~masked_o], yerr=dmago[~masked_o], fmt=".b", label="o")
        ##plt.xlabel("MJD");   plt.ylabel("Magnitude");   plt.gca().invert_yaxis()
        ##plt.savefig("ATLASmag.png")
        ##plt.clf()

        #np.savetxt("ATLAS_c.dat", np.array([MJDc[~masked_c],magc[~masked_c],dmagc[~masked_c]]).T)#,uJyc[~masked_c], duJyc[~masked_c]
        #np.savetxt("ATLAS_o.dat", np.array([MJDo[~masked_o],mago[~masked_o],dmago[~masked_o]]).T)#,uJyo[~masked_o], duJyo[~masked_o] ]).T)
        
        np.savetxt("ATLAS_c.dat", np.array([MJDc,magc,dmagc]).T)#,uJyc[~masked_c], duJyc[~masked_c]
        np.savetxt("ATLAS_o.dat", np.array([MJDo,mago,dmago]).T)#,uJyo[~masked_o], duJyo[~masked_o] ]).T)
        
        # flux is in microjanskys, same for error. RA Dec are in degrees
        #plt.close()
        return len(MJD)
