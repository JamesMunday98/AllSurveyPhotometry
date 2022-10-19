# AllSurveyPhotometry
The package will obtain all time-series photometry associated with the object that is publicly available. All you need is the RA and Dec!

  
# Current list of things that can be accessed with the ease:  
ATLAS  
ASASSN variable catalogue  
ASASSN Web  
Catalina/CRTS  
CDS  
Gaia  
Kepler  
NEOWISE  
Panstarrs  
PTF  
SDSS  
TESS  
WISE  
ZTF  


Extras:  
Gaia HR plotter  
Plot all photometry  
Basic Lomb-Scargle  





# What do you need to do?
- Clone this repository and look at getScripts/getpwds.py  . Here you need to include your [ATLAS](https://fallingstar-data.com/forcedphot/) login details and your [IRSA](https://irsa.ipac.caltech.edu/Missions/ztf.html) account details.
- Read in your RA and Dec (or list of these)
- (Optional) Inspect the search options for each survey
 

## Required (non standard) packages:
[Multiprocess](https://pypi.org/project/multiprocess/)  
[Astropy](https://docs.astropy.org/en/stable/install.html)  
[Beautiful Soup](https://pypi.org/project/beautifulsoup4/)  
[PIL](https://pypi.org/project/Pillow/)  
[Lightkurve](https://docs.lightkurve.org/about/install.html)  
[Astroquery](https://astroquery.readthedocs.io/en/latest/)  
[jdcal](https://pypi.org/project/jdcal/)  
