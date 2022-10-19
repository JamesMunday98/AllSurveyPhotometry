# All Survey Time-Series Photometry
The package will obtain all time-series photometry associated with the object that is publicly available. All you need is the RA and Dec!

  
# Current list of things that can be accessed with the ease:  

| Survey      | Function  | Test Text     |
| :---        |    :----:   |          ---: |
| ATLAS        |    Forced photometry from the ATLAS survey     |          ---: |
| ASASSN        |    Find light curves in the full variable star catalogue     |          ---: |
| ASASSN Web        |    Autogenerate web search criteria for forced photometry   |          ---: |
| Catalina/CRTS        |    Obtain light curves   |          ---: |
| CDS        |    Obtain photometric SED for any search radius. Clickable link to CDS for the object   |          ---: |
| Gaia        |    All Gaia epoch photometry/spectra/RVS   |          ---: |
| :---        |    :----:   |          ---: |
| :---        |    :----:   |          ---: |
| :---        |    :----:   |          ---: |
| :---        |    :----:   |          ---: |
| :---        |    :----:   |          ---: |
| :---        |    :----:   |          ---: |
| :---        |    :----:   |          ---: |
| :---        |    :----:   |          ---: |
| :---        |    :----:   |          ---: |
| Paragraph   | Text        | And more      |  


  
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


# Missing a source that you want included?
Push a script to access the data and/or supply basic details of where the data is obtained.


