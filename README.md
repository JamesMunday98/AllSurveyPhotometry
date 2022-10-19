# All Survey Time-Series Photometry
The package will obtain all time-series photometry associated with the object that is publicly available. All you need is the RA and Dec!

  
# Current list that can be accessed with the ease:  

| Survey      | Function  | Comments     |
| :---        |    :----   |    :----       |
| ATLAS        |    Forced photometry from the ATLAS survey     |    Only query if there are no close contaminants       |
| ASASSN        |    Find light curves in the full variable star catalogue     |    Requires download of the full catalogue (>40 GB), off by default     |
| ASASSN Web        |    Autogenerate web search criteria for forced photometry   |           |
| Catalina/CRTS        |    All epoch photometry   |           |
| CDS        |    Obtain photometric SED for any search radius. Clickable link to CDS for the object   |           |
| Gaia        |    All epoch photometry/spectra/RVS   |           |
| Kepler        |    Query the K2 field, extract lightcurves from tpf files    |       Only query if there are no close contaminants   (in progress, I recommend using your own scripts) |
| NEOWISE        |    Query all photometry with RA/Dec entries around the object   |          The processor to WISE |
| Panstarrs        |    All epoch photometry from DR1   |           |
| PTF        |    All epoch photometry   |          (iPTF to be included) |
| SDSS        |    Spectra, finding charts, mean star magnitudes   |           |
| TESS        |    Query TESS, extract lightcurves from tpf files   |        Only query if there are no close contaminants  (in progress, I recommend using your own scripts) |
| WISE        |    Query all photometry with RA/Dec entries around the object   |           |
| ZTF        |    All epoch photometry   |           |


  


### Extras:  
Gaia HR plotter  
Plot all photometry  
Basic Lomb-Scargle  
checkLocalStars.py - performs a search for Gaia to not bother obtaining TESS/ATLAS/Kepler photometry under certain search radius conditions (read it)  




# What do you need to do?
- Clone this repository and look at getScripts/getpwds.py  . Here you need to include your [ATLAS](https://fallingstar-data.com/forcedphot/) login details and your [IRSA](https://irsa.ipac.caltech.edu/Missions/ztf.html) account details.
- Open the obtainDataMutiprocess.py script and scroll to the __main__ part.  You will see boolean statements for which surveys you want to query and some optional commands (e.g. proper motion) which need to be set... a list of default options are included but may not grant full functionality.
- Read in your RA and Dec (or lists of these)
- (Optional) Inspect the search options for each survey - separated at the top of obtainMultiProcess.py are the search radii for each survey
- If you simply run the script after unpacking the

# Where might the queries break?
- Is your proper motion huge? Not all surveys have a proper motion option that can be accounted for. You should increase your search radius... but be prepared for some junk
- Remember that ground based telescopes will have their declination limits - your source might not be available everywhere
 

## Required (non standard) packages:

| Package      | Tested with version  |
| :---        |    :----   |
|[Multiprocess](https://pypi.org/project/multiprocess/)|  
|[Astropy](https://docs.astropy.org/en/stable/install.html)| 5.1 |   
|[Beautiful Soup](https://pypi.org/project/beautifulsoup4/)| 4.11.1 |  
|[PIL](https://pypi.org/project/Pillow/)| 9.2.0 | 
|[Lightkurve](https://docs.lightkurve.org/about/install.html)| 2.2.0 |  
|[Astroquery](https://astroquery.readthedocs.io/en/latest/)|  0.4.6 |  
|[jdcal](https://pypi.org/project/jdcal/)| 1.4.1 |  


# Missing a source that you want included?
Push a script to access the data and/or supply me basic details of how the data is obtained.




