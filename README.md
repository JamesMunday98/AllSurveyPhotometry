# All Survey Time-Series Photometry
The package will obtain all time-series photometry associated with the object that is publicly available. All you need is the RA and Dec!

  
# Current list that can be accessed with the ease:  

| Survey      | Function  | Comments     |  Output Time Format (UTC, unless special)  |
| :---        |    :----   |    :----       | :----       |
| ATLAS        |    Forced photometry from the ATLAS survey     |    Only query if there are no close contaminants       | MJD | 
| ASASSN        |    Find light curves in the full variable star catalogue     |    Requires download of the full catalogue (>40 GB), off by default  | BJD  |
| ASASSN Web        |    Autogenerate web search criteria for forced photometry   |           |
| Catalina/CRTS        |    All epoch photometry   | |    BJD       |
| CDS        |    Obtain photometric SED for any search radius. Clickable link to CDS for the object   |           |
| Gaia        |    All epoch photometry/spectra/RVS   |    Gaia units... BJD but slightly different       |
| Kepler        |    Query the K2 field, extract lightcurves from tpf files    |       Only query if there are no close contaminants   (in progress, I recommend using your own scripts) | Units from Kepler |
| NEOWISE        |    Query all photometry with RA/Dec entries around the object   |          The processor to WISE |  BJD |
| Panstarrs        |    All epoch photometry from DR1   | |     BJD      |
| PTF        |    All epoch photometry   |          (iPTF to be included) |  MJD   |
| SDSS        |    Spectra, finding charts, mean star magnitudes   |           |
| TESS        |    Query TESS, extract lightcurves from tpf files   |        Only query if there are no close contaminants  (in progress, I recommend using your own scripts) | BJD | 
| WISE        |    Query all photometry with RA/Dec entries around the object  | |    BJD       |
| ZTF        |    All epoch photometry  | |    MJD       |


  


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


# A note on science usage
I recommend using this package to simply to check what data is available on various surveys... and then to check the data outputs and use my individual getScripts to help with that. I take no responsibility for supplying data for published work - for this you should check the original sources.

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



# To-do list
Make all time similar (will not matter much at all for most science needs)  
Include ASAS  
Include APASS
Improve TESS/Kepler


