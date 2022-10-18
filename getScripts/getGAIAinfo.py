#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 14 12:46:43 2022

@author: james
"""



https://www.cosmos.esa.int/documents/29201/5432054/DraftDataModel-DR3.pdf/d2049ffb-2cea-2581-4b66-2ea859686018?t=1648466925132
page 315



7.1 epoch
_photometry
Epoch photometry. Each row in this table contains the light curve for
a given object in bands G, BP and RP as stored in the DataLink Massive data base. This table makes extensive use of array types. It can be
obtained selecting the RAW data structure option. A flat table (sparse
cube), which one photometric point per source per row can be obtained
using INDIVIDUAL or COMBINED.
Note this table is not available through the main archive TAP interface,
but via the Massive Data service, indexed by the VO Datalink protocol,
described in Sect. ??
.







##
there is a flag "HAS_EPOCH_RV" to tell if a target has them or not
also "vari_epoch_radial_velocity"