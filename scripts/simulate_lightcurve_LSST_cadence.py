
# coding: utf-8

# # Get a LSST Simulated Cadence for Arbitrary Sky Position
# ___
#
#
# ___

# In[ ]:

"""

"""
get_ipython().magic('matplotlib inline')
# %matplotlib notebook

from __future__ import print_function

try:
    from importlib import reload
except:
    pass


from matplotlib import pyplot as plt
import matplotlib.colors as mpl_colors

import os
import sys
import warnings

import numpy as np
import pandas as pd
import astropy as ap
import sfdmap

import sqlite3
from sqlalchemy import create_engine

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

import lsst_tools.utils as utils
import lsst_tools.utils.colours as colours

import pyCoCo as pccsims
import pycoco as pcc


# Connect to .db file that contains the opsim output, read into a dataframe, `opsimdf`

# In[ ]:

# %%timeit ## 1 loop, best of 3: 1min 31s per loop
opsimdbpath = os.environ.get('OPSIMDBPATH')
print(opsimdbpath)
# opsimdbpath = "/Users/berto/data/LSST/OpSimOutputDBs/astro_lsst_01_1068_sqlite.db"
opsimdbpath = "/Users/berto/data/LSST/OpSimOutputDBs/minion_1016_sqlite.db"

conn = create_engine('sqlite:///'+opsimdbpath, echo = False)
opsimdf = pd.read_sql_table('Summary', con=conn)

print("OPSIMDF read")
# Check that the db looks as we expect

# In[ ]:

# opsimdf.head()


# Connecting to `.db` takes ages (~1min), and is a pain if you mess up, so create a 'working' instance to fiddle with

# In[ ]:

working_df = opsimdf


# In[ ]:

working_df = utils.find_MW_extinction(working_df)


# ## Choosing Position and Identifying Fields
# ___
#
# Use **`lsst_tools.utils.generate_coordinates()`** to give us a position
#

# In[ ]:

#%%timeit ## The slowest run took 39.04 times longer than the fastest. This could mean that an intermediate result is being cached.
#         ## 1000 loops, best of 3: 246 µs per loop
n = 1
pos_df = utils.generate_coordinates(n)
print(pos_df)
# pos_df["phi"] = -0.122
# pos_df["theta"] = 0.0


# Check the positions are sensible

# In[ ]:

#%%timeit ## 1 loop, best of 3: 235 ms per loop
utils.plot_position_points(pos_df)


# We want to find out which field(s) the position is in. Create table that tells us the field centres, and find the edges.

# In[ ]:

# working_df = utils.get_field_corners(working_df.drop_duplicates("fieldID"))[['fieldRA', 'fieldDec']]
allfields_df = utils.get_field_corners(pd.DataFrame(working_df.drop_duplicates("fieldID"))[['fieldID', 'fieldRA', 'fieldDec']])


# In[ ]:

# allfields_df.head()


# In[ ]:

field_df = allfields_df[(allfields_df["RA_upper"] >= pos_df["theta"][0] + np.pi) &
                (allfields_df["RA_lower"] <= pos_df["theta"][0] + np.pi) &
                (allfields_df["Dec_lower"] >= pos_df["phi"][0])  &
                (allfields_df["Dec_upper"] <= pos_df["phi"][0])]


# This narrows down the fields in which our position appears:

# In[ ]:

# field_df


# ### Deep Drilling Fields:
# ___
#
# * minion2016: 1427,744,2412,290,2786
# * astro_lsst_01_1068: 1427,744,2412,290,2786
# * Fake_Rolling: 290,2786
#

# In[ ]:

# opsimdf.loc[opsimdf["fieldID"] == 505].head() ## WFD
#
#
# # In[ ]:
#
# opsimdf.loc[opsimdf["fieldID"] == 1427].head() ## DDF


#

# In[ ]:




# In[ ]:




# ## Extinction

# Get the Milky Way extinction along the line of site towards the SNe. The working_df contains the EBV_MW at the field centre, but we can do better than that, by using the position of the SN itself.

# In[ ]:

working_df["EBV_MW"][working_df["fieldID"].isin(working_df["fieldID"].values)].values


# In[ ]:

# pos_df


# In[ ]:

# %%timeit ## 1 loop, best of 3: 88.6 ms per loop
pos_df = utils.find_MW_extinction(pos_df)

extinction = pos_df["EBV_MW"].values[0]
print("EBV_MW = ", extinction)


# In[ ]:

# filter_path = "/Users/berto/Code/CoCo/data/filters"
filter_path = pcc._default_filter_dir_path
# coco_root_path = "/Users/berto/Code/CoCo"
coco_root_path = pcc._default_coco_dir_path

coco = pccsims.pyCoCo(utils.b(filter_path), utils.b(coco_root_path))


# inputs:
#
# * **`snname`**
# * **`redshift`**
# * **`absmag offset`**
# * **`EBV MW`**
# * **`EBV Host`**
# * **`Rv`**
# * **`MJD at Peak`**
#
# * **`MJD to simulate`**
# * **`filters to simulate`**

# In[ ]:

z_obs = 0.007
# z_obs = 0.1
# z_obs = 0.2
host_EBV = 0.1
MW_EBV = extinction
# mjdmax = 60307.314753999999
# mjdmax = 59580 + 1.* 365. ## DDF 2786
# mjdmax = 59580 + 1.5* 365. ## WFD 550
mjdmax = 59580 + 1.3* 365. ## WFD 2297


# In[ ]:

mjd_to_sim = working_df["expMJD"][working_df["fieldID"].isin(field_df["fieldID"].values)].values
limiting_mags = working_df["fiveSigmaDepth"][working_df["fieldID"].isin(field_df["fieldID"].values)].values

filters_to_sim = working_df["filter"][working_df["fieldID"].isin(field_df["fieldID"].values)].values
filters_to_sim = np.array([utils.b('LSST_'+x) for x in filters_to_sim])


# In[ ]:

# mjd_to_sim = np.array([90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 105.0, 110.0, 115.0, 120.0])
# filters_to_sim = np.array([b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r', b'SDSS_r'])
# mjd_to_sim =
# filters_to_sim =

print(mjdmax)
print(mjd_to_sim)
print(filters_to_sim)


# In[ ]:

# plt.scatter(mjd_to_sim, np.ones(len(mjd_to_sim)))
# plt.scatter([mjdmax], [1])


# In[ ]:

# filters_to_sim = np.array([i.replace(utils.b("LSST"), utils.b("SDSS")) for i in filters_to_sim])


# In[ ]:

# flux, flux_err = coco.simulate(b"SN2007uy",
#                     z_obs, 0.0, MW_EBV, host_EBV, 3.1,
#                     mjdmax, mjd_to_sim,
#                     filters_to_sim)
snname = b"SN2007uy"
# snname = b"SN2009jf"
flux, flux_err = coco.simulate(snname,
                    z_obs, 0.0, MW_EBV, host_EBV, 3.1,
                    mjdmax, mjd_to_sim,
                    filters_to_sim)
# flux, flux_err = coco.simulate(b"SN2009jf",
#                     z_obs, 0.0, 0.1, 0.1, 3.1,
#                     mjdmax, mjd_to_sim,
#                     filters_to_sim)


# In[ ]:

# flux, flux_err


# In[ ]:

# print(np.nanmax(flux))


# In[ ]:

# reload(pcc.classes)


# In[ ]:

p = pcc.PhotometryClass()
phot_table = pcc.utils.simulate_out_to_ap_table(mjd_to_sim, flux, flux_err, filters_to_sim)
phot_table = phot_table[np.where(phot_table["flux"] > 1e-20)]
p.load_table(phot_table)


# In[ ]:

p.plot(xlim = [59950, 60100])


# In[ ]:

# p.data


# In[ ]:

# t = pcc.utils.simulate_out_to_ap_table(mjd_to_sim, flux, flux_err, filters_to_sim)


# In[ ]:

# t[np.logical_and(t["filter"] != "SDSS_y", t["filter"]!= "SDSS_z")]


# In[ ]:
#
# p = pcc.PhotometryClass()
# pt = pcc.utils.simulate_out_to_ap_table(mjd_to_sim, flux, flux_err, filters_to_sim)
# pt_lim_m = pt
# pt_lim_m["fivesigmadepth"] = limiting_mags
#
#
# # In[ ]:
#
# pt = pt_lim_m[np.where(pt_lim_m["flux"] > 1e-22)]
#
#
# # In[ ]:
#
# pt
#
#
# # In[ ]:
#
# days_before = 40
# days_after = 100
# print(mjdmax)
# w = np.where(np.logical_and(pt["MJD"] > mjdmax - days_before, pt["MJD"] < mjdmax + days_after))
# print(w)
# pt[w]
#
#
# # In[ ]:
#
# print(snname)
# print(chosenfield)
# print(z_obs)
# print(extinction)
# print(opsimdbpath)
# print(mjdmax)
#
#
# # In[ ]:
#
# reload(pcc)
# reload(pcc.classes)
# # reload(pcc.utils)
#
#
# # In[ ]:
#
#
#
#
# # In[ ]:
#
# p.load_table(pt[w], verbose= True)
#
#
# p.plot()
#
#
# # In[ ]:
#
# p.unpack(verbose = True)
#
#
# # In[ ]:
#
# phot_table = p.phot.loc["filter", "LSST_u"]
#
#
# # In[ ]:
#
# phot_table
#
#
# # In[ ]:
#
# phot_table.meta["filter_filename"] = "foo"
#
#
# # In[ ]:
#
# pt[w]


# In[ ]:

# filter_file_type = '.dat'
# filter_names = np.unique(pt[w]["filter"])
# print(filter_names)
# pt.add_index('filter', unique = True)
# for filter_name in filter_names:
# #     phot_table = self.phot.loc["filter", filter_name]
#     filter_filename = filter_name + filter_file_type
#     print(filter_filename)


# In[ ]:

# reload(pcc)
# phot = pcc.PhotometryClass()
# infile = "/Users/berto/projects/LSST/cadence/lightcurves/SNSim_0002_minion_1016_SN2007uy_z=02_EBVMW=0029_EBVHOST=01_fieldID=2297WFD.dat"
# phot.load(path = infile, names = ("MJD", "flux", "flux_err", "filter", "fivesigmadepth"), verbose = True)
# # phot.load_phot_from_file(path = infile, format = "ascii.commented_header")

# phot.unpack(verbose = True)


# In[ ]:

# phot.data


# In[ ]:

# phot_table = Table.read(infile, names = ("MJD", "flux", "flux_err", "filter", "fivesigmadepth"), format = "ascii")


# In[ ]:

# phot_table.meta["filename"] = infile

# phot_table["MJD"].unit = u.day
# phot_table["flux"].unit = u.cgs.erg / u.si.angstrom / u.si.cm ** 2 / u.si.s
# phot_table["flux_err"].unit =  phot_table["flux"].unit


# In[ ]:

# reload(pcc)
# p = pcc.PhotometryClass()
# p.load_table(phot_table, verbose = True)


# In[ ]:
