"""

"""

from __future__ import print_function

try:
    from importlib import reload
except:
    pass


from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches

import os
import sys
import warnings
import numpy as np
import pandas as pd
import astropy as ap

from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic, ICRS, Angle, Latitude, Longitude

from .colours import hex, RGB, RGB255

import sfdmap


def generate_coordinates(nruns, low_ra = 0., high_ra = 360., low_dec = -90., high_dec = 10):
    """
    Function to generate a set of coordinates for simulated SNe.

    Parameters
    ----------

    nruns     :

    low_ra    :

    high_ra   :

    low_dec   :

    high_dec  :


    Returns
    -------

    df     : pandas DataFrame with two columns - theta and phi, which are the
             spherical coordinates in radians

    -------

    RobFirth - with help from Wycombe7 and CFro

    """

    umin = np.radians(low_ra)/(2.*np.pi)
    umax = np.radians(high_ra)/(2.*np.pi)

    vmin = (np.cos(np.radians(90.+low_dec))+1.0)/2.
    vmax = (np.cos(np.radians(90.+high_dec))+1.0)/2.

    df = pd.DataFrame({"theta" : (2. * np.pi * (np.random.uniform(umin, umax, size = nruns) - 0.5)),
                      "phi" : (np.arccos(2. * np.random.uniform(vmin, vmax, size = nruns) - 1.) - np.pi/2.)})

    df["RA"] = df["theta"] + np.pi
    df["Dec"] = df["phi"]
    return df


def find_MW_extinction(df):
    """

    Parameters
    ---

    Returns
    ---

    """

    m = sfdmap.SFDMap()
    if "theta" in df.columns and "phi" in df.columns:
        df["EBV_MW"] = m.ebv(df["theta"].values, df["phi"].values, unit = "radian")
    elif "fieldRA" in df.columns and "fieldDec" in df.columns:
        df["EBV_MW"] = m.ebv(df["fieldRA"].values, df["fieldDec"].values, unit = "radian")
    return df


def plot_position_points(df, galacticplane = True, lw = 2 ):
    """
    for small realisations! If large, use plot_position_heatmap
    """
    fig = plt.figure()
    fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                        right = 0.97, hspace=0, wspace = .1)

    ax_aitoff = fig.add_subplot(111, projection="aitoff")
    ax_aitoff.grid(True)

    if "theta" in df.columns and "phi" in df.columns:
        ax_aitoff.scatter(df["theta"], df["phi"])
    elif "fieldRA" in df.columns and "fieldDec" in df.columns:
        ax_aitoff.scatter(df["fieldRA"] - np.pi, df["fieldDec"])
    else:
        ax_aitoff.scatter(df["RA"] - np.pi, df["Dec"])

    gp_df = get_galactic_plane()

    ax_aitoff.plot(gp_df["RA"], gp_df["Dec"], lw = lw, color = hex['wetasphalt'], alpha = 0.5)

    plt.show()
    pass


def plot_position_heatmap(df):
    """

    """

    if "theta" in df.columns and "phi" in df.columns:
        hist, xedges, yedges = np.histogram2d(df["theta"], df["phi"], bins = 100)
    else:
        hist, xedges, yedges = np.histogram2d(df["fieldRA"], df["fieldDec"], bins = 100)

    X, Y = np.meshgrid(xedges, yedges)

    # fig = plt.figure(figsize=[12, 4])
    fig = plt.figure()
    fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                        right = 0.97, hspace=0, wspace = .1)

    ax_aitoff = fig.add_subplot(111, projection = "aitoff")
    ax_aitoff.grid(True)

    # ax_flat = plt.subplot2grid((1,3), (0,0), colspan = 1, rowspan = 1)
    # ax_aitoff = plt.subplot2grid((1,3), (0,1), colspan = 2, rowspan = 1, projection="aitoff")

    # im_flat = ax_flat.pcolormesh(X, Y, hist.T)
    im_aitoff = ax_aitoff.pcolormesh(X, Y, hist.T)


    fig.colorbar(im_aitoff, ax=ax_aitoff)

    plt.show()
    pass


def get_field_corners(df, fov_degree = 3.5):
    """

    """
    fov = fov_degree*u.degree
    dist_to_edge = fov.to(u.radian).value/2.

    if "fieldRA" in df.keys() and "fieldDec" in df.keys():
        df["RA_upper"] = df["fieldRA"] + dist_to_edge
        df["RA_lower"] = df["fieldRA"] - dist_to_edge
        df["Dec_upper"] = df["fieldDec"] - dist_to_edge
        df["Dec_lower"] = df["fieldDec"] + dist_to_edge

    elif "RA" in df.keys() and "Dec" in df.keys():
        df["RA_upper"] = df["RA"] + dist_to_edge
        df["RA_lower"] = df["RA"] - dist_to_edge
        df["Dec_upper"] = df["Dec"] - dist_to_edge
        df["Dec_lower"] = df["Dec"] + dist_to_edge

    elif "phi" in df.keys() and "theta" in df.keys():
        df["RA_upper"] = df["theta"] + dist_to_edge - np.pi
        df["RA_lower"] = df["theta"] - dist_to_edge - np.pi
        df["Dec_upper"] = df["phi"] - dist_to_edge
        df["Dec_lower"] = df["phi"] + dist_to_edge

    return df


def plot_field(df):
    """
    for small realisations! If large, use plot_position_heatmap
    """


    fig = plt.figure()
    fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                        right = 0.97, hspace=0, wspace = .1)

    ax_aitoff = fig.add_subplot(111, projection="aitoff")
    ax_aitoff.grid(True)

    x = np.array([df.iloc[0,:]["RA_upper"], df.iloc[0,:]["RA_upper"], df.iloc[0,:]["RA_lower"], df.iloc[0,:]["RA_lower"], df.iloc[0,:]["RA_upper"]])
    y = np.array([df.iloc[0,:]["Dec_upper"], df.iloc[0,:]["Dec_lower"], df.iloc[0,:]["Dec_lower"], df.iloc[0,:]["Dec_upper"], df.iloc[0,:]["Dec_upper"]])

    ax_aitoff.plot(x - np.pi, y, lw =2)
    plt.show()
    pass


def plot_field_patch(df, fov_degree = 3.5, projection="hammer"):
    """
    for small realisations! If large, use plot_position_heatmap
    """
    fov = fov_degree*u.degree
    fov_rad = fov.to(u.radian)

    patches = []
    wRA = np.where(df.columns == "RA_lower")[0][0]
    wDec = np.where(df.columns == "Dec_lower")[0][0]

    for row in df.itertuples():
        rect = mpatches.Rectangle(xy = (row[wRA] - np.pi, row[wDec]),
                                  width = fov_rad.value, height = fov_rad.value, ec="Black")
        patches.append(rect)

    collection = PatchCollection(patches)

    fig = plt.figure()
    fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                        right = 0.97, hspace=0, wspace = .1)

    if not projection:
        ax = fig.add_subplot(111)
    else:
        ax = fig.add_subplot(111, projection = projection)

    ax.grid(True)
    ax.add_collection(collection)

    # x = np.array([df.iloc[0,:]["RA_upper"], df.iloc[0,:]["RA_upper"], df.iloc[0,:]["RA_lower"], df.iloc[0,:]["RA_lower"], df.iloc[0,:]["RA_upper"]])
    # y = np.array([df.iloc[0,:]["Dec_upper"], df.iloc[0,:]["Dec_lower"], df.iloc[0,:]["Dec_lower"], df.iloc[0,:]["Dec_upper"], df.iloc[0,:]["Dec_upper"]])
    #
    # ax_aitoff.plot(x - np.pi, y, lw =2)
    if not projection:
        ax.set_xlim(df["RA_lower"].min(skipna = True), df["RA_upper"].max(skipna = True))
        ax.set_ylim(df["Dec_lower"].min(skipna = True), df["Dec_upper"].max(skipna = True))

    plt.show()
    pass


def get_galactic_plane(remake = False, path = "data/galactic_plane_RADec.dat"):
    print(os.path.abspath(os.path.join(__file__, os.pardir, os.pardir)))
    fullpath = os.path.abspath(os.path.join(__file__, os.pardir, os.pardir, "data/galactic_plane_RADec.dat"))

    if remake:
        try:
            ## make line of Galactic = 0 Lat
            long = np.arange(-180., 181., 1.)
            lat = np.zeros_like(long)
            long_ra_rad_arr = np.array([])
            lat_dec_rad_arr = np.array([])

            for i in range(len(long)):
                sc = SkyCoord(long[i], lat[i], unit='deg', frame=Galactic)
                long_ra_rad_arr = np.append(long_ra_rad_arr, sc.icrs.ra.wrap_at(180*u.deg).radian)
                lat_dec_rad_arr = np.append(lat_dec_rad_arr, sc.icrs.dec.radian)

            w = np.argsort(long_ra_rad_arr)

            coord_table = ap.table.Table((long_ra_rad_arr[w], lat_dec_rad_arr[w]), names = ("RA", "Dec" ))

            coord_table.write(path = fullpath, format = "ascii.commented_header")
        except:
            warnings.warn("something went wrong")
    else:
        coord_table = ap.io.ascii.read(fullpath, format = "commented_header")
    return coord_table


def print_path():
    # fullpath = os.path.join(__file__, os.pardir, os.pardir)
    fullpath = os.path.join(__file__, os.pardir, os.pardir, "data/galactic_plane_RADec.dat")
    print(os.path.abspath(fullpath))
    pass


if sys.version_info < (3,):
    def b(x):
        return x
else:
    import codecs
    def b(x):
        return codecs.latin_1_encode(x)[0]

if __name__ == "__main__":

    print(os.pardir(__file__))

else:
    pass
