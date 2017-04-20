"""

"""

from __future__ import print_function

try:
    from importlib import reload
except:
    pass


from matplotlib import pyplot as plt

import os
import warnings

import numpy as np
import pandas as pd
import astropy as ap

from astropy import units as u
from astropy.coordinates import SkyCoord

import sfdmap



def generate_coordinates(nruns, low_ra = 0., high_ra = 360., low_dec = -90., high_dec = -10):
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
    return df


def find_MW_extinction(df):
    """

    Parameters
    ---

    Returns
    ---

    """

    m = sfdmap.SFDMap()
    df["EBV_MW"] = m.ebv(df["theta"], df["phi"], unit = "radian")
    return df

def plot_position_points(df):
    """
    for small realisations! If large, use plot_position_heatmap
    """
    fig = plt.figure()
    fig.subplots_adjust(left = 0.09, bottom = 0.13, top = 0.99,
                        right = 0.97, hspace=0, wspace = .1)

    ax_aitoff = fig.add_subplot(111, projection="aitoff")
    ax_aitoff.grid(True)

    ax_aitoff.scatter(df["theta"], df["phi"])

    plt.show()
    pass

def plot_position_heatmap(df):
    """

    """

    hist, xedges, yedges = np.histogram2d(df["theta"], df["phi"], bins = 100)
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
