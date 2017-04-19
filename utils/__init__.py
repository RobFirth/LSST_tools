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
