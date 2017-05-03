"""
Tools for looking at OpSim run outputs

Author: Rob Firth, University of Southampton, 03/2017
"""

from __future__ import print_function ## Force python3-like printing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
import os
import time

from sqlalchemy import create_engine

opsimdbpath = os.environ.get('OPSIMDBPATH')
