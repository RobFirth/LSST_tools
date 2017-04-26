"""

"""
from __future__ import print_function

try:
    from importlib import reload
except:
    pass


from matplotlib import pyplot as plt
import matplotlib.colors as mpl_colors

import os
import warnings

import numpy as np
import pandas as pd
import astropy as ap
import sfdmap

import sqlite3
from sqlalchemy import create_engine

from astropy import units as u
from astropy.coordinates import SkyCoord

import lsst_tools.utils as utils
import lsst_tools.utils.colours as colours

import pyCoCo as pccsims

from multiprocessing import Pool, cpu_count
from functools import partial

def connect_to_db(opsimdbpath = "/Users/berto/data/LSST/OpSimOutputDBs/minion_1016_sqlite.db"):
    """

    """
    # opsimdbpath = "/Users/berto/data/LSST/OpSimOutputDBs/minion_1016_sqlite.db"
    # opsimdbpath = "/Users/berto/data/LSST/OpSimOutputDBs/astro_lsst_01_1068_sqlite.db"
    # opsimdbpath = "/Users/berto/data/LSST/OpSimOutputDBs/Fake_Rolling/Rolling_3_80/Rolling_3_80.db"

    conn = create_engine('sqlite:///'+opsimdbpath, echo = False)
    opsimdf = pd.read_sql_table('Summary', con=conn)
    return opsimdf

def visualise_sky_coverage(df, ddf_fieldID):
    df_to_plot = df[["fieldID", "fieldRA", "fieldDec", "cmap_visit_value"]]

    # outpath = "/Users/berto/plots/LSST/cadence/minion/Test_skymap_night_"
    outpath = os.path.join(os.environ["HOME"], "plots/LSST/cadence/Fake_Rolling/skymap_night_")

    df_to_plot["zorder"] = df_to_plot["fieldID"]*0
    ###
    # DDF fieldID???
    ddf = [1427,744,2412,290,2786] ## for opsimdb astro_lsst_01_1068, minion_1016
    # ddf = [290,2786]
    ###
    df_to_plot["zorder"].loc[df_to_plot["fieldID"].isin(ddf)] = 99
    norm = mpl_colors.LogNorm(vmin=max_visits, vmax=min_visits)

    for i, group in enumerate(night_group.groups):

        if i == 0:
            fig = plt.figure()
            fig.subplots_adjust(left = 0.05, bottom = 0.05, top = 0.99,
                                right = 0.97, hspace=0, wspace = .1)

            ax_aitoff = fig.add_subplot(111, projection = "aitoff")
            ax_active = fig.add_subplot(111, projection = "aitoff")
            ax_aitoff.grid(True)
            ax_active.grid(True)

            ax_aitoff.set_xticklabels(['$2^h$','$4^h$','$6^h$', '$8^h$', '$10^h$', '$12^h$',
                                      '$14^h$','$16^h$','$18^h$','$20^h$','$22^h$'], fontsize=10)
            ax_active.set_xticklabels(['$2^h$','$4^h$','$6^h$', '$8^h$', '$10^h$', '$12^h$',
                                      '$14^h$','$16^h$','$18^h$','$20^h$','$22^h$'], fontsize=10)

            ax_cbar = fig.add_subplot(20,1,19)
            cb1 = mpl_colorbar.ColorbarBase(ax_cbar, cmap=cmap,
                                    norm=norm, orientation='horizontal')
            cb1.set_label(r'$\textnormal{Number of Visits}$')

            working_df_to_plot = df_to_plot.loc[night_group.groups[group]]

        else:
            working_df_to_plot = pd.concat([working_df_to_plot,df_to_plot.loc[night_group.groups[group]]])

            working_df_to_plot.drop_duplicates(subset=['fieldID'], inplace=True, keep="last")

        s = r"$\textnormal{night} = " + str(i) + "$"
        txt = plt.figtext(0.85, 0.8, s)

        s1 = ax_aitoff.scatter(working_df_to_plot["fieldRA"] - np.pi,
                          working_df_to_plot["fieldDec"],
                          color = cmap(working_df_to_plot["cmap_visit_value"]))


        s2 = ax_aitoff.scatter(working_df_to_plot.loc[df_to_plot["zorder"] == 99,["fieldRA"]] - np.pi,
                          working_df_to_plot.loc[df_to_plot["zorder"] == 99,["fieldDec"]],
                          color = cmap(working_df_to_plot[working_df_to_plot["zorder"] == 99]["cmap_visit_value"])
                          )

        s_a = ax_active.scatter(df_to_plot.loc[night_group.groups[group]]["fieldRA"] - np.pi,
                          df_to_plot.loc[night_group.groups[group]]["fieldDec"],
                          color = cmap(df_to_plot.loc[night_group.groups[group]]["cmap_visit_value"]),
                          edgecolor = colours.hex["black"])


        plt.draw()
        fig.savefig(outpath + str('%05d' % i)+".png", format = 'png', dpi=200)
        s_a.remove()
        txt.remove()

def main1():
    opsimdf= connect_to_db()

    with closing(Pool(processes = cpu_count())) as pool:
        x = pool.map(partial(visualise_sky_coverage_multiproc, ), )
        pool.terminate()
    pass

    pass

def visualise_sky_coverage_multiproc(df, ddf_fieldID):
    df_to_plot = df[["fieldID", "fieldRA", "fieldDec", "cmap_visit_value"]]

    # outpath = "/Users/berto/plots/LSST/cadence/minion/Test_skymap_night_"
    outpath = os.path.join(os.environ["HOME"], "plots/LSST/cadence/Fake_Rolling/skymap_night_")

    df_to_plot["zorder"] = df_to_plot["fieldID"]*0
    ###
    # DDF fieldID???
    ddf = [1427,744,2412,290,2786] ## for opsimdb astro_lsst_01_1068, minion_1016
    # ddf = [290,2786]
    ###
    df_to_plot["zorder"].loc[df_to_plot["fieldID"].isin(ddf)] = 99
    norm = mpl_colors.LogNorm(vmin=max_visits, vmax=min_visits)

    for i, group in enumerate(night_group.groups):

        if i == 0:
            fig = plt.figure()
            fig.subplots_adjust(left = 0.05, bottom = 0.05, top = 0.99,
                                right = 0.97, hspace=0, wspace = .1)

            ax_aitoff = fig.add_subplot(111, projection = "aitoff")
            ax_active = fig.add_subplot(111, projection = "aitoff")
            ax_aitoff.grid(True)
            ax_active.grid(True)

            ax_aitoff.set_xticklabels(['$2^h$','$4^h$','$6^h$', '$8^h$', '$10^h$', '$12^h$',
                                      '$14^h$','$16^h$','$18^h$','$20^h$','$22^h$'], fontsize=10)
            ax_active.set_xticklabels(['$2^h$','$4^h$','$6^h$', '$8^h$', '$10^h$', '$12^h$',
                                      '$14^h$','$16^h$','$18^h$','$20^h$','$22^h$'], fontsize=10)

            ax_cbar = fig.add_subplot(20,1,19)
            cb1 = mpl_colorbar.ColorbarBase(ax_cbar, cmap=cmap,
                                    norm=norm, orientation='horizontal')
            cb1.set_label(r'$\textnormal{Number of Visits}$')

            working_df_to_plot = df_to_plot.loc[night_group.groups[group]]

        else:
            working_df_to_plot = pd.concat([working_df_to_plot,df_to_plot.loc[night_group.groups[group]]])

            working_df_to_plot.drop_duplicates(subset=['fieldID'], inplace=True, keep="last")

        s = r"$\textnormal{night} = " + str(i) + "$"
        txt = plt.figtext(0.85, 0.8, s)

        s1 = ax_aitoff.scatter(working_df_to_plot["fieldRA"] - np.pi,
                          working_df_to_plot["fieldDec"],
                          color = cmap(working_df_to_plot["cmap_visit_value"]))


        s2 = ax_aitoff.scatter(working_df_to_plot.loc[df_to_plot["zorder"] == 99,["fieldRA"]] - np.pi,
                          working_df_to_plot.loc[df_to_plot["zorder"] == 99,["fieldDec"]],
                          color = cmap(working_df_to_plot[working_df_to_plot["zorder"] == 99]["cmap_visit_value"])
                          )

        s_a = ax_active.scatter(df_to_plot.loc[night_group.groups[group]]["fieldRA"] - np.pi,
                          df_to_plot.loc[night_group.groups[group]]["fieldDec"],
                          color = cmap(df_to_plot.loc[night_group.groups[group]]["cmap_visit_value"]),
                          edgecolor = colours.hex["black"])


        plt.draw()
        fig.savefig(outpath + str('%05d' % i)+".png", format = 'png', dpi=200)
        s_a.remove()
        txt.remove()
