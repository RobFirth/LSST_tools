"""

"""
from __future__ import print_function

try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

from matplotlib import pyplot as plt
import matplotlib.colors as mpl_colors

import numpy as np
import pandas as pd

from sqlalchemy import create_engine
from astropy.cosmology import LambdaCDM
from scipy.interpolate import InterpolatedUnivariateSpline

import lsst_tools.utils as utils
import lsst_tools.utils.colours as colours

import pyCoCo as pccsims
import pycoco as pcc


def connect_to_db(opsimdbpath="/Users/berto/data/LSST/OpSimOutputDBs/minion_1016_sqlite.db"):
    """

    """
    # opsimdbpath = "/Users/berto/data/LSST/OpSimOutputDBs/minion_1016_sqlite.db"
    # opsimdbpath = "/Users/berto/data/LSST/OpSimOutputDBs/astro_lsst_01_1068_sqlite.db"
    # opsimdbpath = "/Users/berto/data/LSST/OpSimOutputDBs/Fake_Rolling/Rolling_3_80/Rolling_3_80.db"

    conn = create_engine('sqlite:///'+opsimdbpath, echo=False)
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


def generate_lc(opsimdf):
    """

    """

    n = 1
    pos_df = utils.generate_coordinates(n)

    allfields_df = utils.get_field_corners(pd.DataFrame(working_df.drop_duplicates("fieldID"))[['fieldID', 'fieldRA', 'fieldDec']])

    field_df = allfields_df[(allfields_df["RA_upper"] >= pos_df["theta"][0] + np.pi) &
                (allfields_df["RA_lower"] <= pos_df["theta"][0] + np.pi) &
                (allfields_df["Dec_lower"] >= pos_df["phi"][0])  &
                (allfields_df["Dec_upper"] <= pos_df["phi"][0])]

    pos_df = utils.find_MW_extinction(pos_df)

    extinction = pos_df["EBV_MW"].values[0]
    print(extinction)

    # filter_path = "/Users/berto/Code/CoCo/data/filters"
    filter_path = pcc._default_filter_dir_path
    # coco_root_path = "/Users/berto/Code/CoCo"
    coco_root_path = pcc._default_coco_dir_path

    coco = pccsims.pyCoCo(utils.b(filter_path), utils.b(coco_root_path))

    z_obs = 0.007
    # z_obs = 0.1
    # z_obs = 0.2
    host_EBV = 0.2
    MW_EBV = extinction
    # mjdmax = 60307.314753999999
    # mjdmax = 59580 + 1.* 365. ## DDF 2786
    # mjdmax = 59580 + 1.5* 365. ## WFD 550
    mjdmax = 59580 + 1.3* 365. ## WFD 2297

    mjd_to_sim = working_df["expMJD"][working_df["fieldID"].isin(field_df["fieldID"].values)].values
    limiting_mags = working_df["fiveSigmaDepth"][working_df["fieldID"].isin(field_df["fieldID"].values)].values

    filters_to_sim = working_df["filter"][working_df["fieldID"].isin(field_df["fieldID"].values)].values
    filters_to_sim = np.array([utils.b('LSST_'+x) for x in filters_to_sim])

    snname = b"SN2007uy"
    # snname = b"SN2009jf"
    flux, flux_err = coco.simulate(snname,
                        z_obs, 0.0, MW_EBV, host_EBV, 3.1,
                        mjdmax, mjd_to_sim,
                        filters_to_sim)

    SNp = pcc.PhotometryClass()
    phot_table = pcc.utils.simulate_out_to_ap_table(mjd_to_sim, flux, flux_err, filters_to_sim)
    phot_table = phot_table[np.where(phot_table["flux"] > 1e-20)]
    SNp.load_table(phot_table)

    return SNp

# def main1():
#     opsimdf= connect_to_db()
#
#     with closing(Pool(processes = cpu_count())) as pool:
#         x = pool.map(partial(visualise_sky_coverage_multiproc, ), )
#         pool.terminate()
#     pass
#
#     pass
#
# def visualise_sky_coverage_multiproc(df, ddf_fieldID):
#     df_to_plot = df[["fieldID", "fieldRA", "fieldDec", "cmap_visit_value"]]
#
#     # outpath = "/Users/berto/plots/LSST/cadence/minion/Test_skymap_night_"
#     outpath = os.path.join(os.environ["HOME"], "plots/LSST/cadence/Fake_Rolling/skymap_night_")
#
#     df_to_plot["zorder"] = df_to_plot["fieldID"]*0
#     ###
#     # DDF fieldID???
#     ddf = [1427,744,2412,290,2786] ## for opsimdb astro_lsst_01_1068, minion_1016
#     # ddf = [290,2786]
#     ###
#     df_to_plot["zorder"].loc[df_to_plot["fieldID"].isin(ddf)] = 99
#     norm = mpl_colors.LogNorm(vmin=max_visits, vmax=min_visits)
#
#     for i, group in enumerate(night_group.groups):
#
#         if i == 0:
#             fig = plt.figure()
#             fig.subplots_adjust(left = 0.05, bottom = 0.05, top = 0.99,
#                                 right = 0.97, hspace=0, wspace = .1)
#
#             ax_aitoff = fig.add_subplot(111, projection = "aitoff")
#             ax_active = fig.add_subplot(111, projection = "aitoff")
#             ax_aitoff.grid(True)
#             ax_active.grid(True)
#
#             ax_aitoff.set_xticklabels(['$2^h$','$4^h$','$6^h$', '$8^h$', '$10^h$', '$12^h$',
#                                       '$14^h$','$16^h$','$18^h$','$20^h$','$22^h$'], fontsize=10)
#             ax_active.set_xticklabels(['$2^h$','$4^h$','$6^h$', '$8^h$', '$10^h$', '$12^h$',
#                                       '$14^h$','$16^h$','$18^h$','$20^h$','$22^h$'], fontsize=10)
#
#             ax_cbar = fig.add_subplot(20,1,19)
#             cb1 = mpl_colorbar.ColorbarBase(ax_cbar, cmap=cmap,
#                                     norm=norm, orientation='horizontal')
#             cb1.set_label(r'$\textnormal{Number of Visits}$')
#
#             working_df_to_plot = df_to_plot.loc[night_group.groups[group]]
#
#         else:
#             working_df_to_plot = pd.concat([working_df_to_plot,df_to_plot.loc[night_group.groups[group]]])
#
#             working_df_to_plot.drop_duplicates(subset=['fieldID'], inplace=True, keep="last")
#
#         s = r"$\textnormal{night} = " + str(i) + "$"
#         txt = plt.figtext(0.85, 0.8, s)
#
#         s1 = ax_aitoff.scatter(working_df_to_plot["fieldRA"] - np.pi,
#                           working_df_to_plot["fieldDec"],
#                           color = cmap(working_df_to_plot["cmap_visit_value"]))
#
#
#         s2 = ax_aitoff.scatter(working_df_to_plot.loc[df_to_plot["zorder"] == 99,["fieldRA"]] - np.pi,
#                           working_df_to_plot.loc[df_to_plot["zorder"] == 99,["fieldDec"]],
#                           color = cmap(working_df_to_plot[working_df_to_plot["zorder"] == 99]["cmap_visit_value"])
#                           )
#
#         s_a = ax_active.scatter(df_to_plot.loc[night_group.groups[group]]["fieldRA"] - np.pi,
#                           df_to_plot.loc[night_group.groups[group]]["fieldDec"],
#                           color = cmap(df_to_plot.loc[night_group.groups[group]]["cmap_visit_value"]),
#                           edgecolor = colours.hex["black"])
#
#
#         plt.draw()
#         fig.savefig(outpath + str('%05d' % i)+".png", format = 'png', dpi=200)
#         s_a.remove()
#         txt.remove()
# working_df = utils.get_field_corners(working_df.drop_duplicates("fieldID"))[['fieldRA', 'fieldDec']]


def generate_lc(working_df, pos_df = False, n = 1):

    if not pos_df:
        pos_df = utils.generate_coordinates(n)

    fields_df = utils.get_field_corners(pd.DataFrame(working_df.drop_duplicates("fieldID"))[['fieldID', 'fieldRA', 'fieldDec']])

    field_df = fields_df[(fields_df["RA_upper"] >= df["theta"][0] + np.pi) &
                (fields_df["RA_lower"] <= df["theta"][0] + np.pi) &
                (fields_df["Dec_lower"] >= df["phi"][0])  &
                (fields_df["Dec_upper"] <= df["phi"][0])]

    coco = pccsims.pyCoCo(b(filter_path), b(coco_root_path))

    return


#  # Originally used for SDSS Sample Generation

def choose_subtype(return_string=True):
    """
    1|Ib | 11|0.328|
    2|Ic | 9.7|0.290|
    3|IIb | 12.8|0.382|
    total - 33.5
    """
    n = 33.5 * np.random.random()

    if n <= 11:
        if return_string:
            return "Ib"
        else:
            return 1

    if 11 < n <= 11 + 9.7:
        if return_string:
            return "Ic"
        else:
            return 2

    if 11 + 9.7 < n <= 11 + 9.7 + 12.8:
        if return_string:
            return "IIb"
        else:
            return 3
    pass


def choose_extinction_host(x0=0.0, sigma=0.2, n=10000):
    """
    """
    return np.fabs(np.random.normal(loc=x0, scale=sigma, size=n))


def choose_magoffset(x0=0.0, sigma=0.2, n=10000):
    """
    """
    return np.random.normal(loc=x0, scale=sigma, size=n)


def calculate_SFR(z):
    """

    :param z:
    :return:
    """
    a = 0.0170
    b = 0.13
    c = 3.3
    d = 5.3
    h = 0.7

    sfr = (a + b * z) * h / (1. + (z / c) ** d)

    return sfr


def string_format_for_mapping(x, prefix = "SDSS_"):
    return pcc.utils.b(prefix + x)


def choose_z_flat(z_max=0.6, n=1):
    """

    :param n:
    :param z_max:
    :return:
    """
    return z_max * np.random.random(n)


def choose_z_volume_SFR(n_req=10000, zmax=8.0, binsize=0.01):
    """

    :return:
    """


    z = np.arange(0.0, zmax, binsize)
    z_dz = np.arange(0.0 + binsize, zmax + binsize, binsize)

    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

    v_z = cosmo.comoving_volume(z)
    v_z_dz = cosmo.comoving_volume(z_dz)

    v_dz = v_z_dz - v_z

    norm_v_dz = v_dz / np.nanmax(v_dz)

    sfr_z = calculate_SFR(z)
    sfr_norm = sfr_z / np.nanmax(sfr_z)

    volumetric_rate = norm_v_dz * sfr_norm
    normed_volumetric_rate = volumetric_rate / np.nanmax(volumetric_rate)

    pdf = InterpolatedUnivariateSpline(z, normed_volumetric_rate)


    n = 0
    z_sim = []
    while n < n_req:
        x = np.random.random() * zmax
        y = np.random.random()

        if y <= pdf(x):
            z_sim.append(x)
            n += 1


    return z_sim


def choose_MJDmax(obslog, n=1):
    """

    :param obslog:
    :param n:
    :return:
    """

    random_MJD = (obslog.mjd.max() - obslog.mjd.min()) * np.random.random(n) + obslog.mjd.min()
    return random_MJD
