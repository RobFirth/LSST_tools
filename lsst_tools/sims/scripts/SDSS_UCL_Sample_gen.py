# import matplotlib.pyplot as plt
import json
import pandas as pd
import numpy as np

import pycoco as pcc
import pyCoCo as pccsim

from astropy.cosmology import LambdaCDM
from scipy.interpolate import InterpolatedUnivariateSpline

import matplotlib.pyplot as plt

import lsst_tools as lsstt
from lcsim.simlib import SIMLIBReader
from lcsim.lcsim import LCSim
from datetime import datetime

if __name__ == "__main__":
    print("Running __main__")
    # verbose = True
    verbose = False
    log = True
    # log = False
    # logall = False
    if log:
        gentime = str(datetime.now())
        logvars = ["gentime",
                   "logpath",
                   "field_index",
                   "field",
                   "CCD_index",
                   "z_sim",
                   "MW_EBV",
                   "mag_offset",
                   "host_EBV",
                   "mjdmax",
                   "subtype",
                   "w",
                   "snindex",
                   "snname",
                   "flux",
                   "n",
                   "n_sne",
                   ]

    n_sne_req = 10000
    z_max = 1.0
    plot = False

    ## Initialise pyCoCo
    fltPath = b"/Users/berto/Code/CoCo/data/filters"
    rootPath = b"/Users/berto/Code/CoCo"
    # rootPath = b"/Users/berto/projects/stunt_CoCo"

    coco = pccsim.pyCoCo(fltPath, rootPath)
    lcs = LCSim()

    simlib_file05 = "/Users/berto/projects/SDSS_sims/SDSS_SIMLIB/SDSS_2005.SIMLIB"
    simlib_file06= "/Users/berto/projects/SDSS_sims/SDSS_SIMLIB/SDSS_2006.SIMLIB"
    simlib_file07 = "/Users/berto/projects/SDSS_sims/SDSS_SIMLIB/SDSS_2007.SIMLIB"

    simlib05 = SIMLIBReader(simlib_file=simlib_file05)
    simlib06 = SIMLIBReader(simlib_file=simlib_file06)
    simlib07 = SIMLIBReader(simlib_file=simlib_file07)

    simlib_array = np.array([simlib05, simlib06, simlib07])

    # mag_offsets = lsstt.sims.choose_magoffset(n = n_sne_req)
    # z_obs = lsstt.sims.choose_z_volume_SFR(n_req=n_sne_req, zmax=z_max)

    ## This is faster than choose_z_volume_SFR -
    ## choose_z_volume_SFR is faster to generate large numbers
    binsize = 0.01
    z = np.arange(0.0, z_max, binsize)
    z_dz = np.arange(0.0 + binsize, z_max + binsize, binsize)

    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

    v_z = cosmo.comoving_volume(z)
    v_z_dz = cosmo.comoving_volume(z_dz)

    v_dz = v_z_dz - v_z

    norm_v_dz = v_dz / np.nanmax(v_dz)

    sfr_z = lsstt.sims.calculate_SFR(z)
    sfr_norm = sfr_z / np.nanmax(sfr_z)

    volumetric_rate = norm_v_dz * sfr_norm
    normed_volumetric_rate = volumetric_rate / np.nanmax(volumetric_rate)
    ## Generate probability density function
    pdf = InterpolatedUnivariateSpline(z, normed_volumetric_rate)


    info = pcc.classes.InfoClass()
    # info.load()
    # use only good SNe
    info.load(path = "/Users/berto/Code/CoCo/data/info/info_good.dat")

    filter_names = ["SDSS_u","SDSS_g","SDSS_r","SDSS_i","SDSS_z"]

    zp_dict = {}
    for i in filter_names:
        zp_dict[i] = pcc.kcorr.calc_AB_zp(i)

    n_sne = 4730
    n = 180570
    lc_outdir = "/Users/berto/projects/SDSS_sims/lcs/"
    log_outdir = "/Users/berto/projects/SDSS_sims/logs/"

    outfile = lc_outdir + "SN_"
    logfile = log_outdir + "LOG_SN_"

    while n_sne < n_sne_req:
        # try:
        n += 1
        # print(n)
        if log:
            logpath = logfile + str(n_sne+1).rjust(6, "0") +".json"

        ## Choose simlib
        sl = np.random.choice(simlib_array)

        if verbose: print(sl)

        ## Choose Field
        field_index = np.random.choice(len(sl.get_fields()))
        field = sl.get_fields()[field_index]

        if verbose: print(field_index, field)

        ## Choose CCD
        CCD_index = np.random.choice(sl.get_ccds(field))
        if verbose: print(CCD_index)

        ## Get obslog
        obslog = sl.get_obslog(field=field, ccd=CCD_index)
        if verbose: print(obslog.head())

        ## Convert obslog filters
        obslog["flt_b"] = obslog["flt"].map(lsstt.sims.string_format_for_mapping)
        if verbose: print(obslog)

        filters_to_sim = obslog.flt_b.values
        mjd_to_sim = obslog.mjd.values

        ## choose z
        # z_sim = lsstt.sims.choose_z(z_max = z_max)
        # if len(z_obs) > 1:
        #     w_z = np.random.randint(0, len(z_obs)-1)
        # else:
        #     w_z = 0
        #
        # z_sim = z_obs[w_z]
        # z_sim = lsstt.sims.choose_z_volume_SFR(n_req=1, zmax=z_max)
        n_z = 0
        while n_z < 1:
            x = np.random.random() * z_max
            y = np.random.random()

            if y <= pdf(x):
                z_sim = (x)
                n_z += 1

        if verbose: print("z= ", z_sim)

        ## choose MWEBV
        MW_EBV = obslog.mwebv.mean()
        if verbose: print(MW_EBV)

        ## Choose MagOffset
        # mag_offset = np.random.choice(mag_offsets)
        mag_offset = lsstt.sims.choose_magoffset(n=1)[0]
        if verbose: print("magoffset = ", mag_offset)

        ## Choose HostEBV
        host_EBV = lsstt.sims.choose_extinction_host(n=1)
        if verbose: print(MW_EBV)

        ## Choose MJDmax
        mjdmax = lsstt.sims.choose_MJDmax(obslog, n=1)
        if verbose: print(mjdmax)


        ## Choose SN Type
        subtype = lsstt.sims.choose_subtype()
        if verbose: print(subtype)

        ## Choose SN
        w = np.where(info.table["Type"] == subtype)[0]
        snindex = np.random.choice(w)
        snname = pcc.utils.b(info.table["snname"].data[snindex])
        if verbose: print(w, snname, snindex)

        w_trim = np.logical_and(mjd_to_sim < mjdmax + 100, mjd_to_sim > mjdmax - 50)
        mjd_to_sim = mjd_to_sim[w_trim]
        filters_to_sim = filters_to_sim[w_trim]

        mjd_to_sim = mjd_to_sim - mjdmax
        mjd_to_sim = mjd_to_sim / (1.0 + z_sim)

        # snname = b"SN2011dh"
        # mag_offset = -2.0 ## Make Ia-like
        ## Simulate "Perfect" LC
        # flux, flux_err = coco.simulate(snname,
        #                                z_sim, mag_offset, MW_EBV, host_EBV, 3.1,
        #                                mjdmax, mjd_to_sim,
        #                                filters_to_sim)
        flux, flux_err = coco.simulate(snname,
                                       z_sim, mag_offset, MW_EBV, host_EBV, 3.1,
                                       0.0, mjd_to_sim,
                                       filters_to_sim)

        mjd_to_sim = mjd_to_sim * (1.0 + z_sim)
        mjd_to_sim = mjd_to_sim + mjdmax

        ## TODO NEED TO TIME DILATE THE LCS - DONE, NEED TO TEST

        # flux, flux_err = coco.simulate(snname,
        #                                z_obs, 0.0, 0.0, 0.0, 3.1,
        #                                mjdmax, mjd_to_sim,
        #                                filters_to_sim)

        if verbose: print(flux)
        # if len(flux[~np.isnan(flux)]) > 5: ## pycoco classes don't like nans!
        if len(flux[~np.isnan(flux)]) > 2:  ## pycoco classes don't like nans!
            p = pcc.classes.PhotometryClass()
            p.load_table(pcc.utils.simulate_out_to_ap_table(mjd_to_sim, flux, flux_err, filters_to_sim), verbose=False)
            if plot: p.plot(enforce_zero=True)

            ## calculate zeropoints and convert to mag
            p_df = p.phot.to_pandas()
            p_df["zp"] = p_df["filter"].map(zp_dict)
            if verbose: print(p_df)
            p_df["mag"] = -2.5 * np.log10(p_df.flux) - p_df.zp

            ## Add noise - returns units of
            flux, flux_err = lcs.simulate(p_df["mag"], obslog, unit="ab")

            if plot:
                plt.errorbar(p_df.MJD, flux, yerr=flux_err, fmt="o")
                plt.show()

            w_detected = np.where((~np.isnan(flux.values)) & ((flux.values/flux_err.values) > 5))[0]
            # if verbose:
            if verbose:
                print("Flux, Flux_err")
                print(flux, flux_err)
            if verbose:
                print("where detected, len")
                print(w_detected, len(w_detected))

            # if len(w_detected) >= 6:
            if len(w_detected) >= 2:
                if verbose: print("good sne")
                ## Remove redshift simulated at top of code from the list
                # z_obs = np.delete(z_obs, [w_z])

                p_df["flux"] = flux
                p_df["flux_err"] = flux_err
                p_df["#MJD"] = p_df["MJD"]
                p_df.fillna(0, inplace = True)
                full_out_path = outfile + str(n_sne + 1).rjust(6, "0") + ".dat"
                p_df[["#MJD", "flux", "flux_err", "filter"]].to_csv(full_out_path, sep=" ", index = False, )

                if log:
                    logdict = {}
                    for i in logvars:
                        if type(locals()[i]) == np.ndarray:
                            logdict[i] = locals()[i].tolist()
                        elif type(locals()[i]) == np.int64:
                            logdict[i] = int(locals()[i])
                        elif type(locals()[i]) == pd.Series:
                            logdict[i] = locals()[i].to_json()
                        elif type(locals()[i]) == bytes:
                            logdict[i] = str(locals()[i], "utf-8")
                        else:
                            logdict[i] = locals()[i]

                    with open(logpath, "w") as ofile:
                        json.dumps(logdict, sort_keys=True,
                                   indent=4, separators=(',', ': '))
                        #     for i in logvars:
                        json.dump(logdict, ofile, sort_keys=True,
                                  indent=4, separators=(',', ': '))

                    #         ofile.write(str(i) + " " + str(locals()[i]) + "\n")
                        ofile.close()

                n_sne += 1

            # if logall:
            #     if log:
            #         with open(logpath, "w") as ofile:
            #             for i in logvars:
            #                 ofile.write(str(i) + " " + str(locals()[i]) + "\n")
        else:
            # if logall:
            #     if log:
            #         with open(logpath, "w") as ofile:
            #             for i in logvars:
            #                 ofile.write(str(i) + " " + str(locals()[i]) + "\n")

            pass
        # except:
        #     pass

    print("number passed = ", n_sne)
    print("number simulated = ", n)
    print("rejection rate %:")
    print(round(100.*(n - n_sne)/n, 2))
else:
    pass
