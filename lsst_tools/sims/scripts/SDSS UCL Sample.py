# import matplotlib.pyplot as plt

import numpy as np

import pycoco as pcc
import pyCoCo as pccsim

# from astropy.table import Table

import lsst_tools as lsstt
from lcsim.simlib import SIMLIBReader
from lcsim.lcsim import LCSim

if __name__ == "__main__":
    print("Running __main__")
    # verbose = True
    verbose = False
    n_sne_req = 10000
    z_max = 0.3
    plot = False

    ## Initialise pyCoCo
    fltPath = b"/Users/berto/Code/CoCo/data/filters"
    # rootPath = b"/Users/berto/Code/CoCo"
    rootPath = b"/Users/berto/projects/stunt_CoCo"

    coco = pccsim.pyCoCo(fltPath, rootPath)
    lcs = LCSim()

    simlib_file05 = "/Users/berto/projects/SDSS_sims/SDSS_SIMLIB/SDSS_2005.SIMLIB"
    simlib_file06= "/Users/berto/projects/SDSS_sims/SDSS_SIMLIB/SDSS_2006.SIMLIB"
    simlib_file07 = "/Users/berto/projects/SDSS_sims/SDSS_SIMLIB/SDSS_2007.SIMLIB"

    simlib05 = SIMLIBReader(simlib_file=simlib_file05)
    simlib06 = SIMLIBReader(simlib_file=simlib_file06)
    simlib07 = SIMLIBReader(simlib_file=simlib_file07)

    simlib_array = np.array([simlib05, simlib06, simlib07])

    mag_offsets = lsstt.sims.choose_magoffset(n_sne_req)

    info = pcc.InfoClass()
    info.load()

    filter_names = ["SDSS_u","SDSS_g","SDSS_r","SDSS_i","SDSS_z"]

    zp_dict = {}
    for i in filter_names:
        zp_dict[i] = pcc.kcorr.calc_AB_zp(i)

    n_sne = 0
    n = 0
    outdir = "/Users/berto/projects/SDSS_sims/lcs/"
    outfile = outdir + "SN_"

    while n_sne < n_sne_req:
        try:
            n += 1
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
            z_sim = lsstt.sims.choose_z(z_max = z_max)
            if verbose: print("z= ", z_sim)

            ## choose MWEBV
            MW_EBV = obslog.mwebv.mean()
            if verbose: print(MW_EBV)

            ## Choose MagOffset
            mag_offset = np.random.choice(mag_offsets)
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
            snname = pcc.utils.b(info.table["snname"].data[w][0])
            if verbose: print(w, snname)

            snname = b"SN2011dh"

            ## Simulate "Perfect" LC
            flux, flux_err = coco.simulate(snname,
                                           z_sim, mag_offset, MW_EBV, host_EBV, 3.1,
                                           mjdmax, mjd_to_sim,
                                           filters_to_sim)
            #
            # flux, flux_err = coco.simulate(snname,
            #                                z_obs, 0.0, 0.0, 0.0, 3.1,
            #                                mjdmax, mjd_to_sim,
            #                                filters_to_sim)

            if verbose: print(flux)
            p = pcc.PhotometryClass()
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
                plt.errorbar(p_df.MJD, flux, yerr=flux_error, fmt="o")
                plt.show()

            w_detected = np.where((~np.isnan(flux.values)) & ((flux.values/flux_err.values) > 5))[0]
            # if verbose:
            if verbose:
                print("Flux, Flux_err")
                print(flux, flux_err)
            if verbose:
                print("where detected, len")
                print(w_detected, len(w_detected))

            if len(w_detected) >= 6:
                p_df["flux"] = flux
                p_df["flux_err"] = flux_err
                p_df["#MJD"] = p_df["MJD"]
                p_df.fillna(0, inplace = True)
                p_df[["#MJD", "flux", "flux_err", "filter"]].to_csv(outfile + str(n_sne+1) +".dat", sep=" ", index = False, )

                n_sne += 1
        except:
            pass
    print("rejection rate %"")
    print(n_sne/n)
else:
    pass
