#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CSVPATH = "./out"
PLOTPATH = "./plot"
XSCALE = "log"
YSCALE = "log"

def main():
    '''Main method'''
    nonlinear_exists = True
    try:
        nonlinear = pd.read_csv(f"{CSVPATH}/nonlinear.csv")
    except FileNotFoundError:
        nonlinear_exists = False
    nonlinearmc_exists = True
    try:
        nonlinear_mc = pd.read_csv(f"{CSVPATH}/nonlinearmc.csv")
    except FileNotFoundError:
        nonlinearmc_exists = False
    if nonlinear_exists:
        # Standard deviation computing
        nl_std_intensity_1 = nonlinear['varintensity1'].apply(np.abs).apply(np.sqrt)
        nl_std_intensity_2 = nonlinear['varintensity2'].apply(np.abs).apply(np.sqrt)
        nl_std_temp_1 = nonlinear['vartemperature1'].apply(np.abs).apply(np.sqrt)
        nl_std_temp_2 = nonlinear['vartemperature2'].apply(np.abs).apply(np.sqrt)
        nl_lb_intensity_1 = nonlinear['intensity1'] - nl_std_intensity_1
        nl_ub_intensity_1 = nonlinear['intensity1'] + nl_std_intensity_1
        nl_lb_intensity_2 = nonlinear['intensity2'] - nl_std_intensity_2
        nl_ub_intensity_2 = nonlinear['intensity2'] + nl_std_intensity_2
        nl_lb_temp_1 = nonlinear['temperature1'] - nl_std_temp_1
        nl_ub_temp_1 = nonlinear['temperature1'] + nl_std_temp_1
        nl_lb_temp_2 = nonlinear['temperature2'] - nl_std_temp_2
        nl_ub_temp_2 = nonlinear['temperature2'] + nl_std_temp_2

        # Intensity
        plt.plot(nonlinear['time'], nonlinear['intensity1'], color='b', label="Material 1")
        plt.plot(nonlinear['time'], nonlinear['intensity2'], color='r', label="Material 2")
        plt.plot(nonlinear['time'], nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_intensity_1, nl_ub_intensity_1, color='b', alpha=0.5)
        plt.plot(nonlinear['time'], nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_intensity_2, nl_ub_intensity_2, color='r', alpha=0.5)
        #plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm$^2$-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/nonlinear_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear['time'], nonlinear['temperature1'], color='b', label="Material 1")
        plt.plot(nonlinear['time'], nonlinear['temperature2'], color='r', label="Material 2")
        plt.plot(nonlinear['time'], nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_temp_1, nl_ub_temp_1, color='b', alpha=0.5)
        plt.plot(nonlinear['time'], nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear['time'], nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear['time'], nl_lb_temp_2, nl_ub_temp_2, color='r', alpha=0.5)
        #plt.title("Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/nonlinear_temperature.png")
        plt.cla()
        plt.clf()

    if nonlinearmc_exists:
        # Standard deviation computing
        mc_nl_std_intensity_1 = nonlinear_mc['varintensity1'].apply(np.abs).apply(np.sqrt)
        mc_nl_std_intensity_2 = nonlinear_mc['varintensity2'].apply(np.abs).apply(np.sqrt)
        mc_nl_std_temp_1 = nonlinear_mc['vartemperature1'].apply(np.abs).apply(np.sqrt)
        mc_nl_std_temp_2 = nonlinear_mc['vartemperature2'].apply(np.abs).apply(np.sqrt)
        mc_nl_lb_intensity_1 = nonlinear_mc['intensity1'] - mc_nl_std_intensity_1
        mc_nl_ub_intensity_1 = nonlinear_mc['intensity1'] + mc_nl_std_intensity_1
        mc_nl_lb_intensity_2 = nonlinear_mc['intensity2'] - mc_nl_std_intensity_2
        mc_nl_ub_intensity_2 = nonlinear_mc['intensity2'] + mc_nl_std_intensity_2
        mc_nl_lb_temp_1 = nonlinear_mc['temperature1'] - mc_nl_std_temp_1
        mc_nl_ub_temp_1 = nonlinear_mc['temperature1'] + mc_nl_std_temp_1
        mc_nl_lb_temp_2 = nonlinear_mc['temperature2'] - mc_nl_std_temp_2
        mc_nl_ub_temp_2 = nonlinear_mc['temperature2'] + mc_nl_std_temp_2

        # Intensity
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity1'], color='b', label="Material 1")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['intensity2'], color='r', label="Material 2")
        plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear_mc['time'], mc_nl_lb_intensity_1, mc_nl_ub_intensity_1, color='b', alpha=0.5)
        plt.plot(nonlinear_mc['time'], mc_nl_lb_intensity_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_intensity_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear_mc['time'], mc_nl_lb_intensity_2, mc_nl_ub_intensity_2, color='r', alpha=0.5)
        #plt.title("Monte Carlo Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm$^2$-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/nonlinear_mc_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature1'], color='b', label="Material 1")
        plt.plot(nonlinear_mc['time'], nonlinear_mc['temperature2'], color='r', label="Material 2")
        plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_1, color='b', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_1, color='b', linestyle=':', label=None)
        plt.fill_between(nonlinear_mc['time'], mc_nl_lb_temp_1, mc_nl_ub_temp_1, color='b', alpha=0.5)
        plt.plot(nonlinear_mc['time'], mc_nl_lb_temp_2, color='r', linestyle=':', label=None)
        plt.plot(nonlinear_mc['time'], mc_nl_ub_temp_2, color='r', linestyle=':', label=None)
        plt.fill_between(nonlinear_mc['time'], mc_nl_lb_temp_2, mc_nl_ub_temp_2, color='r', alpha=0.5)
        #plt.title("Monte Carlo Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/nonlinear_mc_temperature.png")
        plt.cla()
        plt.clf()

if __name__ == '__main__':
    print("Plotting...")
    main()
    print("Done")
