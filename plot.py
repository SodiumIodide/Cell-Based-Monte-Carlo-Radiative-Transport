#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CSVPATH = "./out"
PLOTPATH = "./plot"

def main():
    try:
        infcell_1 = pd.read_csv(f"{CSVPATH}/infcell_1.csv")

        # Intensity
        plt.plot(infcell_1['time'], infcell_1['intensity1'], color='b', label="Material 1")
        plt.plot(infcell_1['time'], infcell_1['intensity2'], color='r', label="Material 2")
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_1_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(infcell_1['time'], infcell_1['temp1'], color='b', label="Material 1")
        plt.plot(infcell_1['time'], infcell_1['temp2'], color='r', label="Material 2")
        plt.title("Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_1_temperature.png")
        plt.cla()
        plt.clf()
    except:
        pass

    try:
        infcell_3 = pd.read_csv(f"{CSVPATH}/infcell_3.csv")

        # Intensity
        plt.plot(infcell_3['time'], infcell_3['intensity1'], color='b', label="Material 1")
        plt.plot(infcell_3['time'], infcell_3['intensity2'], color='r', label="Material 2")
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_3_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(infcell_3['time'], infcell_3['temp1'], color='b', label="Material 1")
        plt.plot(infcell_3['time'], infcell_3['temp2'], color='r', label="Material 2")
        plt.title("Temperature Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_3_temperature.png")
        plt.cla()
        plt.clf()
    except:
        pass

if __name__ == '__main__':
    main()
