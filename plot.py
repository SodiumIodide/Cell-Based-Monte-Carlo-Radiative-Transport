#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CSVPATH = "./out"
PLOTPATH = "./plot"

def main():
    infcell = pd.read_csv(f"{CSVPATH}/infcell.csv")

    # Intensity
    plt.plot(infcell['time'], infcell['intensity1'], color='b', label="Material 1")
    plt.plot(infcell['time'], infcell['intensity2'], color='r', label="Material 2")
    plt.title("Intensity Plot")
    plt.xlabel("Time - ct (cm)")
    plt.ylabel("Intensity (erg/cm^2-s)")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(b=True, which="both", axis="both")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(f"{PLOTPATH}/infcell_intensity.png")
    plt.cla()
    plt.clf()

    # Temperature
    plt.plot(infcell['time'], infcell['temp1'], color='b', label="Material 1")
    plt.plot(infcell['time'], infcell['temp2'], color='r', label="Material 2")
    plt.title("Temperature Plot")
    plt.xlabel("Time - ct (cm)")
    plt.ylabel("Temperature (eV)")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(b=True, which="both", axis="both")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(f"{PLOTPATH}/infcell_temperature.png")
    plt.cla()
    plt.clf()

if __name__ == '__main__':
    main()
