#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CSVPATH = "./out"
PLOTPATH = "./plot"
ERRPATH = f"{PLOTPATH}/error"
XSCALE = "log"
YSCALE = "log"
LINXSCALE = "log"
LINYSCALE = "log"
ERRXSCALE = "log"
ERRYSCALE = "log"

def main():
    try:
        infcell = pd.read_csv(f"{CSVPATH}/infcell.csv")

        # Intensity
        plt.plot(infcell['time'], infcell['intensity1'], color='b', label="Material 1")
        plt.plot(infcell['time'], infcell['intensity2'], color='r', label="Material 2")
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(infcell['time'], infcell['temp1'], color='b', label="Material 1")
        plt.plot(infcell['time'], infcell['temp2'], color='r', label="Material 2")
        plt.title("Energy Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Energy (erg/cm^3)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_temperature.png")
        plt.cla()
        plt.clf()

        # Energy Difference
        plt.plot(infcell['time'], infcell['energydiff'], color='g')
        plt.title("Energy Conservation")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Delta Energy (erg)")
        plt.xscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/energy_conservation.png")
        plt.cla()
        plt.clf()

        # Temperature Scale Figure
        plt.plot(infcell['time'], infcell['tintensity1'], color='b', label="Radiation 1")
        plt.plot(infcell['time'], infcell['ttemp1'], color='c', label="Material 1")
        plt.plot(infcell['time'], infcell['tintensity2'], color='r', label="Radiation 2")
        plt.plot(infcell['time'], infcell['ttemp2'], color='m', label="Material 2")
        plt.title("Temperature Domain - Radiation and Material")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Temperature (eV)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/temperature_space.png")
        plt.cla()
        plt.clf()
    except Exception as e:
        pass

    try:
        linear = pd.read_csv(f"{CSVPATH}/analytic_linear.csv")

        # Intensity
        plt.plot(linear['time'], linear['intensity1'], color='b', label="Material 1")
        plt.plot(linear['time'], linear['intensity2'], color='r', label="Material 2")
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/linear_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(linear['time'], linear['energy1'], color='b', label="Material 1")
        plt.plot(linear['time'], linear['energy2'], color='r', label="Material 2")
        plt.title("Energy Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Energy (erg/cm^3)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/linear_temperature.png")
        plt.cla()
        plt.clf()
    except Exception as e:
        pass

    try:
        infcelllinear = pd.read_csv(f"{CSVPATH}/infcelllinear.csv")

        # Intensity
        plt.plot(infcelllinear['time'], infcelllinear['intensity1'], color='b', label="Material 1")
        plt.plot(infcelllinear['time'], infcelllinear['intensity2'], color='r', label="Material 2")
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/iclinear_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(infcelllinear['time'], infcelllinear['energy1'], color='b', label="Material 1")
        plt.plot(infcelllinear['time'], infcelllinear['energy2'], color='r', label="Material 2")
        plt.title("Energy Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Energy (erg/cm^3)")
        plt.xscale(XSCALE)
        plt.yscale(YSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/iclinear_energy.png")
        plt.cla()
        plt.clf()
    except Exception as e:
        pass

    try:
        infcell = pd.read_csv(f"{CSVPATH}/infcell.csv")
        linear = pd.read_csv(f"{CSVPATH}/analytic_linear.csv")

        # Intensity
        plt.plot(infcell['time'], infcell['intensity1'], color='b', label="Material 1 - MC")
        plt.plot(infcell['time'], infcell['intensity2'], color='r', label="Material 2 - MC")
        plt.plot(linear['time'], linear['intensity1'], color='c', label="Material 1 - Linear", linestyle='--')
        plt.plot(linear['time'], linear['intensity2'], color='m', label="Material 2 - Linear", linestyle='--')
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale(LINXSCALE)
        plt.yscale(LINYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_intensity_linear.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(infcell['time'], infcell['temp1'], color='b', label="Material 1 - MC")
        plt.plot(infcell['time'], infcell['temp2'], color='r', label="Material 2 - MC")
        plt.plot(linear['time'], linear['energy1'], color='c', label="Material 1 - Linear", linestyle='--')
        plt.plot(linear['time'], linear['energy2'], color='m', label="Material 2 - Linear", linestyle='--')
        plt.title("Energy Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Energy (erg/cm^3)")
        plt.xscale(LINXSCALE)
        plt.yscale(LINYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_temperature_linear.png")
        plt.cla()
        plt.clf()

        arr_len = len(linear['time'])

        re_intensity_1 = np.zeros(arr_len)
        re_intensity_2 = np.zeros(arr_len)
        re_energy_1 = np.zeros(arr_len)
        re_energy_2 = np.zeros(arr_len)
        for index, intensity in enumerate(linear['intensity1']):
            re_intensity_1[index] = np.abs(infcell['intensity1'][index] - intensity) / intensity
        for index, intensity in enumerate(linear['intensity2']):
            re_intensity_2[index] = np.abs(infcell['intensity2'][index] - intensity) / intensity
        for index, energy in enumerate(linear['energy1']):
            re_energy_1[index] = np.abs(infcell['temp1'][index] - energy) / energy
        for index, energy in enumerate(linear['energy2']):
            re_energy_2[index] = np.abs(infcell['temp2'][index] - energy) / energy

        print("INFCELL:")
        print(f"Max Relative Error Intensity 1: {np.max(re_intensity_1)}")
        print(f"Mean Relative Error Intensity 1: {np.mean(re_intensity_1)}")
        print(f"Max Relative Error Intensity 2: {np.max(re_intensity_2)}")
        print(f"Mean Relative Error Intensity 2: {np.mean(re_intensity_2)}")
        print(f"Max Relative Error Energy 1: {np.max(re_energy_1)}")
        print(f"Mean Relative Error Energy 1: {np.mean(re_energy_1)}")
        print(f"Max Relative Error Energy 2: {np.max(re_energy_2)}")
        print(f"Mean Relative Error Energy 2: {np.mean(re_energy_2)}\n")

        plt.plot(linear['time'], re_intensity_1, label="Relative Error Intensity 1")
        plt.plot(linear['time'], re_intensity_2, label="Relative Error Intensity 2")
        plt.plot(linear['time'], re_energy_1, label="Relative Error Energy 1")
        plt.plot(linear['time'], re_energy_2, label="Relative Error Energy 2")
        plt.title("Relative Error Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Relative Error")
        plt.xscale(ERRXSCALE)
        plt.yscale(ERRYSCALE)
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{ERRPATH}/infcell_relative_error.png")
        plt.cla()
        plt.clf()
    except Exception as e:
        raise e

if __name__ == '__main__':
    print("Plotting...")
    main()
    print("Done")
