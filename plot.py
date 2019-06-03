#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CSVPATH = "./out"
PLOTPATH = "./plot"
ERRPATH = f"{PLOTPATH}/error_v2"

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
    except Exception as e:
        pass

    try:
        infcell_2 = pd.read_csv(f"{CSVPATH}/infcell_2.csv")

        # Intensity
        plt.plot(infcell_2['time'], infcell_2['intensity1'], color='b', label="Material 1")
        plt.plot(infcell_2['time'], infcell_2['intensity2'], color='r', label="Material 2")
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_2_intensity.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(infcell_2['time'], infcell_2['temp1'], color='b', label="Material 1")
        plt.plot(infcell_2['time'], infcell_2['temp2'], color='r', label="Material 2")
        plt.title("Energy Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Energy (erg/cm^3)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_2_temperature.png")
        plt.cla()
        plt.clf()
    except Exception as e:
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
        plt.title("Energy Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Energy (erg/cm^3)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_3_temperature.png")
        plt.cla()
        plt.clf()
    except Exception as e:
        pass

    try:
        infcell_1 = pd.read_csv(f"{CSVPATH}/infcell_1.csv")
        linear = pd.read_csv(f"{CSVPATH}/analytic_linear.csv")

        # Intensity
        plt.plot(infcell_1['time'], infcell_1['intensity1'], color='b', label="Material 1 - MC")
        plt.plot(infcell_1['time'], infcell_1['intensity2'], color='r', label="Material 2 - MC")
        plt.plot(linear['time'], linear['intensity1'], color='c', label="Material 1 - Linear", linestyle='--')
        plt.plot(linear['time'], linear['intensity2'], color='m', label="Material 2 - Linear", linestyle='--')
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_1_intensity_linear.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(infcell_1['time'], infcell_1['temp1'], color='b', label="Material 1 - MC")
        plt.plot(infcell_1['time'], infcell_1['temp2'], color='r', label="Material 2 - MC")
        plt.plot(linear['time'], linear['energy1'], color='c', label="Material 1 - Linear", linestyle='--')
        plt.plot(linear['time'], linear['energy2'], color='m', label="Material 2 - Linear", linestyle='--')
        plt.title("Energy Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Energy (erg/cm^3)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_1_temperature_linear.png")
        plt.cla()
        plt.clf()

        """ arr_len = len(linear['time'])

        re_intensity_1 = np.zeros(arr_len)
        re_intensity_2 = np.zeros(arr_len)
        re_energy_1 = np.zeros(arr_len)
        re_energy_2 = np.zeros(arr_len)
        for index, intensity in enumerate(linear['intensity1']):
            re_intensity_1[index] = np.abs(infcell_1['intensity1'][index] - intensity) / intensity
        for index, intensity in enumerate(linear['intensity2']):
            re_intensity_2[index] = np.abs(infcell_1['intensity2'][index] - intensity) / intensity
        for index, energy in enumerate(linear['energy1']):
            re_energy_1[index] = np.abs(infcell_1['temp1'][index] - energy) / energy
        for index, energy in enumerate(linear['energy2']):
            re_energy_2[index] = np.abs(infcell_1['temp2'][index] - energy) / energy

        print("INFCELL_1:")
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
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{ERRPATH}/infcell_1_relative_error.png")
        plt.cla()
        plt.clf() """
    except Exception as e:
        raise e

    try:
        infcell_2 = pd.read_csv(f"{CSVPATH}/infcell_2.csv")
        linear = pd.read_csv(f"{CSVPATH}/analytic_linear.csv")

        # Intensity
        plt.plot(infcell_2['time'], infcell_2['intensity1'], color='b', label="Material 1 - MC")
        plt.plot(infcell_2['time'], infcell_2['intensity2'], color='r', label="Material 2 - MC")
        plt.plot(linear['time'], linear['intensity1'], color='c', label="Material 1 - Linear", linestyle='--')
        plt.plot(linear['time'], linear['intensity2'], color='m', label="Material 2 - Linear", linestyle='--')
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_2_intensity_linear.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(infcell_2['time'], infcell_2['temp1'], color='b', label="Material 1 - MC")
        plt.plot(infcell_2['time'], infcell_2['temp2'], color='r', label="Material 2 - MC")
        plt.plot(linear['time'], linear['energy1'], color='c', label="Material 1 - Linear", linestyle='--')
        plt.plot(linear['time'], linear['energy2'], color='m', label="Material 2 - Linear", linestyle='--')
        plt.title("Energy Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Energy (erg/cm^3)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_2_temperature_linear.png")
        plt.cla()
        plt.clf()

        """ arr_len = len(linear['time'])

        re_intensity_1 = np.zeros(arr_len)
        re_intensity_2 = np.zeros(arr_len)
        re_energy_1 = np.zeros(arr_len)
        re_energy_2 = np.zeros(arr_len)
        for index, intensity in enumerate(linear['intensity1']):
            re_intensity_1[index] = np.abs(infcell_2['intensity1'][index] - intensity) / intensity
        for index, intensity in enumerate(linear['intensity2']):
            re_intensity_2[index] = np.abs(infcell_2['intensity2'][index] - intensity) / intensity
        for index, energy in enumerate(linear['energy1']):
            re_energy_1[index] = np.abs(infcell_2['temp1'][index] - energy) / energy
        for index, energy in enumerate(linear['energy2']):
            re_energy_2[index] = np.abs(infcell_2['temp2'][index] - energy) / energy

        print("INFCELL_2:")
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
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{ERRPATH}/infcell_2_relative_error.png")
        plt.cla()
        plt.clf() """
    except Exception as e:
        raise e

    try:
        infcell_3 = pd.read_csv(f"{CSVPATH}/infcell_3.csv")
        linear = pd.read_csv(f"{CSVPATH}/analytic_linear.csv")

        # Intensity
        plt.plot(infcell_3['time'], infcell_3['intensity1'], color='b', label="Material 1 - MC")
        plt.plot(infcell_3['time'], infcell_3['intensity2'], color='r', label="Material 2 - MC")
        plt.plot(linear['time'], linear['intensity1'], color='c', label="Material 1 - Linear", linestyle='--')
        plt.plot(linear['time'], linear['intensity2'], color='m', label="Material 2 - Linear", linestyle='--')
        plt.title("Intensity Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Intensity (erg/cm^2-s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_3_intensity_linear.png")
        plt.cla()
        plt.clf()

        # Temperature
        plt.plot(infcell_3['time'], infcell_3['temp1'], color='b', label="Material 1 - MC")
        plt.plot(infcell_3['time'], infcell_3['temp2'], color='r', label="Material 2 - MC")
        plt.plot(linear['time'], linear['energy1'], color='c', label="Material 1 - Linear", linestyle='--')
        plt.plot(linear['time'], linear['energy2'], color='m', label="Material 2 - Linear", linestyle='--')
        plt.title("Energy Plot")
        plt.xlabel("Time - ct (cm)")
        plt.ylabel("Energy (erg/cm^3)")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{PLOTPATH}/infcell_3_temperature_linear.png")
        plt.cla()
        plt.clf()

        """ arr_len = len(linear['time'])

        re_intensity_1 = np.zeros(arr_len)
        re_intensity_2 = np.zeros(arr_len)
        re_energy_1 = np.zeros(arr_len)
        re_energy_2 = np.zeros(arr_len)
        for index, intensity in enumerate(linear['intensity1']):
            re_intensity_1[index] = np.abs(infcell_3['intensity1'][index] - intensity) / intensity
        for index, intensity in enumerate(linear['intensity2']):
            re_intensity_2[index] = np.abs(infcell_3['intensity2'][index] - intensity) / intensity
        for index, energy in enumerate(linear['energy1']):
            re_energy_1[index] = np.abs(infcell_3['temp1'][index] - energy) / energy
        for index, energy in enumerate(linear['energy2']):
            re_energy_2[index] = np.abs(infcell_3['temp2'][index] - energy) / energy

        print("INFCELL_3:")
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
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(b=True, which="both", axis="both")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(f"{ERRPATH}/infcell_3_relative_error.png")
        plt.cla()
        plt.clf() """
    except Exception as e:
        raise e

if __name__ == '__main__':
    print("Plotting...")
    main()
    print("Done")
