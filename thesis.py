#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import CUSUM_suplib as csl
from basicCUSUM import *

def CUSUM_NoUpdate_Plot():
    csv = pd.read_csv("seq00.csv")
    cusumlistupper = []
    cusumlistlower = []
    mean = csl.FirstNMean(csv["seq0"], 20)
    segments = 3
    k = csl.TopNPercent(csv["seq0"], 5) - csl.BotNPercent(csv["seq0"], 5)
    k /= segments
    delta = k*2
    h = k*10
    cusum = CUSUM(mean, delta, h)

    for each in csv["seq0"]:
        cusum.Update(each)
        cusumlistupper.append(cusum.GetUpperCusum())
        cusumlistlower.append(cusum.GetLowerCusum())

    plt.figure(1)
    plt.figure(1).suptitle('CUSUM without update')
    plt.subplot(211)
    plt.title('upper')
    plt.plot(cusumlistupper)

    plt.subplot(212)
    plt.title('data')
    plt.plot(csv["seq0"])

    plt.figure(2)
    plt.figure(2).suptitle('CUSUM without update')
    plt.subplot(211)
    plt.title('lower')
    plt.plot(cusumlistlower)

    plt.subplot(212)
    plt.title('data')
    plt.plot(csv["seq0"])

def CUSUM_WUpdate_Plot():
    csv = pd.read_csv("seq00.csv")
    cusumlistupper = []
    cusumlistlower = []
    mean = csl.FirstNMean(csv["seq0"], 20)
    segments = 3
    k = csl.TopNPercent(csv["seq0"], 5) - csl.BotNPercent(csv["seq0"], 5)
    k /= segments
    delta = k*2
    h = k*10
    cusum = CUSUM2(mean, delta, h)

    for each in csv["seq0"]:
        cusum.Update2(each)
        cusumlistupper.append(cusum.GetUpperCusum())
        cusumlistlower.append(cusum.GetLowerCusum())

    plt.figure(3)
    plt.figure(3).suptitle('CUSUM with update')
    plt.subplot(211)
    plt.title('upper')
    plt.plot(cusumlistupper)

    plt.subplot(212)
    plt.title('data')
    plt.plot(csv["seq0"])

    plt.figure(4)
    plt.figure(4).suptitle('CUSUM with update')
    plt.subplot(211)
    plt.title('lower')
    plt.plot(cusumlistlower)

    plt.subplot(212)
    plt.title('data')
    plt.plot(csv["seq0"])

def Distances_vs_Angles_Plot():
    distances = pd.read_csv("seq00.csv")
    angles = pd.read_csv("angle.csv")

    plt.figure(5)
    plt.figure(5).suptitle('Sequence 0')
    plt.subplot(211)
    plt.title('Distance')
    plt.plot(distances)

    plt.subplot(212)
    plt.title('Angle')
    plt.plot(angles)

def Angles_Without_Spikes_Plot():
    distances = pd.read_csv("seq00.csv")
    angles = pd.read_csv("angle.csv")
    angles_wo_spikes = csl.RemoveSpikes(list(angles["angle"]), -0.744)

    plt.figure(6)
    plt.figure(6).suptitle('Sequence 0')
    plt.subplot(211)
    plt.title('Angles')
    plt.plot(angles)

    plt.subplot(212)
    plt.title('Angles without spikes')
    plt.plot(angles_wo_spikes)

def CUSUM_Angles_Plot():
    csv = pd.read_csv("angle.csv")
    cusumlistupper = []
    cusumlistlower = []
    mean = csl.FirstNMean(csv["angle"], 20)
    segments = 3
    k = csl.TopNPercent(csv["angle"], 5) - csl.BotNPercent(csv["angle"], 5)
    k /= segments
    delta = k*2
    h = k*10
    cusum = CUSUM2(mean, delta, h)
    for each in csv["angle"]:
        cusum.Update2(each)
        cusumlistupper.append(cusum.GetUpperCusum())
        cusumlistlower.append(cusum.GetLowerCusum())

    plt.figure(7)
    plt.figure(7).suptitle('CUSUM of angles')
    plt.subplot(211)
    plt.title('upper')
    plt.plot(cusumlistupper)

    plt.subplot(212)
    plt.title('data')
    plt.plot(csv["angle"])

    plt.figure(8)
    plt.figure(8).suptitle('CUSUM of angles')
    plt.subplot(211)
    plt.title('lower')
    plt.plot(cusumlistlower)

    plt.subplot(212)
    plt.title('data')
    plt.plot(csv["angle"])

def CUSUM_Angles_wo_Spikes_Plot():
    csv = pd.read_csv("angle.csv")
    angles = csl.RemoveSpikes(list(csv["angle"]), -0.744)
    cusumlistupper = []
    cusumlistlower = []
    mean = csl.FirstNMean(angles, 20)
    segments = 3
    k = csl.TopNPercent(angles, 5) - csl.BotNPercent(angles, 5)
    k /= segments
    delta = k*2
    h = k*10
    cusum = CUSUM2(mean, delta, h)
    for each in angles:
        cusum.Update2(each)
        cusumlistupper.append(cusum.GetUpperCusum())
        cusumlistlower.append(cusum.GetLowerCusum())

    plt.figure(9)
    plt.figure(9).suptitle('CUSUM of angles')
    plt.subplot(211)
    plt.title('upper')
    plt.plot(cusumlistupper)

    plt.subplot(212)
    plt.title('data')
    plt.plot(angles)

    plt.figure(10)
    plt.figure(10).suptitle('CUSUM of angles')
    plt.subplot(211)
    plt.title('lower')
    plt.plot(cusumlistlower)

    plt.subplot(212)
    plt.title('data')
    plt.plot(angles)



def main():
    # Distances_vs_Angles_Plot()
    # Angles_Without_Spikes_Plot()
    # CUSUM_NoUpdate_Plot()
    # CUSUM_WUpdate_Plot()
    # CUSUM_Angles_Plot()
    # CUSUM_Angles_wo_Spikes_Plot()

    plt.show()

if __name__ == "__main__":
    main()
