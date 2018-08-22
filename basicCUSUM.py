#!/usr/bin/env python3

# Gwenaël Gatto
# CUSUM module

# This module contains a basic CUSUM implementation
# A CUSUM object can be created that takes in:
# a mean, change to be detected, threshold to raise alarm, and
# optionally a CUSUM value (default to 0 otherwise)

# d0 = e0 = 0
# d1, d2, ... and e1, e2, ... are calculated recursively
# dl = max[0, dl-1 + (xbar - (mu0 + k))]    upper limit
# el = max[0, el-1 - (xbar - (mu0 - k))]    lower limit
# k = delta/2, where delta is the change to be detected
# alarm is raised if at time r, dr > h or er > h
# dr > h means process has shifted to greater value
# er > h means process has shifted to lower value

import matplotlib.pyplot as plt
import pandas as pd

class CUSUM(object):

    def __init__(self, mean, delta, h, dlcusum=0, elcusum=0):
        # from arguments
        self.mean = mean
        self.delta = delta
        self.h = h
        self.dlcusum = dlcusum
        self.elcusum = elcusum

        # derived
        self.k = delta/2
        self.dlvariance = self.mean + self.k
        self.elvariance = self.mean - self.k

        self.dlalarm = False
        self.elalarm = False

    # Core methods
    def Update(self, xbar):
        # CUSUM values
        self.dlcusum = max(0, self.dlcusum + (xbar - self.dlvariance))
        self.elcusum = max(0, self.elcusum - (xbar - self.elvariance))

        # Alarm state
        if self.h < self.dlcusum:
            self.dlalarm = True
        else:
            self.dlalarm = False
        if self.h < self.elcusum:
            self.elalarm = True
        else:
            self.elalarm = False


    def WhichAlarm(self):
        if self.dlalarm:
            return 1
        else:
            return 0

    def Alarm(self):
        return self.dlalarm or self.elalarm


    # Utility methods
    def Print(self, todo=None):
        if todo!=None:
            print("CUSUM+: ", self.dlcusum)
            print("CUSUM-: ", self.elcusum)
            print("Mean: ", self.mean)
            print("Delta: ", self.delta)
            print("Threshold: ", self.h)
        else:
            print("CUSUM +/-: ", self.dlcusum, "/", self.elcusum)

    def __str__(self):
        return ("+CUSUM: " + str(self.dlcusum) + \
                "\n-CUSUM: " + str(self.elcusum))

    # Accessors/Mutators
    def GetUpperCusum(self):
        return self.dlcusum

    def GetLowerCusum(self):
        return self.elcusum

###############################################################################
# CUSUM class that resets its mean
# New mean is m0+k or m0-k depending which alarm has been raised
###############################################################################
class CUSUM2(CUSUM):
    def __init__(self, mean, delta, h, dlcusum=0, elcusum=0):
        CUSUM.__init__(self, mean, delta, h, dlcusum, elcusum)


    def Update2(self, xbar):
        CUSUM.Update(self, xbar)
        if CUSUM.Alarm(self):
            print("Alarm Raised")
            self.Reinitialize(CUSUM.WhichAlarm(self))
            print("New mean: ", self.mean)


    def Reinitialize(self, alarm):
        if alarm:
            self.mean = self.dlvariance
        else:
            self.mean = self.elvariance

        self.dlcusum = 0
        self.elcusum = 0

        # derived
        self.dlvariance = self.mean + self.k
        self.elvariance = self.mean - self.k

        self.dlalarm = False
        self.elalarm = False
