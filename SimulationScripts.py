#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import nBodySimulator as nbs
import time
import numpy as np


nBodies=800;
accelerationThreshold=0.0002;
Frequency=50
deltaT=60;
nTimesteps=1000;
timeTaken=np.zeros((10,1))
for i in range(8,11):
    t = time.time()
    for j in range(10):
        #nbs.executeSimulationAssuming(i*50,deltaT,nTimesteps)
        nbs.executeApproximateSimulationAssuming(i*50,deltaT,nTimesteps,accelerationThreshold,Frequency)
    timeTaken[i-1] = time.time() - t
    timeTaken/=10
    print(timeTaken[i-1])