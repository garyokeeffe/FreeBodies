#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import iqr
from scipy import sparse
gravitationalConstant=6.67408*(10**(-11))

def createRandomInitialPositionsFor(nBodies):
    return 0.01*np.random.randn(nBodies,3)

def createRandomMassesFor(nBodies):
    masses=1000*np.tile(np.random.rand(nBodies,1), [1,nBodies])
    return masses.transpose()

def createRandomInitialVelocitiesFor(nBodies):
    return np.random.randn(nBodies,3)

def calculateEdgeLengthsUsing(currentPositions):
    xTiled,yTiled,zTiled=tilePositionCoordinates(currentPositions)
    xEdgeLength=xTiled-xTiled.transpose()
    yEdgeLength=yTiled-yTiled.transpose()
    zEdgeLength=zTiled-zTiled.transpose()
    eucledianLengthOfEdgeSquared=np.array(np.square(xEdgeLength)+np.square(yEdgeLength)+np.square(zEdgeLength))
    return xEdgeLength,yEdgeLength,zEdgeLength,eucledianLengthOfEdgeSquared

def tilePositionCoordinates(currentPositions):
    nBodies=len(currentPositions)
    xTiled=np.tile(currentPositions[:,0], [nBodies,1])
    yTiled=np.tile(currentPositions[:,1], [nBodies,1])
    zTiled=np.tile(currentPositions[:,2], [nBodies,1])
    return xTiled,yTiled,zTiled

def calculateForcesOnBodiesUsingThe(currentPositions):
    xEdgeLength,yEdgeLength,zEdgeLength,eucledianLengthOfEdgeSquared=calculateEdgeLengthsUsing(currentPositions)
    xForce=gravitationalConstant*np.divide(xEdgeLength,(eucledianLengthOfEdgeSquared*np.sqrt(eucledianLengthOfEdgeSquared)))
    yForce=gravitationalConstant*np.divide(yEdgeLength,(eucledianLengthOfEdgeSquared*np.sqrt(eucledianLengthOfEdgeSquared)))
    zForce=gravitationalConstant*np.divide(zEdgeLength,(eucledianLengthOfEdgeSquared*np.sqrt(eucledianLengthOfEdgeSquared)))
    return removeSingularitiesFrom(xForce,yForce,zForce)
    
def removeSingularitiesFrom(xForce,yForce,zForce):
    xPositionsOfNaNs = np.isnan(xForce)
    yPositionsOfNaNs = np.isnan(yForce)
    zPositionsOfNaNs = np.isnan(zForce)
    xForce[xPositionsOfNaNs] = 0
    yForce[yPositionsOfNaNs] = 0
    zForce[zPositionsOfNaNs] = 0
    return xForce,yForce,zForce

def calculateAccelerationsOnEacyBodyUsingThe(currentPositions,masses):
    xForce,yForce,zForce=calculateForcesOnBodiesUsingThe(currentPositions)
    xAcceleration=np.multiply(xForce,masses)
    yAcceleration=np.multiply(yForce,masses)
    zAcceleration=np.multiply(zForce,masses)
    return xAcceleration,yAcceleration,zAcceleration
    
def CalculateNewVelocitiesAndPositionsUsing(currentPositions,currentVelocities,deltaT,masses):
    xAcceleration,yAcceleration,zAcceleration=calculateAccelerationsOnEacyBodyUsingThe(currentPositions,masses)
    connectivityMatrix=CalculateConnectivityMatrix(len(currentPositions))
    xChangeInVelocities=deltaT*np.matmul(xAcceleration,connectivityMatrix)
    yChangeInVelocities=deltaT*np.matmul(yAcceleration,connectivityMatrix)
    zChangeInVelocities=deltaT*np.matmul(zAcceleration,connectivityMatrix)
    currentVelocities[:,0]+=xChangeInVelocities
    currentVelocities[:,1]+=yChangeInVelocities
    currentVelocities[:,2]+=zChangeInVelocities
    changeInPositions=deltaT*currentVelocities
    currentPositions+=changeInPositions
    return currentVelocities,currentPositions

def CalculateApproximateNewVelocitiesAndPositionsUsing(currentPositions,currentVelocities,deltaT,masses,accelerationTreshold,Frequency,timestep,overallAcceleration):
    xAcceleration,yAcceleration,zAcceleration=calculateAccelerationsOnEacyBodyUsingThe(currentPositions,masses)
    if timestep % Frequency ==1:
        overallAcceleration=np.sqrt(np.square(xAcceleration)+np.square(yAcceleration)+np.square(zAcceleration))
    xAcceleration[overallAcceleration<accelerationTreshold]=0
    yAcceleration[overallAcceleration<accelerationTreshold]=0
    zAcceleration[overallAcceleration<accelerationTreshold]=0
    xAcceleration=sparse.csr_matrix(xAcceleration)
    yAcceleration=sparse.csr_matrix(yAcceleration)
    zAcceleration=sparse.csr_matrix(zAcceleration)
    connectivityMatrix=CalculateConnectivityMatrix(len(currentPositions))
    xChangeInVelocities=deltaT*(xAcceleration*connectivityMatrix)
    yChangeInVelocities=deltaT*(yAcceleration*connectivityMatrix)
    zChangeInVelocities=deltaT*(zAcceleration*connectivityMatrix)
    currentVelocities[:,0]+=xChangeInVelocities
    currentVelocities[:,1]+=yChangeInVelocities
    currentVelocities[:,2]+=zChangeInVelocities
    changeInPositions=deltaT*currentVelocities
    currentPositions+=changeInPositions
    return currentVelocities,currentPositions,overallAcceleration

def CalculateConnectivityMatrix(nBodies):
    return np.ones((nBodies), dtype=int)

def executeSimulationAssuming(nBodies,deltaT,nTimesteps):
    currentPositions=createRandomInitialPositionsFor(nBodies)
    currentVelocities=createRandomInitialVelocitiesFor(nBodies)
    masses=createRandomMassesFor(nBodies)
    i=0
    while True:
        i=i+1
        currentVelocities,currentPositions=CalculateNewVelocitiesAndPositionsUsing(currentPositions,currentVelocities,deltaT,masses)
        if i==nTimesteps:
            break
def executeApproximateSimulationAssuming(nBodies,deltaT,nTimesteps,accelerationTreshold,Frequency):
    currentPositions=createRandomInitialPositionsFor(nBodies)
    currentVelocities=createRandomInitialVelocitiesFor(nBodies)
    masses=createRandomMassesFor(nBodies)
    currentTimestep=0
    overallAcceleration=np.ones([nBodies, nBodies])
    while True:
        currentTimestep=currentTimestep+1
        currentVelocities,currentPositions,overallAcceleration=CalculateApproximateNewVelocitiesAndPositionsUsing(currentPositions,currentVelocities,deltaT,masses,accelerationTreshold,Frequency,currentTimestep,overallAcceleration)
        if currentTimestep==nTimesteps:
            break

#xAcceleration,yAcceleration,zAcceleration=calculateAccelerationsOnEacyBodyUsingThe(createRandomInitialPositionsFor(100),createRandomMassesFor(100))
#print(np.percentile((np.sqrt(np.square(xAcceleration)+np.square(yAcceleration)+np.square(zAcceleration))), 25))
#print(np.percentile((np.sqrt(np.square(xAcceleration)+np.square(yAcceleration)+np.square(zAcceleration))), 75))