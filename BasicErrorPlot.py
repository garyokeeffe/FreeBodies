#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 15:37:15 2018

"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import copy


# Setting up the system parameter values
n_freebodies=10
n_dimensions=3
n_timesteps=10000
G=6.67408*(10**(-11))

delta_t=1*60;
masses=np.array([[5*(10 ** 7)],
[6.1*(10 ** 3)],
[5.5*(10 ** 4)],
[7.6*(10 ** 4)],
[3.4*(10 ** 4)],
[9.9*(10 ** 4)],
[1*(10 ** 4)],
[7.2*(10 ** 4)],
[6.7*(10 ** 4)],
[1.7*(10 ** 4)]])
#masses[2]=2000000
#masses[1]=100000
masses=np.tile(masses, [1,n_freebodies])
masses=masses.transpose()

velocities=1.5/delta_t*(-np.random.rand(n_freebodies,n_dimensions))

velocities[0,:]=[0,0,0]

velocities=np.array([[0.0000,	0.0000,0.0000],
[-0.0003	,-0.0007	,-0.0015],
[0.0009	,-0.0003	,0.0012],
[-0.0013	,0.0002	,-0.0004],
[0.0014	,0.0009	,0.0008],
[-0.0002	,-0.0011	,-0.0005],
[0.0002	,-0.0012	,0.0013],
[-0.0004	,0.0016	,-0.0003],
[0.0007	,0.0013	,0.0010],
[-0.0009	,-0.0003,-0.0004]])
velocities=2*2*velocities


positions=np.array([[0.00,0.00,0.00],
[-71.266	,85.619,38.380],
[-134.393,-34.303,17.161],
[30.281,-101.933	,28.618],
[-20.780,-10.504	,-22.029],
[-51.725	,-94.240	,59.237],
[86.451,31.849,-20.616],
[154.503,-98.275,-.468],
[40.681,15.197,40.934],
[31.796,-17.506,.165]])

positions1= copy.deepcopy(positions)
positions2= copy.deepcopy(positions)
positions3= copy.deepcopy(positions)

#velocities[:,:,0]=0.01/delta_t*(1000-positions[:,:,0])
velocities1=copy.deepcopy(velocities)
velocities2=copy.deepcopy(velocities)
velocities3=copy.deepcopy(velocities)
i=0
loop=True
fig = plt.figure()
fig.patch.set_facecolor('xkcd:black')
ax = fig.add_subplot(111, projection='3d')
imax=1000
error=np.zeros(imax)
error2=np.zeros(imax)
error3=np.zeros(imax)
while loop==True:
    i=i+1
    #loop=False
    previousposition=positions
    previousvelocity=velocities
    previousposition1=positions1
    previousvelocity1=velocities1
    previousposition2=positions2
    previousvelocity2=velocities2
    previousposition3=positions3
    previousvelocity3=velocities3
    
    posx=np.tile(positions[:,0], [n_freebodies,1])
    posy=np.tile(positions[:,1], [n_freebodies,1])
    posz=np.tile(positions[:,2], [n_freebodies,1])
    posx1=np.tile(positions1[:,0], [n_freebodies,1])
    posy1=np.tile(positions1[:,1], [n_freebodies,1])
    posz1=np.tile(positions1[:,2], [n_freebodies,1])
    posx2=np.tile(positions2[:,0], [n_freebodies,1])
    posy2=np.tile(positions2[:,1], [n_freebodies,1])
    posz2=np.tile(positions2[:,2], [n_freebodies,1])
    posx3=np.tile(positions3[:,0], [n_freebodies,1])
    posy3=np.tile(positions3[:,1], [n_freebodies,1])
    posz3=np.tile(positions3[:,2], [n_freebodies,1])

    edgelengthx=posx-posx.transpose()
    edgelengthy=posy-posy.transpose()
    edgelengthz=posz-posz.transpose()
    edgelengthx1=posx1-posx1.transpose()
    edgelengthy1=posy1-posy1.transpose()
    edgelengthz1=posz1-posz1.transpose()
    edgelengthx2=posx2-posx2.transpose()
    edgelengthy2=posy2-posy2.transpose()
    edgelengthz2=posz2-posz2.transpose()
    edgelengthx3=posx3-posx3.transpose()
    edgelengthy3=posy3-posy3.transpose()
    edgelengthz3=posz3-posz3.transpose()
    
    edgesquares=np.array(np.square(edgelengthx)+np.square(edgelengthy)+np.square(edgelengthz))
    edgelengths=np.array(np.sqrt(edgesquares))
    edgesquares1=np.array(np.square(edgelengthx1)+np.square(edgelengthy1)+np.square(edgelengthz1))
    edgelengths1=np.array(np.sqrt(edgesquares1))
    edgesquares2=np.array(np.square(edgelengthx2)+np.square(edgelengthy2)+np.square(edgelengthz2))
    edgelengths2=np.array(np.sqrt(edgesquares2))
    edgesquares3=np.array(np.square(edgelengthx3)+np.square(edgelengthy3)+np.square(edgelengthz3))
    edgelengths3=np.array(np.sqrt(edgesquares3))
    
    Qx=np.divide(edgelengthx,(edgesquares*edgelengths))
    Qy=np.divide(edgelengthy,(edgesquares*edgelengths))
    Qz=np.divide(edgelengthz,(edgesquares*edgelengths))
    Qx1=np.divide(edgelengthx1,(edgesquares1*edgelengths1))
    Qy1=np.divide(edgelengthy1,(edgesquares1*edgelengths1))
    Qz1=np.divide(edgelengthz1,(edgesquares1*edgelengths1))
    Qx2=np.divide(edgelengthx2,(edgesquares2*edgelengths2))
    Qy2=np.divide(edgelengthy2,(edgesquares2*edgelengths2))
    Qz2=np.divide(edgelengthz2,(edgesquares2*edgelengths2))
    Qx3=np.divide(edgelengthx3,(edgesquares3*edgelengths3))
    Qy3=np.divide(edgelengthy3,(edgesquares3*edgelengths3))
    Qz3=np.divide(edgelengthz3,(edgesquares3*edgelengths3))
    
    Fx=G*Qx
    Fy=G*Qy
    Fz=G*Qz
    
    Fx1=G*Qx1
    Fy1=G*Qy1
    Fz1=G*Qz1
    
    Fx2=G*Qx2
    Fy2=G*Qy2
    Fz2=G*Qz2
    
    Fx3=G*Qx3
    Fy3=G*Qy3
    Fz3=G*Qz3
    
    where_are_NaNsx = np.isnan(Fx)
    where_are_NaNsy = np.isnan(Fy)
    where_are_NaNsz = np.isnan(Fz)
    Fx[where_are_NaNsx] = 0
    Fy[where_are_NaNsy] = 0
    Fz[where_are_NaNsz] = 0
    where_are_NaNsx1 = np.isnan(Fx1)
    where_are_NaNsy1 = np.isnan(Fy1)
    where_are_NaNsz1 = np.isnan(Fz1)
    Fx1[where_are_NaNsx1] = 0
    Fy1[where_are_NaNsy1] = 0
    Fz1[where_are_NaNsz1] = 0
    where_are_NaNsx2 = np.isnan(Fx2)
    where_are_NaNsy2 = np.isnan(Fy2)
    where_are_NaNsz2 = np.isnan(Fz2)
    Fx2[where_are_NaNsx2] = 0
    Fy2[where_are_NaNsy2] = 0
    Fz2[where_are_NaNsz2] = 0
    where_are_NaNsx3 = np.isnan(Fx3)
    where_are_NaNsy3 = np.isnan(Fy3)
    where_are_NaNsz3 = np.isnan(Fz3)
    Fx3[where_are_NaNsx3] = 0
    Fy3[where_are_NaNsy3] = 0
    Fz3[where_are_NaNsz3] = 0
    
    Fx=np.multiply(Fx,masses)
    Fy=np.multiply(Fy,masses)
    Fz=np.multiply(Fz,masses)
    Fx1=np.multiply(Fx1,masses)
    Fy1=np.multiply(Fy1,masses)
    Fz1=np.multiply(Fz1,masses)
    Fx2=np.multiply(Fx2,masses)
    Fy2=np.multiply(Fy2,masses)
    Fz2=np.multiply(Fz2,masses)
    Fx3=np.multiply(Fx3,masses)
    Fy3=np.multiply(Fy3,masses)
    Fz3=np.multiply(Fz3,masses)

    FXZZ=np.sqrt(np.square(Fx1)+np.square(Fy1)+np.square(Fz1))
    if i % 10 ==1:
        FXZZ2=np.sqrt(np.square(Fx2)+np.square(Fy2)+np.square(Fz2))
    if i % 100 ==1:
        FXZZ3=np.sqrt(np.square(Fx3)+np.square(Fy3)+np.square(Fz3))
        
    Fx1[FXZZ<.00000001]=0
    Fy1[FXZZ<.00000001]=0
    Fz1[FXZZ<.00000001]=0
    
    Fx2[FXZZ2<.00000001]=0
    Fy2[FXZZ2<.00000001]=0
    Fz2[FXZZ2<.00000001]=0
    
    Fx3[FXZZ3<.00000001]=0
    Fy3[FXZZ3<.00000001]=0
    Fz3[FXZZ3<.00000001]=0
    
    W=np.ones((n_freebodies), dtype=int)

    velocities[:,0]=previousvelocity[:,0]+delta_t*np.matmul(Fx,W)
    velocities[:,1]=previousvelocity[:,1]+delta_t*np.matmul(Fy,W)
    velocities[:,2]=previousvelocity[:,2]+delta_t*np.matmul(Fz,W)
    velocities1[:,0]=previousvelocity1[:,0]+delta_t*np.matmul(Fx1,W)
    velocities1[:,1]=previousvelocity1[:,1]+delta_t*np.matmul(Fy1,W)
    velocities1[:,2]=previousvelocity1[:,2]+delta_t*np.matmul(Fz1,W)  
    velocities2[:,0]=previousvelocity2[:,0]+delta_t*np.matmul(Fx2,W)
    velocities2[:,1]=previousvelocity2[:,1]+delta_t*np.matmul(Fy2,W)
    velocities2[:,2]=previousvelocity2[:,2]+delta_t*np.matmul(Fz2,W)
    velocities3[:,0]=previousvelocity3[:,0]+delta_t*np.matmul(Fx3,W)
    velocities3[:,1]=previousvelocity3[:,1]+delta_t*np.matmul(Fy3,W)
    velocities3[:,2]=previousvelocity3[:,2]+delta_t*np.matmul(Fz3,W)

    positions[:,0]=previousposition[:,0]+delta_t*velocities[:,0]
    positions[:,1]=previousposition[:,1]+delta_t*velocities[:,1]
    positions[:,2]=previousposition[:,2]+delta_t*velocities[:,2]
    positions1[:,0]=previousposition1[:,0]+delta_t*velocities1[:,0]
    positions1[:,1]=previousposition1[:,1]+delta_t*velocities1[:,1]
    positions1[:,2]=previousposition1[:,2]+delta_t*velocities1[:,2]
    positions2[:,0]=previousposition2[:,0]+delta_t*velocities2[:,0]
    positions2[:,1]=previousposition2[:,1]+delta_t*velocities2[:,1]
    positions2[:,2]=previousposition2[:,2]+delta_t*velocities2[:,2]
    positions3[:,0]=previousposition3[:,0]+delta_t*velocities3[:,0]
    positions3[:,1]=previousposition3[:,1]+delta_t*velocities3[:,1]
    positions3[:,2]=previousposition3[:,2]+delta_t*velocities3[:,2]
    error[i-1]=sum(np.sqrt(np.square(positions1[:,0]-positions[:,0])+np.square(positions1[:,1]-positions[:,1])+np.square(positions1[:,2]-positions[:,2])))/n_freebodies
    error2[i-1]=sum(np.sqrt(np.square(positions2[:,0]-positions[:,0])+np.square(positions2[:,1]-positions[:,1])+np.square(positions2[:,2]-positions[:,2])))/n_freebodies
    error3[i-1]=sum(np.sqrt(np.square(positions3[:,0]-positions[:,0])+np.square(positions3[:,1]-positions[:,1])+np.square(positions3[:,2]-positions[:,2])))/n_freebodies
    ax.clear()
    ax.set_facecolor('xkcd:black')
    ax.scatter(positions[:,0], positions[:,1], positions[:,2], s=masses**(1/3), c='w',marker="o")
    ax.scatter(positions1[:,0], positions1[:,1], positions1[:,2], s=masses**(1/3), c='y',marker="o")
    ax.scatter(positions2[:,0], positions2[:,1], positions2[:,2], s=masses**(1/3), c='g',marker="o")
    ax.scatter(positions3[:,0], positions3[:,1], positions3[:,2], s=masses**(1/3), c='b',marker="o")
    for j in range(0, n_freebodies):
        ax.plot3D([positions[j,0],positions1[j,0]],[positions[j,1],positions1[j,1]],[positions[j,2],positions1[j,2]], 'r-')
        ax.plot3D([positions[j,0],positions2[j,0]],[positions[j,1],positions2[j,1]],[positions[j,2],positions2[j,2]], 'm-')
        ax.plot3D([positions[j,0],positions3[j,0]],[positions[j,1],positions3[j,1]],[positions[j,2],positions3[j,2]], 'p-')
    ax.set_xlim((positions[0,0]+positions1[0,0])/2-200, (positions[0,0]+positions1[0,0])/2+200)
    ax.set_ylim((positions[0,1]+positions1[0,1])/2-200, (positions[0,1]+positions1[0,1])/2+200)
    ax.set_zlim((positions[0,2]+positions1[0,2])/2-200, (positions[0,2]+positions1[0,2])/2+200)
    plt.axis('off')
    plt.pause(0.001)
    if i==imax:
        loop=False
print("The final error is", error[i-1])
plt.figure()
plt.plot(error, 'k-')
plt.plot(error2, 'k--')
plt.plot(error3, 'k-.')
#plt.set_xlim(0,imax)
plt.xlabel('Timesteps')
plt.ylabel('Error (m)')

plt.show()
