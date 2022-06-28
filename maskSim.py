# -*- coding: utf-8 -*-
"""
Created on Fri May  6 10:32:13 2022

@author: tbart
"""

import numpy as np
import cmath
from cmath import phase
import math
from math import pi
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib import image
from PIL import ImageOps, Image
from numpy import asarray
import random

"""

First experiment:
    Do a single pixel correlation

"""
plt.close('all')
#2 photon beams, with some intensity distribution. In this case
#its just a 1, because its only 1 photon
#beamA = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
#beamB = [0.5, 0.4, 0.3, 0.2, 0.1, 0]



#Process:
# Model light as gaussian array. Left half is beamA, right half is beamB
# Create SLM phase reflection from image (normalize 0-1)

#testImage = Image.open("byu_test.png")
testImage = Image.open("byu_test.png")
testImage = ImageOps.grayscale(testImage)
imageArray = asarray(testImage)

#numRows must be even, so I can split the beam in half
maskRows = imageArray.shape[0]
maskCols = imageArray.shape[1]



#apply beam A to an slm with a constant phase shift of
A = 1 # magnitude of slm complex 

# Initializing value of x-axis and y-axis
# in the range -2 to +2
gx, gy = np.meshgrid(np.linspace(-1,1,maskCols*2), np.linspace(-1,1,maskRows))
dst = np.sqrt(gx*gx+gy*gy)
 
# Initializing sigma and muu
sigma = 1
muu = 0.000
 
# Calculating Gaussian array
gauss = np.exp(-( (dst-muu)**2 / ( 2.0 * sigma**2 ) ) ) #models light from laser
plt.imshow( gauss, cmap = cm.gray)
plt.show() #shows gaussian

beamA = np.hsplit(gauss,2)[0] # left side of gaussian array
beamB = np.hsplit(gauss,2)[1] # right side of gaussian array


#phase_mask =  [[1 for i in range(maskCols)] for j in range(maskRows)] # maskCols x maskRows array of 1's

#pm_test = np.divide(pi,1+np.arange(maskRows*maskCols).reshape((maskRows,maskCols))) # array, each entry is pi divided by 1 through maskRows*maskCols
#pm_test2 = np.multiply(1/pi, np.random.randint(-0, high=8, size=(maskRows,maskCols))) #array, nums 0-7 divided by pi (0 - 7/pi)
pm_test3 = np.multiply(imageArray, 1/255) #scaling 0-1 instead of 0-255
phi_mask = pm_test3

slm_mask = [[0 for i in range(maskCols)] for j in range(maskRows)]
phaseOfSlm = [[0 for i in range(maskCols)] for j in range(maskRows)]

for i in range(maskRows):
    for k in range(maskCols):
        pm = phi_mask[i,k]
        slm_mask[i][k] = complex(A, pm)
        phaseOfSlm[i][k] = phase(complex(A, pm))
        #print(str(i) + ", " + str(k) + ": " + str(slm_mask[i][k]) + ", " + str(phaseOfSlm[i][k]))
    
B_post_slm = np.multiply(beamB, slm_mask) # B portion of laser after SLM reflection
#B_post_slm = slm_mask #because beam intensity is 1

numPhi0 = 4
phi0 = np.multiply(math.pi,[0, 1/2, 1, 3/2]) # 4 phases

row_data = [[[0 for i in range(maskCols)] for j in range(maskRows)] for k in range(numPhi0)] # room for image at each phase rotation
for p0 in range(len(phi0)): # repeats for 0, pi/2, pi, 3pi/2
    
    slm = np.multiply(np.ones((maskRows, maskCols)), complex(A, phi0[p0])) #full image of each phase
    #print(complex(A, phi0[p0]))
    #print(phase(complex(A, phi0[p0])) * (180 / (2*math.pi)))

    A_post_slm = np.multiply(beamA, slm) # A side after slm with specific phase
    #A_post_slm = slm #because beam intensity is 1
    
    
        
    for i in range(maskRows):
        for k in range(maskCols):
            phasea = phase(A_post_slm[i][k]) 
            phaseb = phase(B_post_slm[i][k])
            rd = 0.5 * (1 + math.cos(phaseb + phasea)) #find equation to remember what it's for
            row_data[p0][i][k] = rd #saving calculation for each phase

correlation = [[0 for i in range(maskCols)] for j in range(maskRows)]
for i in range(maskRows):
    for k in range(maskCols):
        
        #simulate amplification noise with a random multiplier
        ampNoiseLow = 0.98
        ampNoiseHigh = 1.02        
        
        real = row_data[0][i][k] - row_data[2][i][k] # 0 - pi
        imag = row_data[1][i][k] - row_data[3][i][k] # pi/2 - 3pi/2
        real = real * random.uniform(ampNoiseLow, ampNoiseHigh)
        imag = imag * random.uniform(ampNoiseLow, ampNoiseHigh)
        c = complex(real, imag)
        correlation[i][k] = phase(c) #saves resulting phase
    
axes=[]
fig = plt.figure()

axes.append(fig.add_subplot(1, 3, 1))
plt.imshow( correlation, cmap = cm.gray)
subplot_title=("Reconstructed")
axes[-1].set_title(subplot_title)  


axes.append(fig.add_subplot(1, 3, 2))
plt.imshow( phaseOfSlm, cmap = cm.gray)
subplot_title=("Reference Image")
axes[-1].set_title(subplot_title)  

fig.tight_layout() 
plt.show()



#plt.plot(phaseOfSlm)
#factor = np.divide(correlation, phaseOfSlm)
    



