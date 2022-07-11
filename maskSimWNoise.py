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
from random import choices

def NormalizeData(data): #returns array normalized 0-1 (0 is min, 1 is max)
    return (data - np.min(data)) / (np.max(data) - np.min(data))


plt.close('all') #reseting plots

#2 photon beams, with some intensity distribution. In this case
#its just a 1, because its only 1 photon
#beamA = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
#beamB = [0.5, 0.4, 0.3, 0.2, 0.1, 0]

testImage = Image.open("diagonal_test2.png")
testImage = ImageOps.grayscale(testImage)
imageArray = asarray(testImage)

#numRows must be even, so I can split the beam in half
maskRows = imageArray.shape[0] #diagonal_test2 = 10
maskCols = imageArray.shape[1] #diagonal_test2 = 5
#print("Rows: " + str(maskRows) + "\nCols: " + str(maskCols))


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
gauss = np.exp(-( (dst-muu)**2 / ( 2.0 * sigma**2 ) ) )
#plt.imshow( gauss, cmap = cm.gray)
#plt.show()

beamA = np.hsplit(gauss,2)[0]
beamB = np.hsplit(gauss,2)[1]


#phase_mask =  [[1 for i in range(maskCols)] for j in range(maskRows)]

#pm_test = np.divide(pi,1+np.arange(maskRows*maskCols).reshape((maskRows,maskCols))) 
#pm_test2 = np.multiply(1/pi, np.random.randint(-0, high=8, size=(maskRows,maskCols)))
pm_test3 = np.multiply(imageArray, 1/255)
phi_mask = pm_test3

slm_mask = [[0 for i in range(maskCols)] for j in range(maskRows)]
phaseOfSlm = [[0 for i in range(maskCols)] for j in range(maskRows)]

for i in range(maskRows):
    for k in range(maskCols):
        pm = phi_mask[i,k]
        slm_mask[i][k] = complex(A, pm)
        phaseOfSlm[i][k] = phase(complex(A, pm))

    
#B_post_slm = np.multiply(beamB, slm_mask) #gaussian array intensity
B_post_slm = slm_mask #because beam intensity is 1

numPhi0 = 4
phi0 = np.multiply(math.pi,[0, 1/2, 1, 3/2])

row_data = [[[0 for i in range(maskCols)] for j in range(maskRows)] for k in range(numPhi0)]
row_counts = [[[0 for i in range(maskCols)] for j in range(maskRows)] for k in range(numPhi0)]
row_counts_100 = [[[0 for i in range(maskCols)] for j in range(maskRows)] for k in range(numPhi0)]

beamDatai = []
beamDatak = []


for p0 in range(len(phi0)):
    
    slm = np.multiply(np.ones((maskRows, maskCols)), complex(A, phi0[p0]))
    #A_post_slm = np.multiply(beamA, slm)
    A_post_slm = slm #because beam intensity is 1
    
    #Run a bunch photon pais per reference phase shift
    numPhotons = 50000
    cPhoton = 1    

    while cPhoton < numPhotons:

        
        #Randomly Generate a pair of photons with complimentary coordinates
        normStDev = min(maskRows, maskCols) / 1
        randIa = np.random.normal(maskRows/2, normStDev, 1)
        randI = round(randIa[0])
        beamDatai.append(randI)
        randKa = np.random.normal(0, normStDev, 1)
        randK = abs(round(randKa[0]))
        beamDatak.append(randK)
        
        
        
        ##above code does a normal distribition 
        ##below code does a uniform distribution
        
        randI = random.randint(0, maskRows-1)
        randK = random.randint(0, maskCols-1)
        
        if randI >= maskRows or randK >= maskCols:
            cPhoton = cPhoton + 1
            continue
        if randI < 0 or randK < 0:
            cPhoton = cPhoton + 1
            continue
       
       
        #calculate the phase shift based on the coords and mask info
        #this allows me to calculate R, the probability of getting a click
        phasea = phase(A_post_slm[randI][randK]) 
        phaseb = phase(B_post_slm[randI][randK])
        rd = 0.5 * (1 + math.cos(phaseb + phasea))
        if(row_data[p0][randI][randK] == 0):
            row_data[p0][randI][randK] = rd
        
        #detector efficiency is about 48%
        deff = 0.48
        rd_x_deff = rd * deff
        population = [0,1]
        weights = [rd, 1-rd]
        clickOrNot = choices(population,weights)
        
         
        #create false coincidence x% of the time
        false_coin_rate = 0.005 #1 %
        population = [0,1]
        weights = [1-false_coin_rate, false_coin_rate]
        fc_or_not = choices(population,weights)
        
        if fc_or_not[0] == 1:
            row_counts[p0][randI][randK] = row_counts[p0][randI][randK] + 1
        elif clickOrNot[0] == 1:
            row_counts[p0][randI][randK] = row_counts[p0][randI][randK] + 1
        
        
        row_counts_100[p0][randI][randK] = row_counts_100[p0][randI][randK] + 1
        
        ##above code does a uniform distribution

        cPhoton = cPhoton + 1
                

#print(row_counts[0][randI][randK])
#print(row_counts_100[0][randI][randK])

#plt.hist(beamDatai)
#plt.figure()
#plt.hist(beamDatak)

#extract R (The probability that a photon pair will click given a phase shift)
#from how many photons were aimed at a cell / how many clicke
#in reality we can't know this but we may be able to estimate it
row_data = np.array(row_data)
extractedR = np.divide(row_counts, row_counts_100)

probDif = np.subtract(extractedR, row_data)

correlation = [[0 for i in range(maskCols)] for j in range(maskRows)]
correlation_recon = [[0 for i in range(maskCols)] for j in range(maskRows)]
for i in range(maskRows):
    for k in range(maskCols):
        real = row_data[0][i][k] - row_data[2][i][k]
        imag = row_data[1][i][k] - row_data[3][i][k]
        
        #simulate amplification noise with a random multiplier
        ampNoiseLow = 0.999
        ampNoiseHigh = 1.001
        
        real = real * random.uniform(ampNoiseLow, ampNoiseHigh)
        imag = imag * random.uniform(ampNoiseLow, ampNoiseHigh)
        c = complex(real, imag)
        correlation[i][k] = 2*math.pi -phase(c)
        
        real_recon = extractedR[0][i][k] - extractedR[2][i][k]
        imag_recon = extractedR[1][i][k] - extractedR[3][i][k]
        
        #divide by 100 then mulitply by not quite 100 to simulate amplification noise
        real_recon = real_recon * random.uniform(ampNoiseLow, ampNoiseHigh)
        imag_recon = imag_recon * random.uniform(ampNoiseLow, ampNoiseHigh)
        c_recon = complex(real_recon, imag_recon)
        correlation_recon[i][k] = -phase(c_recon)


correlation = np.array(correlation)
correlation_recon = np.array(correlation_recon)

cr_norm = NormalizeData(correlation_recon)
cr_norm_idea = NormalizeData(correlation)
im_norm = phi_mask
    
axes=[]
fig = plt.figure()

#try to calc SNR of reconstruction



axes.append(fig.add_subplot(2, 3, 1))
plt.imshow( cr_norm_idea, cmap = cm.gray)
subplot_title=("Reconstructed")
axes[-1].set_title(subplot_title)  


axes.append(fig.add_subplot(2, 3, 2))
plt.imshow( phaseOfSlm, cmap = cm.gray)
subplot_title=("Reference Image")
axes[-1].set_title(subplot_title)  

axes.append(fig.add_subplot(2, 3, 3))
plt.imshow( row_counts[0], cmap = cm.gray)
subplot_title=("Clicks per cell")
axes[-1].set_title(subplot_title)  

axes.append(fig.add_subplot(2, 3, 4))
plt.imshow( cr_norm, cmap = cm.gray)
subplot_title=("Reconstructed 2")
axes[-1].set_title(subplot_title)  


axes.append(fig.add_subplot(2, 3, 5))
plt.imshow( phaseOfSlm, cmap = cm.gray)
subplot_title=("Reference Image")
axes[-1].set_title(subplot_title)  

axes.append(fig.add_subplot(2, 3, 6))
plt.imshow( row_counts_100[0], cmap = cm.gray)
subplot_title=("Total photons per cell")
axes[-1].set_title(subplot_title)  

fig.tight_layout() 
plt.show()



#plt.plot(phaseOfSlm)
#factor = np.divide(correlation, phaseOfSlm)
    



