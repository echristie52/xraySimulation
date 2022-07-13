# Eric Christie
# BYU, 6/23/22
# Simulation of SLM and Image Reconstruction
# based on paper, Polarization entanglement-enabled quantum holography
# slmSim uses equation 6 to theoretically determine R(k) for reconstruction
# V1: no variation of intensity and initial phase change is 0 (psi0)


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

#can plot any number of arrays with titles, scaled booleans select if graphed 0-2pi or not
def plot(arrays, titles, scaled):
    axes=[]
    fig = plt.figure()
    
    for plotNum in range(len(arrays)): #repeats for each array given
        axes.append(fig.add_subplot(1, len(arrays), plotNum+1))

        if(scaled[plotNum]): #bool to see if scaled or not
            plt.imshow(arrays[plotNum], cmap = cm.twilight, vmin = 0, vmax = 2*math.pi)
            
        else:   
            plt.imshow(arrays[plotNum], cmap = cm.gray)
        plt.colorbar(shrink = 0.5)
        
        if plotNum < len(titles):
            subplot_title = titles[plotNum]
            axes[-1].set_title(subplot_title)

    fig.tight_layout()
    plt.show()  #outputs plots


###### Begin Simulation ######
plt.close('all')


phase0 = 0.0

ampNoiseLow = 0.98 #amplifier noise values
ampNoiseHigh = 1.02 

#open SLM image
testImage = Image.open("byu_test.png")
testImage = ImageOps.grayscale(testImage)
imageArray = asarray(testImage)

maskRows = imageArray.shape[0] 
maskCols = imageArray.shape[1]

slm = np.multiply(imageArray, 2*math.pi/255.0) # scales image to 0-2pi phase changes

# Make array of photons (dimensions for both beams)
#photons = [[Photon(magnitude, phase0) for i in range(maskRows)] for j in range(2*maskCols)]
photons = []
for i in range(maskRows):
    photons.append([])
    for j in range(2 * maskCols):
        photons[i].append(phase0) #can replace these with more accurate initial conditions
photons = np.asarray(photons)

beamA = np.hsplit(photons,2)[0] #left half of photons, reflect Image
beamB = np.hsplit(photons,2)[1] #right half of photons, reflect phase shifts



# Photon reflection on SLM and Reconstruction
# Does each position k/-k one at a time to see full calculation

A_post_slm = beamA #image
phases = np.multiply(math.pi, [0, 1/2, 1, 3/2])
numPhases = 4
reconstruction = [[0 for i in range(maskCols)] for j in range(maskRows)] #results from intensity reconstruction (w/o coincidence counting)
accuracy = [[0 for i in range(maskCols)] for j in range(maskRows)]

for i in range(maskRows): # i,j are position of k/-k photons
    for j in range(maskCols): 
        #phaseChange from slm to beamA (reflection onto A_post_slm)
        A_post_slm[i][j] = (beamA[i][j] + slm[i][j]) / 1.0 #adds phase change onto beamA photon 

        #Bob phase changes
        phaseIntensity = [0 for p in range(numPhases)] #will save calculated values
        for p in range(numPhases): #goes through all 4 phase shifts and calculates R's
            #phaseA = beamA + slmzz
            #phaseB = beamB + phases
            phaseIntensity[p] = 0.5 * (1 + math.cos(beamB[i][j] + phases[p] + A_post_slm[i][j])) #takes inital beamB, adds this phase change, and then A
            
        
        #do reconstruction math with noise
        real = phaseIntensity[0] - phaseIntensity[2] # values for Eq 1
        imag = phaseIntensity[1] - phaseIntensity[3]
        real = real * random.uniform(ampNoiseLow, ampNoiseHigh) #low-level detector noise
        imag = imag * random.uniform(ampNoiseLow, ampNoiseHigh)
        c = phase(complex(real, imag)) # extracts phase
        
        if c < 0: #wraps negative phase shift to positive
            c += 2*math.pi

        reconstruction[i][j] = 2*math.pi - c #(2pi - results) yields original phase shift
        
        accuracy[i][j] = np.abs(reconstruction[i][j] - slm[i][j]) # records difference for comparison
 

# Begin Plotting of Results
plot([reconstruction, slm, accuracy], ["Reconstructed", "Reference Image", "Accuracy"], [True, True, False])
        
        

