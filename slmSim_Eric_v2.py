# Eric Christie
# BYU, 6/23/22
# Simulation of SLM and Image Reconstruction
# based on paper, Polarization entanglement-enabled quantum holography
# slmSim uses equation 6 to theoretically determine R(k) for reconstruction
# V2: Photon object different intensity and phase based on gaussian and paper distributions


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

#Constants
ampNoiseLow = 0.98 #amplifier noise values
ampNoiseHigh = 1.02 

class Photon:
    def __init__(self, intensity, phase): #initialized Photon
        self.intensity = intensity
        self.phase = phase
    
    def setIntensity(self, intensity): #setters/getters
        self.intensity = intensity

    def getIntensity(self):
        return self.intensity
    
    def setPhase(self, phase):
        self.phase = phase

    def getPhase(self):
        return self.phase

    def changePhase(self, phase):
        self.phase += phase
    
# scales an array to a new max/min
def scaleArray(array, newMax, newMin):
    max = np.max(array)
    min = np.min(array)
    newArray = array - min # newVal = (oldVal - min) * (newRange / oldRange) + newMin
    newArray = np.multiply((newMax-newMin)/(max-min), newArray)
    newArray += newMin
    return newArray



###### Begin Simulation ######
plt.close('all')

#open SLM image
testImage = Image.open("byu_test.png")
testImage = ImageOps.grayscale(testImage)
imageArray = asarray(testImage)

maskRows = imageArray.shape[0] #numRows must be even, so I can split the beam in half
maskCols = imageArray.shape[1]

slm = np.multiply(imageArray, 2*math.pi/255.0) # scales image to 0-2pi phase changes

#Gaussian array to simulate laser intensity
# Initializing value of x-axis and y-axis
# in the range -2 to +2
gx, gy = np.meshgrid(np.linspace(-1,1,maskCols*2), np.linspace(-1,1,maskRows))
dst = np.sqrt(gx*gx+gy*gy)
 
sigma = 1 # Initializing sigma and muu
muu = 0.000
 
gauss = np.exp(-( (dst-muu)**2 / ( 2.0 * sigma**2 ) ) ) # Calculating Gaussian array

#scale gauss for phase distortion
psi0_correction = scaleArray(gauss, 1, -8) #scales gauss to psi0_correction mask

psi0 = np.multiply(-1, psi0_correction) # correction is opposite of initial psi0

# Plotting of incoming laser characteristics
axes=[]
fig = plt.figure()

axes.append(fig.add_subplot(1, 3, 1)) # reconstruction
plt.imshow( gauss, cmap = cm.gray)
subplot_title=("Intensity")
axes[-1].set_title(subplot_title)  

axes.append(fig.add_subplot(1, 3, 2)) #image
plt.imshow( psi0, cmap = cm.gray)
subplot_title=("Phase")
axes[-1].set_title(subplot_title)  

axes.append(fig.add_subplot(1, 3, 3)) #image
plt.imshow( psi0_correction, cmap = cm.gray)
subplot_title=("Phase Correction")
axes[-1].set_title(subplot_title) 

fig.tight_layout() 
plt.show()




# Make array of photons (dimensions for both beams)
photons = []
for i in range(maskRows):
    photons.append([])
    for j in range(2 * maskCols):
        photons[i].append(Photon(gauss[i][j], (psi0[i][j]))) #adds photon based on gaussian intensity and scaled phase distortion
photons = np.asarray(photons)

beamA = np.hsplit(photons,2)[0] #left half of photons, reflect Image
beamB = np.hsplit(photons,2)[1] #right half of photons, reflect phase shifts

slm_B_mask = np.hsplit(psi0_correction,2)[1]



# Photon reflection on SLM and Reconstruction
# Does each position k/-k one at a time to see full calculation
phases = np.multiply(math.pi, [0, 1/2, 1, 3/2])
numPhases = 4
reconstruction = [[0 for i in range(maskCols)] for j in range(maskRows)] #results from intensity reconstruction (w/o coincidence counting)
accuracy = [[0 for i in range(maskCols)] for j in range(maskRows)] #compares reults to original phase

for i in range(maskRows): # i,j are position of k/-k photons
    for j in range(maskCols): 
        #phaseChange from slm to beamA
        beamA[i][j].changePhase(slm[i][j])

        #Bob phase changes
        phaseIntensity = [0.0 for p in range(numPhases)] #will save calculated values
        for p in range(numPhases): #goes through all 4 phase shifts and calculates R's  # + slm_B_mask[i][j]
            phaseIntensity[p] = 0.5 * (1 + math.cos(beamB[i][j].getPhase() + slm_B_mask[i][j] + phases[p] + beamA[i][j].getPhase())) #takes inital beamB, adds this phase change, and then A            
        
        #do reconstruction math with noise
        real = phaseIntensity[0] - phaseIntensity[2]
        imag = phaseIntensity[1] - phaseIntensity[3]
        real = real * random.uniform(ampNoiseLow, ampNoiseHigh)
        imag = imag * random.uniform(ampNoiseLow, ampNoiseHigh)
        c = phase(complex(real, imag))
        if c < 0:
            c += 2*math.pi #wraps negative phase shift to positive
        #reconstruction[maskRows-1-i][maskCols-1-j] = 2*math.pi - c #c is -theta_A, therefore 2pi-c = actual theta_A
        reconstruction[i][j] = 2*math.pi - c # takes away rotated aspect of reconstruction, for easier
        accuracy[i][j] = reconstruction[i][j] - slm[i][j]
        if accuracy[i][j] > math.pi:
            accuracy[i][j] -= 2*math.pi
        elif accuracy[i][j] < -1*math.pi:
            accuracy[i][j] += 2*math.pi

# summarizes accuracy
print("Mean Accuracy: "+ str(np.mean(np.abs(accuracy))) + "\t[" + str(np.min(accuracy)) + ", " + str(np.max(accuracy)) + "]")

# Begin Plotting of Results
axes=[]
fig = plt.figure()

axes.append(fig.add_subplot(1, 3, 1)) # reconstruction
plt.imshow( reconstruction, cmap = cm.gray)
subplot_title=("Reconstructed")
axes[-1].set_title(subplot_title)  

axes.append(fig.add_subplot(1, 3, 2)) #image
plt.imshow( slm, cmap = cm.gray)
subplot_title=("Reference Image")
axes[-1].set_title(subplot_title)  

axes.append(fig.add_subplot(1, 3, 3)) #image
plt.imshow( accuracy, cmap = cm.gray)
subplot_title=("Accuracy")
axes[-1].set_title(subplot_title) 

fig.tight_layout() 
plt.show()


        

