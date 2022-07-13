# Eric Christie
# BYU, 7/11/22
# Simulation of SLM and Image Reconstruction
# based on paper, Polarization entanglement-enabled quantum holography
# Reconstruct SLM Alice image given intensity correlation frames
# V1: Generate "frames" and reconstruct the image without noise
# Comments: maskRows/maskCols --> imageRows/imageCols

from matplotlib.colors import SymLogNorm
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




#Globals
global imageRows, imageCols #dimensions of image
global slm


#photon class to store magnitude and phase
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
        self.phase %= 2*math.pi # keeps in range [0, 2pi]


###### Helper Functions ######
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

#Pulls in image and updates globals
def uploadImage(filename):
    global imageRows, imageCols
    global slm

    testImage = Image.open(filename)
    testImage = ImageOps.grayscale(testImage)
    imageArray = asarray(testImage)
    slm = np.multiply(imageArray, 2*math.pi/255.0) # scales image to 0-2pi phase changes

    imageRows = imageArray.shape[0] 
    imageCols = imageArray.shape[1]




###### Beginning of Simulation ######
uploadImage("diagonal_test4.png")

# Constants
numPhases = 4
phases = np.multiply(math.pi, [0, 1/2, 1, 3/2])
totalFrames = 8#1000 # in paper, either 2.5e6 or 5e6
numFrames = int(totalFrames / numPhases)
detector_efficiency = 0.48




# Generate Photons and Intensity Frames

frames = [[[[0 for j in range(2*imageCols)] for i in range(imageRows)] for l in range(numFrames)] for p in range(numPhases)]






# Reconstruct from Intensity Frames

R_value = [0 for p in range(numPhases)] #rewritten for each pixel
reconstruction = [[0 for i in range(imageRows)] for j in range(imageCols)]

# Loops through i,j for +/- k
# for each phase shift, calculates intensity correlation from all frames
# Reconstructs with Eq 1
for i in range(imageRows):
    for j in range(imageCols):
        for p in range(numPhases):
            sum1 = 0
            sum2 = 0

            #Equation from Supplementary Information 1 (similar to Methods, Eq 4)
            for l in range(numFrames):
                # I_l(k) * I_l(-k)
                sum1 += frames[p][l][i][imageCols+j] * frames[p][l][imageRows-1-i][imageCols-1-j]
            
            for l in range(numFrames - 1):
                # I_l(k) * I_(l+1)(-k)
                sum2 += frames[p][l][i][imageCols+j] * frames[p][l+1][imageRows-1-i][imageCols-1-j]
            
            R_value[p] = ( sum1 / numFrames ) - ( sum2 / (numFrames-1) )
        
        real = R_value[0] - R_value[2] #reconstruction values
        imag = R_value[1] - R_value[3]

        reconstruction[i][j] = phase(complex(real, imag)) #extracts phase
        # reconstruction[i][j] -= psi0_correction[i][j] # apply static phase correction mask

# plot Reconstruction
plot([reconstruction, slm], ["Reconstruction", "Reference"], [True, True])
