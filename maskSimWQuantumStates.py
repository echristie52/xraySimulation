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
import enum
from statistics import mean, stdev


testImage = Image.open("diagonal_test.png")
testImage = ImageOps.grayscale(testImage)
imageArray = asarray(testImage)

#numRows must be even, so I can split the beam in half
global maskRows, maskCols
maskRows = imageArray.shape[0]
maskCols = imageArray.shape[1]

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def getTheoreticalR(phi):
    return 0.5 * (A**2 + B**2 + A*B*math.cos(phi))

def getSinglePhotonR(phi):
    return 0.5 * (Ao**2 + Bo**2 + 2*Ao*Bo*math.cos(phi))
    #return getTheoreticalR(phi)**0.5

#Bob and Alice are references to the Nature paper
#Bob is a constant phase shift SLM
#Alice is the image SLM
class Frame:
    def __init__(self, phi0, w, h):
        self.phi0 = phi0
        self.bob = [[0 for i in range(w)] for j in range(h)]
        self.alice = [[0 for i in range(w)] for j in range(h)]
        
    def getBobFlip(self):
        return np.flipud(np.fliplr(self.bob))

class Polarizations(enum.Enum):
    H = 0
    V = 1

class Photon:
        
    def __init__(self, polarization):
        self.polarization = polarization
        self.intensity = 1
        
        if self.polarization == Polarizations.H:
           self.phase = 0
        else:
           self.phase = pi/2
           
    def setP(self, p: Polarizations):
        self.polarization = p
        if p == Polarizations.H:
            self.phase = 0
        else:
            self.phase = pi/2
          

class PhotonPair:
    global maskRows, maskCols
    def __init__(self, x, y):
        self.x = x
        self.y = y
        ## vv not convinced these are right - dont think it matters tho
        self.minusX = maskRows - x - 1
        self.minusY = maskCols - y - 1
        
        HV_rand = random.randint(0,1);
        if HV_rand == 0:
            self.pairP = Polarizations.H
        elif HV_rand == 1:
            self.pairP = Polarizations.V
            
        #Photon 1 goes through BOB
        self.P1 = Photon(self.pairP)
        #Photon 2 goes through ALICE
        self.P2 = Photon(self.pairP)
        
def doSlm(pp: PhotonPair, bobSLM, aliceSLM):
    pp.P1.phase = pp.P1.phase + bobSLM[pp.x][pp.y]
    pp.P2.phase = pp.P2.phase + aliceSLM[pp.x][pp.y]
        
def polarize(pp: PhotonPair):
    phi = pp.P1.phase + pp.P2.phase;
    rd = getTheoreticalR(phi)
    population = [0,1]
    weights = [rd, 1-rd]
    clickOrNot = choices(population,weights)
    if clickOrNot == 0:
        pp.P1.intensity = 0
        pp.P2.intensity = 0
        
def polarizeQ(pp: PhotonPair):
    rd1 = math.cos(pp.P1.phase - polarizer_angle_h)**2
    population = [0,1]
    weights = [1-rd1, rd1]
    p1Pass = choices(population,weights)[0]
    pp.P1.intensity = p1Pass
    rd2 = math.cos(pp.P2.phase - polarizer_angle_h)**2
    weights = [1-rd2, rd2]
    p2Pass = choices(population,weights)[0]
    pp.P2.intensity = p2Pass    
    
def fakePolarizeQ(pp: PhotonPair):

    p1_10 = math.cos(-pp.P1.phase + polarizer_angle_h)**2
    p2_10 = math.cos(-pp.P2.phase + polarizer_angle_h)**2
    
    pp.P1.intensity = p1_10
    pp.P2.intensity = p2_10
    



#Setup phase mask properties
polarizer_angle_h = pi/4
polarizer_angle_v = pi/4
CS = math.cos(polarizer_angle_h)*math.cos(polarizer_angle_v) + math.sin(polarizer_angle_h)*math.sin(polarizer_angle_v)
A = (math.cos(polarizer_angle_v)**2)/(CS**2)
B = (math.sin(polarizer_angle_v)**2)/(CS**2)
Ao = (math.cos(polarizer_angle_v))/(CS)
Bo = (math.sin(polarizer_angle_v))/(CS) 

#For generating a gausian probability over the frame
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


phase_mask =  [[1 for i in range(maskCols)] for j in range(maskRows)]

pm_test = np.divide(pi,1+np.arange(maskRows*maskCols).reshape((maskRows,maskCols))) 
pm_test2 = np.multiply(1/pi, np.random.randint(-0, high=8, size=(maskRows,maskCols)))
pm_test3 = np.multiply(imageArray, pi/255)
pm_test4 = np.array([[pi/4, pi/4],[pi/4, pi/4]])
phi_mask = pm_test3

#maskRows = 2;
#maskCols = 2;

aliceSLM = phi_mask
    

numPhi0 = 4
phi0 = np.multiply(math.pi,[0, 1/2, 1, 3/2])
#phi0 = np.multiply(math.pi,[1])

row_data = [[[0 for i in range(maskCols)] for j in range(maskRows)] for k in range(numPhi0)]
row_countsA = [[[0 for i in range(maskCols)] for j in range(maskRows)] for k in range(numPhi0)]
row_countsB = [[[0 for i in range(maskCols)] for j in range(maskRows)] for k in range(numPhi0)]
row_counts_100 = [[[0 for i in range(maskCols)] for j in range(maskRows)] for k in range(numPhi0)]

beamDatai = []
beamDatak = []

#Run a bunch photon pais per reference phase shift
numFrames = 100
frames = [[Frame(0, maskCols, maskRows) for i in range(numFrames)] for k in range(numPhi0)]

for p0 in range(len(phi0)):
    
    bobSLM = np.multiply(np.ones((maskRows, maskCols)), phi0[p0])
    

    cFrame = 0    

    while cFrame < numFrames:

        photonsPerFrame = 10
        cPhoton = 0
        #get current frame from array
        frame = Frame(phi0[p0], maskCols, maskRows)
        
        while cPhoton < photonsPerFrame:

            """
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
                cFrame = cFrame + 1
                continue
            if randI < 0 or randK < 0:
                cFrame = cFrame + 1
                continue
           
            #lets control the rand I rand K pattern shall we?
            randI = 0
            randK = 0
           
            #calculate the phase shift based on the coords and mask info
            #this allows me to calculate R, the probability of getting a click
            phasea = phase(A_post_slm[randI][randK]) 
            phaseb = phase(B_post_slm[randI][randK])
            phi = phaseb + phasea #+phase distortion?
            rd = getTheoreticalR(phi)
            if(row_data[p0][randI][randK] == 0):
                row_data[p0][randI][randK] = rd
            
            #detector efficiency is about 48%
            deff = 0.99#.48
            rd_x_deff = rd * deff
            population = [0,1]
            weights = [rd, 1-rd]
            clickOrNot = choices(population,weights)
            
             
            #create false coincidence x% of the time
            false_coin_rate = 0.0005 #1 %
            population = [0,1]
            weights = [1-false_coin_rate, false_coin_rate]
            fc_or_not = choices(population,weights)
            
            
            #create an entangled photon pair
            pp = PhotonPair(randI, randK)
            
            
            if pp.P1.polarization == Polarizations.V:
                #pass the photons through the SLM
                doSlm(pp, slm, slm_mask)
            
            #pass the photon pair through the polarizer
            fakePolarizeQ(pp)
            
            if fc_or_not[0] == 1:
                row_counts[p0][randI][randK] = row_counts[p0][randI][randK] + 1
           
            
            row_counts[p0][randI][randK] = row_counts[p0][randI][randK] + 1
            row_counts_100[p0][randI][randK] = row_counts_100[p0][randI][randK] + 1
            
            #capture data in the frame
            frame.bob[pp.minusX][pp.minusY] += pp.P1.intensity
            frame.alice[pp.x][pp.y] += pp.P2.intensity
            
            cPhoton += 1
            
            ##above code does a uniform distribution
            ##below code for loops through all the cells
            """
            
            for i in range(maskRows):
                for k in range(maskCols):
                    phasea = bobSLM[i][k] 
                    phaseb = aliceSLM[i][k]
                    rd = 0.25 * (1+ math.cos(phaseb + phasea))
                    row_data[p0][i][k] = rd
                

                    
                    #create an entangled photon pair
                    pp = PhotonPair(i, k)
                    #pp.P1.setP(Polarizations.H)
                    #pp.P2.setP(Polarizations.H)
                    if pp.P1.polarization == Polarizations.H:
                        #pass the photons through the SLM
                        doSlm(pp, bobSLM, aliceSLM)
                    
                    #pass the photon pair through the polarizer
                    fakePolarizeQ(pp)
                    
                    #capture data in the frame
                    #print(pp.P1.intensity)
                    #print(pp.P2.intensity)
                    frame.bob[pp.minusX][pp.minusY] += pp.P1.intensity
                    frame.alice[pp.x][pp.y] += pp.P2.intensity
                    #print(frame.bob[pp.minusX][pp.minusY])
                    #print(frame.alice[pp.x][pp.y])
                    
                    row_countsB[p0][i][k] += pp.P1.intensity #bob counts
                    row_countsA[p0][i][k] += pp.P2.intensity #alice counts
                    row_counts_100[p0][i][k] += + 1
                    
            cPhoton += 1
            
            
        frames[p0][cFrame] = frame
        cFrame = cFrame + 1
                

#print(row_counts[0][randI][randK])
#print(row_counts_100[0][randI][randK])

#plt.hist(beamDatai)
#plt.figure()
#plt.hist(beamDatak)

#extract R (The probability that a photon pair will click given a phase shift)
#from how many photons were aimed at a cell / how many clicke
#in reality we can't know this but we may be able to estimate it
row_data = np.array(row_data)
extractedR = np.divide(row_countsA, row_counts_100)

probDif = np.subtract(extractedR, row_data)

##NEXT STEP - use the frame data to create the R function as shown in the paper


correlation_Q = np.zeros((numPhi0,maskRows,maskCols))     
for p0 in range(numPhi0):
    for i in range(numFrames):
        if i == numFrames - 1:
            continue
        frame = frames[p0][i]
        nextFrame = frames[p0][i+1]
        mult1 = np.multiply(frame.alice, frame.getBobFlip()) 
        mult2 = np.multiply(nextFrame.getBobFlip(), frame.alice)
        accum = np.subtract(mult1,mult2)
        correlation_Q[p0] = correlation_Q[p0] + accum


emccdNoiseLow = 1
emccdNoiseHigh = 1



correlation_Q = np.divide(correlation_Q, 4*numFrames)
    

extractedR = np.divide(row_countsA, row_counts_100)
correlation_Q = extractedR
correlation = [[0 for i in range(maskCols)] for j in range(maskRows)]
correlation_recon = [[0 for i in range(maskCols)] for j in range(maskRows)]


for i in range(maskRows):
    for k in range(maskCols):
        real = row_data[0][i][k] - row_data[2][i][k]
        imag = row_data[1][i][k] - row_data[3][i][k]
        
        #simulate amplification noise with a random multiplier
        ampNoiseLow = 1
        ampNoiseHigh = 1
        
        real = real * random.uniform(ampNoiseLow, ampNoiseHigh)
        imag = imag * random.uniform(ampNoiseLow, ampNoiseHigh)
        c = complex(real, imag)
        correlation[i][k] = -phase(c)

        real_recon = correlation_Q[0][i][k] - correlation_Q[2][i][k]
        imag_recon = correlation_Q[1][i][k] - correlation_Q[3][i][k]
        
        #divide by 100 then mulitply by not quite 100 to simulate amplification noise
        real_recon = real_recon * random.uniform(ampNoiseLow, ampNoiseHigh)
        imag_recon = imag_recon * random.uniform(ampNoiseLow, ampNoiseHigh)
        c_recon = complex(real_recon, imag_recon)
        correlation_recon[i][k] = -phase(c_recon)

                
correlation_recon = np.absolute(correlation_recon)
        



correlation = np.array(correlation)
correlation_recon = np.array(correlation_recon)
row_countsA = np.array(row_countsA)
row_countsB = np.array(row_countsB)



cr_norm = NormalizeData(correlation_recon)
cr_norm_idea = NormalizeData(correlation)
im_norm = phi_mask
    
axes=[]
fig = plt.figure()

#try to calc SNR of reconstruction
ap = np.mean(correlation, dtype=np.float64)
noise = np.subtract(correlation, ap)
std_noise = np.std(noise)
SNR = 20*math.log10( ap / std_noise )

#print(cr_norm)
#print(cr_norm_idea)
#print(SNR)

axes.append(fig.add_subplot(2, 2, 2))
plt.imshow(correlation, vmin=0, vmax = math.pi, cmap = cm.gray)
subplot_title=("Ideal Reconstructed")
axes[-1].set_title(subplot_title)  


axes.append(fig.add_subplot(2, 2, 1))
plt.imshow( aliceSLM, vmin=0, vmax = math.pi, cmap = cm.gray)
subplot_title=("Reference Image")
axes[-1].set_title(subplot_title)  

axes.append(fig.add_subplot(2, 2, 3))
plt.imshow( row_countsA[0], cmap = cm.gray)
subplot_title=("Clicks per cell")
axes[-1].set_title(subplot_title)  

axes.append(fig.add_subplot(2, 2, 4))
plt.imshow( correlation_recon, cmap = cm.gray)
subplot_title=("Reconstructed")
axes[-1].set_title(subplot_title)  



fig.tight_layout() 
plt.show()



#plt.plot(phaseOfSlm)
#factor = np.divide(correlation, phaseOfSlm)
    



