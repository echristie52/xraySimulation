# Eric Christie
# BYU, 8/17/22
# Reconstructing images from the paper using original data
# Paper: Polarization entanglement-enabled quantum holography
# Data: https://researchdata.gla.ac.uk/1093/

import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from cmath import phase

'''
Indexing:   for data frames,    [col][row]
            for lists,          [row][col]
'''

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

numPhase = 4

reference = pd.read_excel("Fig1_reference.xlsx")
reconstruction = pd.read_excel("Fig1_reconstruction.xlsx")
phaseCorrelations = []
for p in range(numPhase):
    phaseCorrelations.append(pd.read_excel("Fig1_phase" + str(p) + ".xlsx"))

dimensions = phaseCorrelations[0].shape
#plot(phaseCorrelations, ["0", "pi/2", "pi", "3pi/2"], [False, False, False, False])





myReconstruction = []
accuracy = []
for i in range(dimensions[0]):
    myReconstruction.append([])
    for j in range(dimensions[1]):
        real = phaseCorrelations[0][j][i] - phaseCorrelations[2][j][i]
        imag = phaseCorrelations[1][j][i] - phaseCorrelations[3][j][i]

        c = -1*phase(complex(real, imag))
        #if abs(c) > 1:
        c %= 2*math.pi
        myReconstruction[i].append(c)


for i in range(dimensions[0]):
    accuracy.append([])
    for j in range(dimensions[1]):
        reconstruction[j][i] %= 2*math.pi
        accuracy[i].append(reconstruction[j][i] - myReconstruction[i][j])



plot([myReconstruction, reconstruction, accuracy], ["My Reconstruction", "Paper Reconstruction", "Accuracy"], [True, True, False])
#plot([myReconstruction, reconstruction], ["My Reconstruction", "Paper Reconstruction"], [False, False])

