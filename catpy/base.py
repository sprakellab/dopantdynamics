"""This module contains functions for plotting dat based on particle position."""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import art3d

import catpy.IO.xml as xml

def Scatterplot2D(Positions,images,a,S,show=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    Slice = np.equal(Positions[:,2],round(a*S + 0.5,8))
    Selection = Positions[Slice,:]	
    ax.scatter(Selection[:,0], Selection[:,1],s=10000)
    ax.plot(images[:,0],images[:,1])
    ax.set_xlim(S*a,(S+1)*a+1)
    ax.set_ylim(S*a,(S+1)*a+1)	
    if show is True:	
        plt.show()

def Scatterplot3D(Input, data, images=None, show=False):
    """Plots particles in a 3D frame."""

    if type(data) == str:
        filename = '%s%s' % (Input['systempath'], data)
        f = open(filename)
        xml.ReadHeader(f)
        positions = xml.ReadCoordinates(f,Input['ALLparticles'])
        f.close()	
    else:
        positions = data
    # CHOOSING LIMITS	
    MinLim = np.zeros([3])	
    MaxLim = np.zeros([3])	
    Differences = np.zeros([3])	
    diff = 0

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(positions[:,0], positions[:,1], positions[:,2],s=10)
    if images is not None:
        ax.plot(images[:,0],images[:,1],images[:,2],color = 'r')
        for i in range(3):	
            MinLim[i] = min(images[:,i]) - 0.1
            MaxLim[i] = max(images[:,i]) + 0.1
            Differences[i] = MaxLim[i] - MinLim[i]
            if Differences[i] > diff:
                diff = Differences[i]
        ax.set_xlim(MinLim[0],MinLim[0]+diff)
        ax.set_ylim(MinLim[1],MinLim[1]+diff)
        ax.set_zlim(MinLim[2],MinLim[2]+diff)
    else:
        for i in range(3):	
            MinLim[i] = min(positions[:,i]) - 0.1
            MaxLim[i] = max(positions[:,i]) + 0.1
            Differences[i] = MaxLim[i] - MinLim[i]
            if Differences[i] > diff:
                diff = Differences[i]
        ax.set_xlim(MinLim[0],MinLim[0]+diff)
        ax.set_ylim(MinLim[1],MinLim[1]+diff)
        ax.set_zlim(MinLim[2],MinLim[2]+diff)		
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.axis('equal')
    if show is True:	
        plt.show()

def Scatterplot3DInterstitial(Input, data, show=False):
    """Plots particles in a 3D frame."""
    if type(data) == str:
        filename = '%s%s' % (Input['systempath'], data)
        f = open(filename)
        xml.ReadHeader(f)
        positions = xml.ReadCoordinates(f,Input['ALLparticles'])
        f.close()	
    else:
        positions = data
    # CHOOSING LIMITS	
    MinLim = np.zeros([3])	
    MaxLim = np.zeros([3])	
    Differences = np.zeros([3])	
    diff = 0
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(positions[Input['BASEparticles']:,0], positions[Input['BASEparticles']:,1], positions[Input['BASEparticles']:,2],s=10)
    for i in range(3):	
        MinLim[i] = min(positions[:,i]) - 0.1
        MaxLim[i] = max(positions[:,i]) + 0.1
        Differences[i] = MaxLim[i] - MinLim[i]
        if Differences[i] > diff:
            diff = Differences[i]
    ax.set_xlim(MinLim[0],MinLim[0]+diff)
    ax.set_ylim(MinLim[1],MinLim[1]+diff)
    ax.set_zlim(MinLim[2],MinLim[2]+diff)		
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.axis('equal')
    if show is True:	
        plt.show()

def PlotIPD(SaveName,Bin_Edges_Norm,TotalNormHist,show=False):
    """Plots the pair correlation function from calculated data.

    SaveName -- Name of the file to be saved without extension
    """
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.plot(Bin_Edges_Norm[1:(Bin_Edges_Norm.size-1)], TotalNormHist)
    plt.xlabel('r')
    plt.ylabel('G(r)')
    plt.xlim(0,10)
    pdfSaveName = SaveName + '.pdf'
    pngSaveName = SaveName + '.png'
    fig.savefig(pdfSaveName)
    fig.savefig(pngSaveName)
    if show is True:
        plt.show()
  
