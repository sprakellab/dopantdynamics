"""This module contains functions for analyzing data based on 
particle positions. Here specifically Pair correlationfunctions.
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy.io import savemat

import catpy.IO.xml as xml

def SinglePairCorrelationFunction(Input, filename, binsize=0.05):
    """
    This script calculates the pair correlation functions from simulation
    data.

    Keyword arguments:
    Filename -- Filename from script.
    binsize  -- The size of a bin in the histogram (default 0.05). 
    calculated (default 10.0). 
    """
    numberofbins = np.ceil(Input['max_range_IPD']/binsize)	
    bin_edges, normalizedshellvolume = InitializationEdgesVolumes(Input, binsize, numberofbins)
    singlenormhist = SingleHist(Input, filename, bin_edges, normalizedshellvolume)
    bin_edges_norm = RNormalization(singlenormhist, bin_edges)
    return singlenormhist, bin_edges_norm

def MultiplePairCorrelationFunction(Input, binsize=0.05):
    """Calculates the pair correlation functions from simulation data. 

    Input is a whole folder of simulation data called 
    FILENAME.xxxxxxxxx.xml. Where xxxxxxxx stands for the stepNumber 
    of the data.
    
    Keyword arguments:
    Input -- Dict with all settings for simulation and analysis.
    binsize -- The size of a bin in the histogram (default 0.05).
    maxrange -- Maximum range in which the Pair Correlation Function 
    is calculated (default 10.0). 

    Output:	
    TotalNormHist -- Normalized intensities of the pair correlation 
    histogram average of all histograms 	
    Bin_Edges_Norm -- Normalized values for bin_edges 
    """
    numberofIPD = ((Input['last_step_IPD']-Input['first_step_IPD']) / Input['XMLsaveperiod'])
    numberofbins = np.ceil(Input['max_range_IPD']/binsize)
    bin_edges, normalizedshellvolume = InitializationEdgesVolumes(Input, binsize, numberofbins)
    singlenormhist = np.zeros([numberofbins-1, numberofIPD])
    totalnormhist = np.zeros([numberofbins-1])
    rowcounter = 0
    for i in range(Input['first_step_IPD'],Input['last_step_IPD'],
                   Input['XMLsaveperiod']):
        filename = '%s%s.%010d.%s' % (Input['folder_snapshot_data'],
                                      Input['snapshot_data_name'],
                                      i,
                                      Input['snapshot_data_extension'])	
        singlenormhist[:,rowcounter] = SingleHist(Input, filename, bin_edges, normalizedshellvolume)
        totalnormhist += singlenormhist[:,rowcounter]
        rowcounter += rowcounter
    totalnormhist = totalnormhist / numberofIPD
    bin_edges_norm = RNormalization(totalnormhist, bin_edges)
    return totalnormhist,bin_edges_norm

# SUPPORTING FUNCTIONS

def SetSmallBox(Input, lx, ly, lz):
    """
    Determining Dimensions of the inner box.

    Keyword arguments:
    lx -- size x dimension complete box.
    ly -- size y dimension complete box.
    lz -- size z dimension complete box.
    """
    l = np.array([lx,ly,lz])
    s = (l/2.0) - Input['max_range_IPD']
    return s[0], s[1], s[2]

def IdentifySmallBoxParticles(positions,sx,sy,sz):
    """Analyzing which particles are inside or outside the box.

    Keyword arguments:	
    positions -- Position list of all particle positions.
    sx -- size x dimension inner box.
    sy -- size y dimension inner box. 
    sz -- size z dimension inner box. 
    
    Output: 
    refparticles -- Position list of all particles within the smallbox
    """
    ref_particle_indx_x = np.abs(positions[:,0]) < sx
    ref_particle_indx_y = np.abs(positions[:,1]) < sy
    ref_particle_indx_z = np.abs(positions[:,2]) < sz
    ref_particle_indx_tot = np.logical_and(np.logical_and(ref_particle_indx_x, 
                                                          ref_particle_indx_y), 
                                                          ref_particle_indx_z)
    refparticles = positions[ref_particle_indx_tot, :]
    return refparticles

def IPDCalculation(Input, refparticles, positions):
    """Calculating the Inter Particle Distances.

    Keyword arguments:
    refparticles -- Position List of all particles in the smallbox 
    positions -- Position List of all the particles in the whole box
    MaxRange -- Maximum range in which the Pair Correlation Function is
    calculated 
    
    Output:	
    distances -- A list of all the inter particle distances between the
    particles in the smallbox and the wholebox 
    """	
    interparticledistance = cdist(refparticles, positions, 'euclidean')
    distances_indx = interparticledistance < Input['max_range_IPD']
    distances = interparticledistance[distances_indx]
    return distances

def NormalizedShellVolume(bin_edges, binsize):
    """Normalization for the Shellvolume.
    
    Keyword arguments:	
    bin_edges -- The edges of the histogram bins 
    binsize -- The size of a bin in the histogram 
    
    Output:
    normshellvolumes -- A list with the normalized volumes of shells
    used for normalizing the Pair Correlation Function
    """
    shellvolumes = (((4.0*np.pi/3.0)*(bin_edges + binsize)**3.0) - 
                    ((4.0*np.pi/3.0)*(bin_edges)**3.0))
    shellvolumes = shellvolumes[:(shellvolumes.size-1)]
    normshellvolumes = shellvolumes / np.sum(shellvolumes)
    return normshellvolumes

def RandomIPDDistribution(distances,normalizedshellvolume):
    """Calculating the distribution of IPD's for a random distribution.

    Keyword arguments:	
    distances -- A list of all the inter particle distances between
    the particles in the smallbox and the wholebox.
    normalizedshellvolume -- The shell volumes normalized to the
    complete volume. 
    
    Output: 
    randomipd -- The pair correlation function of randomly divided 
    particles.
    """		
    randomipd = normalizedshellvolume * distances.size
    return randomipd	

def NormalizedIPD(distances,hist,normalizedshellvolume):
    """Normalizing the Pair Correlation plot.

    Keyword arguments: 	
    distances -- A list of all the inter particle distances between 
    the particles in the smallbox and the wholebox 
    hist -- The histogram produced from the distances list 
    normalizedshellvolume -- The shell volumes normalized to the
    complete volume 
    
    Output:	
    hist -- Histogram normalized for a random distribution of the same
    amount of particles
    """		
    RandomIPD = RandomIPDDistribution(distances,normalizedshellvolume)
    hist = hist[1:] / RandomIPD[1:(RandomIPD.size)] # Delete autocorrelations
    return hist

def InitializationEdgesVolumes(Input, binsize, numberofbins):
    """Setting the bin edges and normalizing shellvolumes.

    Keyword arguments:
    Bindize -- The size of a bin in the histogram 
    numberofbins -- number of bins

    Output:	
    bin_edges -- The r values for the edges of the histogram bins 
    normshellvolume -- The shell volumes normalized to the complete volume 
    """
    # Setting bin edges
    Dummy, bin_edges = np.histogram(0, bins=numberofbins, range=(0.0,Input['max_range_IPD']))
    # calculating shell volumes for normalization
    normshellvolume = NormalizedShellVolume(bin_edges, binsize) 
    return bin_edges, normshellvolume

def SingleHist(Input, Filename, Bin_Edges, normalizedshellvolume):
    """Calculating a single histogram and normalizing over G.
    
    Keyword arguments:
    Filename -- Filename of the xml positions file 
    Bin_Edges -- The r values for the edges of the histogram bins
    normalizedshellvolume -- The shell volumes normalized to the 
    complete volume 
    
    Output:	
    singlenormhist -- Histogram intesities produced from a single 
    xml file 
    Positions -- A list with the positions of all particles 
    """
    F = open(Filename, 'r')
    TimeStep,Dimensions,NParticles,Lx,Ly,Lz = xml.ReadHeader(F)
    Positions = xml.ReadCoordinates(F,Input['BASEparticles'])
    F.close()
    Sx,Sy,Sz = SetSmallBox(Input, Lx, Ly, Lz)
    RefParticles = IdentifySmallBoxParticles(Positions, Sx, Sy, Sz)
    Distances = IPDCalculation(Input, RefParticles, Positions)
    Hist, Dump = np.histogram(Distances, bins=Bin_Edges)
    singlenormhist = NormalizedIPD(Distances, Hist, normalizedshellvolume)
    return singlenormhist

def RNormalization(NormHist, Bin_Edges):
    """Normalizes the histogram over r, based on the first largest peak.

    Keyword arguments:
    NormHist -- Histogram normalized for the volume 
    bin_edges -- The r values for the edges of the histogram bins 

    Output:
    bin_edges_norm -- Normalized values for bin_edges 
    """
    bin_edges_norm = Bin_Edges[:]
    return bin_edges_norm

def SaveIPD(Input,Bin_Edges_Norm,NormHist):
    """Exports the pair correlation function as a .mat and or .txt file.

    Keyword_arguments:
    Bin_Edges_Norm -- Normalized values for bin_edges 
    NormHist -- Normalized intensities of the pair correlation histogram 
    """
    if Input['savetxtIPD'] is True:
        txtSaveName = Input['savenameIPD'] + '.txt'
        s = open(txtSaveName, "w")
        for i in range(NormHist.size):
            s.write(str(Bin_Edges_Norm[i+1]) + '\t' + str(NormHist[i])+ '\n')
        s.close()
    if Input['savematIPD'] is True:	
        matSaveName = Input['savenameIPD'] + 'mat'
        savemat(matSaveName, {'pos':positions})

def Sxy_operator(Input, cylR, sigma):
    filename = '%sWigner.0000000030.xml' % (Input['folder_snapshot_data'])
    f = open(filename, 'r')
    xml.ReadHeader(f)
    data = xml.ReadCoordinates(f, Input['BASEparticles'])
    f.close()

    #### tmp is a random frame (x,y,z) 
    #### Get rotation matrix:
    #rotMat,fullMat = rotation_matrix(0.25*np.pi, [0, 1, 0], point=None)
    #data = np.dot(data, rotMat) # This is where the magic happens

    grid, S = Sxy(data, cylR, sigma)
    x = np.reshape(grid[:,:,0],(-1,1))    
    y = np.reshape(grid[:,:,1],(-1,1))
    sofq = np.reshape(S,(-1,1))
    output = np.concatenate((x,y,sofq),axis=1)
    savename = '%s%ssofqsnapshot.txt' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, output, delimiter='\t')

def Sxy(positions, cylR, sigma):
    """Calculates S(q) in the x,y plane.

    Keyword_arguments:
    positions -- Positions of particles (x, y, z). 
    cylR -- Amount of shells per period.
    sigma -- Width of the qrange. 'sigma=1' gives one period,
             'sigma = 2' gives half a period etc.

    output: 
    grid -- Numpy array of square matrices with qx and qy values.
    S -- Square matrix corresponding to grid with S values.
    """
    qrange = np.arange(-4.0*np.pi/sigma,4.0*np.pi/sigma, 4.0*np.pi/cylR)
    qlen = len(qrange)
    beadsinframe = len(positions)
    xgrid, ygrid = np.meshgrid(qrange, qrange)
    grid = np.concatenate((np.reshape(xgrid,(qlen,qlen,1)),
                           np.reshape(ygrid, (qlen,qlen,1))),
                           axis=2)
    flatarray = np.zeros([beadsinframe, 2])
    flatarray[:,:] = positions[:,:2]
    centerofmass = np.mean(flatarray, axis=0)
    differencevectors = flatarray - centerofmass
    
    #compute the particle by particle contribution to the S(q)   
    gridperparticle = np.repeat(np.reshape(grid,(qlen,qlen,2,1)), 
                                beadsinframe, 
                                axis=3)
    exponent = np.sum(gridperparticle*differencevectors.T, axis=2)
    sq_raw = np.sum(np.exp(1j*(exponent)), axis=2)
    
    #normalize by number of beads and remove self contribution 
    S = np.abs(np.square(sq_raw)) / beadsinframe - 1.0
    
    # if q=0 is computed then minor corrections need 
    # to be applied to avoid spurious delta function.
    q0 = np.nonzero(np.equal(qrange, 0.0))
    S[q0,q0] = S[q0,q0] - (beadsinframe-1)
    
    # take out delta function at origin from fourier 
    # transforming g(r) rather than h(r)
    midpoint = int(qlen/2.0)
    S[midpoint-1:midpoint+1,midpoint-1:midpoint+1] = 0
    return grid, S

### ROTATIE CODE:
def rotation_matrix(angle, direction, point=None):
    """Return matrix to rotate about axis defined by point and direction.
       Goed gejat is beter dan slecht verzonnen, Ruben 2016.
    """
    def unit_vector(data, axis=None, out=None):
        """Return ndarray normalized by length, i.e. Euclidean norm, along axis.
        """
        if out is None:
            data = np.array(data, dtype=np.float64, copy=True)
            if data.ndim == 1:
                data /= np.sqrt(np.dot(data, data))
                return data
        else:
            if out is not data:
                out[:] = np.array(data, copy=False)
            data = out
        length = np.atleast_1d(np.sum(data*data, axis))
        np.sqrt(length, length)
        if axis is not None:
            length = np.expand_dims(length, axis)
        data /= length
        if out is None:
            return data

    sina = np.sin(angle)
    cosa = np.cos(angle)
    direction = unit_vector(direction[:3])
    # rotation matrix around unit vector
    R = np.diag([cosa, cosa, cosa])
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array([[ 0.0,         -direction[2],  direction[1]],
                      [ direction[2], 0.0,          -direction[0]],
                      [-direction[1], direction[0],  0.0]])
    M = np.identity(4)
    M[:3, :3] = R
    if point is not None:
        # rotation not around origin
        point = np.array(point[:3], dtype=ny.float64, copy=False)
        M[:3, 3] = point - np.dot(R, point)
    return (R,M)
