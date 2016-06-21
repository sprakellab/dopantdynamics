"""Functions for plotting potential fields in crystal structures"""
from __future__ import division, absolute_import, print_function
from mayavi import mlab

import numpy as np
import catpy.IO.xml as xml
from colormaps import _viridis_data

class BaseCrystal:
    """crystal that defines the potential field"""

    def __init__(self, volumefraction, Nunitcells):
        self.vf = volumefraction
        self.n = Nunitcells # Number of unitcells
        self.m = self.n + 1
        self.n_particles = 2.0*self.m**3
        self.particle_volume = (4.0/3.0)*np.pi*((0.5)**3.0)
        self.Lbox = (self.n_particles*self.particle_volume/self.vf)**(1.0/3.0)
        self.a = (self.Lbox)/(self.m)
        self.BASEkappa = 1.8
        self.BASEepsilon = 713
        self.INTERkappa = 1.8
        self.INTERepsilon = 227
        self.coordinates = np.zeros([0,0])
        self.rcutoff = 10
        self.set_potential()

    def BCClatticesetup(self):
        """sets up periodic BCC crystal"""
        Offset = np.ones(3)*0.25*self.a
        positions = np.zeros([self.n_particles,3])
        count = 0
        count = xml.CubicLattice(positions, 
                                 count, 
                                 self.a, 
                                 Offset, 
                                 self.m, 
                                 self.m, 
                                 self.m)
        count = xml.CubicLattice(positions, 
                                 count, 
                                 self.a, 
                                 Offset, 
                                 self.m, 
                                 self.m, 
                                 self.m, 
                                 0.5, 
                                 0.5, 
                                 0.5)
        self.coordinates = positions

    def set_potential(self):
        self.epsilon = (self.BASEepsilon + self.INTERepsilon)/2.0
        self.kappa = (self.BASEkappa + self.INTERkappa)/2.0


def V(X, Y, Z, coordinates, exp=False, crystal_=None):
    size = len(X)
    potential = np.zeros([size, size, size])
    for i in range(size):
        for j in range(size):
            for k in range(size):
                voxel = np.array([X[i,j,k],Y[i,j,k],Z[i,j,k]])
                r = np.sqrt(np.sum(np.square(coordinates-voxel), axis=1))
                if exp:
                    cutoff = np.less(r,10)	
                else:
                    cutoff = np.less(r,crystal_.rcutoff)
                E = np.zeros(len(r))
                if exp:
                    E[:] = 447.0*(np.exp(-1.8*r)/r)
                else:
                    E[:] = crystal_.epsilon*(np.exp(-crystal_.kappa*r)/r)  
                potential[i,j,k] = np.sum(E*cutoff)
    return potential          

def calculate_isosurface(density,m):
    crystal = BaseCrystal(density, m)
    crystal.BCClatticesetup()
    average = np.average(crystal.coordinates[:,0])
    start = average - crystal.a
    end = average + crystal.a
    spacing = 100j
    X, Y, Z = np.mgrid[start:end:spacing, start:end:spacing, start:end:spacing]
    U = V(X, Y, Z, crystal.coordinates, crystal_=crystal)
    DeltaU = U - np.min(U)
    filename = './output/potentialfield.txt'
    data = np.hstack((np.reshape(X, (-1,1)),
                      np.reshape(Y, (-1,1)),
                      np.reshape(Z, (-1,1)),
                      np.reshape(DeltaU, (-1,1))
                    ))
    np.savetxt(filename, data, delimiter='\t')
    return crystal

def plot_isosurface(crystal):
    filename = './output/potentialfield.txt'
    data = np.genfromtxt(filename, delimiter='\t')
    size = np.round((len(data))**(1/3))
    X = np.reshape(data[:,0], (size,size,size))
    Y = np.reshape(data[:,1], (size,size,size))
    Z = np.reshape(data[:,2], (size,size,size))
    DeltaU = np.reshape(data[:,3], (size,size,size))
    average = np.average(crystal.coordinates[:,0])
    start = average - crystal.a
    end = average + crystal.a
    coords1 = np.array([[start, start, start]])
    coords2 = np.array([[end, end, end]])
    array1 = np.repeat(coords1,len(crystal.coordinates),axis=0)
    array2 = np.repeat(coords2,len(crystal.coordinates),axis=0)
    basefilter1 = np.greater(crystal.coordinates,array1)
    basefilter2 = np.less(crystal.coordinates,array2)
    basefilter = np.nonzero(np.all(basefilter1*basefilter2, axis=1))
    base = crystal.coordinates[basefilter]   

    mlab.figure(bgcolor=(1, 1, 1), fgcolor=(1, 1, 1), size=(2048,2048))
    dataset = mlab.contour3d(X, Y, Z, DeltaU, contours=[3.50],color=(1,0.25,0))
    scatter = mlab.points3d(base[:,0], 
                            base[:,1], 
                            base[:,2],
                            color=(0.255,0.647,0.88),
                            resolution=24, 
                            scale_factor=1.0, 
                            opacity=0.40)
    mlab.view(azimuth=17, elevation=90, distance=10, focalpoint=[average,average-0.2,average])
    mlab.draw()
    savename = './output/3Dpotential.png'
    mlab.savefig(savename, size=(2048,2048))
    mlab.show()

def calculate_isosurface_exp(coordinates, a):
    average = np.average(coordinates[:,0])
    start = average - a
    end = average + a
    spacing = 100j
    X, Y, Z = np.mgrid[start:end:spacing, start:end:spacing, start:end:spacing]
    U = V(X, Y, Z, coordinates, exp=True)
    DeltaU = U - np.min(U)
    print(np.max(DeltaU))
    filename = './output/exppotentialfield.txt'
    data = np.hstack((np.reshape(X, (-1,1)),
                      np.reshape(Y, (-1,1)),
                      np.reshape(Z, (-1,1)),
                      np.reshape(DeltaU, (-1,1))
                    ))
    np.savetxt(filename, data, delimiter='\t')
 
def plot_isosurface_exp(coordinates, a):
    filename = './output/exppotentialfield.txt'
    data = np.genfromtxt(filename, delimiter='\t')
    size = np.round((len(data))**(1/3))
    X = np.reshape(data[:,0], (size,size,size))
    Y = np.reshape(data[:,1], (size,size,size))
    Z = np.reshape(data[:,2], (size,size,size))
    DeltaU = np.reshape(data[:,3], (size,size,size))
    average = np.average(coordinates[:,0])
    start = average - a
    end = average + a
    coords1 = np.array([[start, start, start]])
    coords2 = np.array([[end, end, end]])
    array1 = np.repeat(coords1,len(coordinates),axis=0)
    array2 = np.repeat(coords2,len(coordinates),axis=0)
    basefilter1 = np.greater(coordinates,array1)
    basefilter2 = np.less(coordinates,array2)
    basefilter = np.nonzero(np.all(basefilter1*basefilter2, axis=1))
    base = coordinates[basefilter]   

    mlab.figure(bgcolor=(1, 1, 1), fgcolor=(1, 1, 1), size=(2048,2048))
    dataset = mlab.contour3d(X, Y, Z, DeltaU, contours=[42],color=(1,0.25,0))
    scatter = mlab.points3d(base[:,0], 
                            base[:,1], 
                            base[:,2],
                            color=(0.255,0.647,0.88),
                            resolution=24, 
                            scale_factor=1.0, 
                            opacity=0.40)
    mlab.view(azimuth=17, elevation=90, distance=11, focalpoint=[0,-0.2,0.2])
    mlab.draw()
    savename = './output/3Dpotentialsim.png'
    mlab.savefig(savename, size=(2048,2048))
    mlab.show()

