"""Module for calculating variation in the Activation Energy of a crystal."""

import numpy as np
import matplotlib.pyplot as plt

import catpy.IO.xml as xml
import catpy.basisfunctions as bf

def calculate_variation(Input, resolution):
    """Calculates the variation in the activation Energy"""
    filelist = bf.lookupfiles(path=Input['folder_snapshot_data'], 
                              filechr='Wigner')
    filenumbers = np.zeros((len(filelist)))
    i = 0
    for filename in filelist:
        filepath = '%s%s' % (Input['folder_snapshot_data'], filename)
        f = open(filepath, 'r')
        xml.ReadHeader(f)
        positions = xml.ReadCoordinates(f, Input['BASEparticles'])
        f.close()
        coordinates, potential = potential_slice(Input, positions, resolution)
        if i == 0:
            all_potential = potential.reshape((1,len(potential),len(potential)))
        else:
            all_potential = np.append(all_potential,
                                      potential.reshape((1,len(potential),len(potential))),
                                      axis=0)
        
        filenumber = int(filename.strip('Wigner.').strip('.xml'))
        filenumbers[i] = filenumber
        if i == 0:
             save_potential_field(Input, coordinates, potential, 'testje')
        i += 1
    sorted_filenumbers = np.argsort(filenumbers,axis=0)
    sorted_potentials = all_potential[sorted_filenumbers]
    average_potential = np.average(sorted_potentials,axis=0)
    variance_potential = np.var(sorted_potentials,axis=0)
    coordx, coordy = np.nonzero(np.less(variance_potential,100)) 
    print len(coordx)
    if len(coordx) > 0:
        time = np.sort(sorted_filenumbers)*Input['dt']*Input['BASEtau']
    #eplt.plotsliceheatmap(Input, coordinates, potential, positions, resolution, plotpart=False)
        eplt.plot_energy_track(time, sorted_potentials[:,coordx[0],coordy[0]] - np.min(sorted_potentials[:,coordx[0],coordy[0]]))
    #eplt.plotsliceheatmap_variance(Input, coordinates, variance_potential, positions, resolution, Max=100, plotpart=False)

def calculate_potential_slice(Input, resolution, experimental=False):
    filepath = '%sWigner.0099000000.xml' % (Input['folder_snapshot_data'])
    if experimental is True:
        filepath = '%sWigner.0000000025.xml' % (Input['folder_snapshot_data'])
    f = open(filepath, 'r')
    xml.ReadHeader(f)
    positions = xml.ReadCoordinates(f, Input['BASEparticles'])
    f.close()
    coordinates, potential = potential_slice(Input, positions, resolution)
    print coordinates
    save_potential_field(Input, coordinates, potential)

def potential_slice(Input, positions, resolution):
    """Calculates the potential energy field in a slice of crystal"""
    unitcell = (Input['n'] / 2) + 1
    start = -0.5*Input['BOXsize'] + 0.25*Input['a'] + unitcell*Input['a']
    end = -0.5*Input['BOXsize'] + 0.25*Input['a'] + (unitcell+1)*Input['a'] 
    z = start
    coordinates = np.mgrid[start:end:(1j*resolution),start:end:(1j*resolution),z:(z+1)].reshape(3,resolution,resolution)
    potential = calculate_potential(Input, positions, coordinates)
    return coordinates, potential

def calculate_deviations(Input):
    """Calculates the deviations from the average particle site"""
    filelist = bf.lookupfiles(path=Input['folder_snapshot_data'], 
                              filechr='Wigner')
    i = 0
    for filename in filelist:
        if int(filename.strip('Wigner').strip('.xml')) > 10000000:
            filepath = '%s%s' % (Input['folder_snapshot_data'], filename)
            f = open(filepath, 'r')
            xml.ReadHeader(f)
            positions = xml.ReadCoordinates(f, Input['ALLparticles'])
            f.close()
            if i == 0:
                all_positions = np.reshape(positions, (1, len(positions), 3))
            else:
                positions = unfold_boundary_crossings(Input, all_positions, positions, i)
                positions = np.reshape(positions, (1, len(positions), 3))
                all_positions = np.append(all_positions, positions, axis=0)
            i += 1
    average_positions = np.average(all_positions, axis=0)
    deviations = np.std(np.sqrt(np.sum(np.square(all_positions - average_positions), axis=2)),axis=0)
    filepath = '%sWigner.0099000000.xml' % (Input['folder_snapshot_data'])
    f = open(filepath, 'r')
    xml.ReadHeader(f)
    positions = xml.ReadCoordinates(f, Input['ALLparticles'])
    f.close()
    positions = unfold_boundary_crossings(Input, all_positions, positions, i)
    snapdeviations = np.sqrt(np.sum(np.square(average_positions - positions), axis=1))
    posanddev = np.zeros((Input['ALLparticles'],8))
    posanddev[:,:3] = average_positions
    posanddev[:,3] = deviations / Input['a']
    posanddev[:,4:7] = positions
    posanddev[:,7] = snapdeviations / Input['a']
    head = '%-18s\t%-18s\t%-18s\t%-18s\t%-18s\t%-18s\t%-18s\t%-18s' % ('x[sigma]', 'y[sigma]', 'z[sigma]', 'lindemann','snapx[sigma]', 'snapy[sigma]', 'snapz[sigma]', 'snapdev[sigma]')
    savename = '%s%s.deviations' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, posanddev, fmt='%18.11e', header=head, delimiter='\t')
    return positions, deviations

def unfold_boundary_crossings(Input, all_positions, positions, i):
    test1 = all_positions[i-1] - positions
    test2 = -1.0*(all_positions[i-1] - positions)
    test_i1 = np.nonzero(np.greater(test1, 0.80*Input['BOXsize']))
    test_i2 = np.nonzero(np.greater(test2, 0.80*Input['BOXsize']))
    positions[test_i1[0],test_i1[1]] = (all_positions[i-1,test_i1[0],test_i1[1]] + 
                                        (Input['BOXsize'] - test1[test_i1[0],test_i1[1]]))
    positions[test_i2[0],test_i2[1]] = (all_positions[i-1,test_i2[0],test_i2[1]] - 
                                        (Input['BOXsize'] - test2[test_i2[0],test_i2[1]]))
    return positions

def calculate_potential(Input, positions, coordinates):
    dlen = np.tile(coordinates.T, (len(positions),1,1,1)).T - positions.reshape(-1,1,1,3).T
    r = np.sqrt(np.sum(np.square(dlen),axis=0))
    cutoff = np.less(r, Input['rcut'])
    Epairs = Eyp(r, Input['IBepsilon'], Input['IBkappa'])
    potential = np.sum(Epairs*cutoff, axis=2).T # correction for transposing in the above steps
    return potential

def Eyp(r, epsilon, kappa):
    """Definition for the Yukawa potential.
    
    Keyword arguments:
    r -- The inter particle distance.
    epsilon -- A yukawa paremeter with the dimensions energy*distance.
    kappa -- The inverse screening length.
    
    Output:
    epsilon*(np.exp(-kappa*r)/r) -- The yukawa potential energy	
    """
    return  epsilon*(np.exp(-kappa*r)/r)

def save_potential_field(Input, coordinates, potential):
    x = np.reshape(coordinates[0],(-1,1))    
    y = np.reshape(coordinates[1],(-1,1)) 
    z = np.reshape(coordinates[2],(-1,1))
    p = np.reshape(potential,(-1,1)) 
    data = np.concatenate((x,y,z,p),axis=1)
    text = '%-18s\t%-18s\t%-18s\t%-18s' % ('x[sigma]', 'y[sigma]', 'z[sigma]', 'potential[kT]')
    savename = '%s%spotentialfieldsnapshot.txt' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, data, delimiter='\t', header=text)

def read_potential_field(Input):
    filepath = '%s%spotentialfieldsnapshot.txt' % (Input['folder_processed_data'], Input['system'])
    data = np.genfromtxt(filepath, delimiter='\t', skip_header=1)
    side_length = np.sqrt(len(data))
    x = np.reshape(data[:,0],(1,side_length,side_length))    
    y = np.reshape(data[:,1],(1,side_length,side_length))
    z = np.reshape(data[:,2],(1,side_length,side_length))
    coordinates = np.concatenate((x,y,z), axis=0)    
    potential = np.reshape(data[:,3],(side_length,side_length))    
    return coordinates, potential    


