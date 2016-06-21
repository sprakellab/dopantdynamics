"""Calculate mean square displacement for particle tracks."""

import logging
import os
import numpy as np
import datetime
import multiprocessing
from lmfit import Model
from scipy.spatial.distance import cdist

import catpy.IO.xml as xml
import catpy.basisfunctions as bf

def msd_wholefolder(Input, nprocs=25):
    """Calculates the msd profiles of a folder with tracks"""

    def worker(Input, filelist, out_q):
        outdict = {}
        for filename in filelist:
            outdict[i] = calculate_msd(Input, filename)			
        out_q.put(outdict)

    os.mkdir(Input['folder_msd_data'])
    filelist = bf.lookupfiles(path=Input['folder_interstitial_track'], filechr='INTER')
    out_q = multiprocessing.Queue()
    chunksize = int(np.ceil(len(filelist) / float(nprocs)))
    procs = []
    for i in range(nprocs):
        p = multiprocessing.Process(target=worker, args=(Input,filelist[chunksize*i:chunksize*(i + 1)], out_q))
        procs.append(p)
        p.start()
    for p in procs:
        p.join()
    logging.info('finished calculating all independent msd profiles.')

def calculate_msd(Input, filename):
    """Gathering tracks, calling msd, saving msd"""
    path = '%s%s' % (Input['folder_interstitial_track'], filename)
    rawdata = np.genfromtxt(path, delimiter='\t')
    data = rawdata[(Input['analysisstart']/Input['XMLoutputperiod']):]
    msdarray = msd(Input, data)
    savepath = '%s%smsd' % (Input['folder_msd_data'], filename.strip('pos'))
    np.savetxt(savepath, msdarray, delimiter='\t')

def calculate_msd_corr(Input, filename):
    """Gathering tracks, calling msd, saving msd"""
    path = '%s%s' % (Input['folder_interstitial_track'], filename)
    rawdata = np.genfromtxt(path, delimiter='\t')
    data = rawdata[(Input['analysisstart']/Input['XMLoutputperiod']):]
    msdarray = msd_corr(Input, data)
    savepath = '/home/hoomd/Desktop/%smsd' % (filename.strip('pos'))
    np.savetxt(savepath, msdarray, delimiter='\t')

def msd(Input, data, nstep=200):
    """Core for msd calculation of particle track"""
    nlogspace = nstep - 48
    index_scalarrange = np.arange(1,49)
    index_logrange = np.logspace(np.log10(50),np.log10(len(data)-1),num=nlogspace,dtype=int)
    index_range = np.hstack((index_scalarrange, index_logrange))
    msd = np.zeros(nstep)
    for i in range(nstep):
        index = int(index_range[i])
        r_0 = data[:-index,1:]
        r = data[index:,1:]
        dr = r - r_0
        dr_sq = np.sum(np.square(dr), axis=1)
        msd[i] = np.average(dr_sq)
    time = index_range*Input['XMLoutputperiod']*Input['dt']
    time = np.reshape(time,(len(time),1))
    msd = np.reshape(msd,(len(msd),1))
    dataarray = np.concatenate((time, msd), axis=1)
    return dataarray

def msd_corr(Input, data, nstep=200):
    """Core for msd calculation of particle track"""
    nlogspace = nstep - 48
    index_scalarrange = np.arange(1,49)
    index_logrange = np.logspace(np.log10(50),np.log10(len(data)-1),num=nlogspace,dtype=int)
    index_range = np.hstack((index_scalarrange, index_logrange))
    msd = np.zeros(nstep)
    for i in range(nstep):
        index = int(index_range[i])
        r_0 = data[:-index,1:]
        t_0 = data[:-index,0]
        r = data[index:,1:]
        t = data[index:,0]
        dr = r - r_0
        dt = t - t_0
        t_select = np.nonzero(np.equal(dt, index*Input['XMLoutputperiod']))
        print str(len(dt)) + '\t'+  str(len(dr[t_select]))
        dr_sq = np.sum(np.square(dr[t_select]), axis=1)
        msd[i] = np.average(dr_sq)
    time = index_range*Input['XMLoutputperiod']*Input['dt']
    time = np.reshape(time,(len(time),1))
    msd = np.reshape(msd,(len(msd),1))
    dataarray = np.concatenate((time, msd), axis=1)
    return dataarray

def meanmsd_wholefolder(Input):
    filelist = bf.lookupfiles(path=Input['folder_msd_data'], filechr='INTER')
    path = '%s%s' % (Input['folder_msd_data'], filelist[0])
    data = np.genfromtxt(path, delimiter='\t')        
    total_msd = np.zeros(len(data))
    mean_msd = np.zeros((len(data),2))
    mean_msd[:,0] = data[:,0]
    for filename in filelist:
        path = '%s%s' % (Input['folder_msd_data'], filename)
        data = np.genfromtxt(path, delimiter='\t')
        total_msd += data[:,1]
    mean_msd[:,1] = total_msd / len(filelist)
    savepath = '%s%smean.msd' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savepath, mean_msd, delimiter='\t')

def gradient_msd(Input, data):
    gradient = (data[1:,1] - data[:-1,1]) / (data[1:,0] - data[:-1,0])
    average_gradient = np.zeros(np.shape(data))
    average_gradient[:,0] = data[:,0]
    average_gradient[0,1] = gradient[0]
    average_gradient[-1,1] = gradient[-1]
    average_gradient[1:-1,1] = (gradient[1:] + gradient[:-1]) / 2.0
    return average_gradient

def diffusion_msd(Input):
    filename = '%s%smean.msd' % (Input['folder_processed_data'], Input['system'])
    data = np.genfromtxt(filename, delimiter='\t')
    gradient = (data[1:,1] - data[:-1,1]) / (data[1:,0] - data[:-1,0])
    average_gradient = np.zeros(len(data))
    average_gradient[0] = gradient[0]
    average_gradient[-1] = gradient[-1]
    average_gradient[1:-1] = (gradient[1:] + gradient[:-1]) / 2.0
    diffusion = np.zeros(np.shape(data))
    diffusion[:,0] = data[:,0]
    diffusion[:,1] = (1.0/6.0)*average_gradient
    savename = '%s%smean.diff' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, diffusion, delimiter='\t')
    return diffusion

def fit_D_phi(phi_, D, c1, c2):
    """fit stuff"""
    Dmod = Model(Dhop)
    params = Dmod.make_params()
    params['c1'].set(c1, vary=False)
    params['c2'].set(c2, vary=False)
    params['c3'].set(1, vary=False)
    params['beta'].set(0.70, min=0)
    result = Dmod.fit(D, params, phi=phi_)
    print result.fit_report(min_correl=0.5)
    beta = result.params['beta'].value
    c3 = result.params['c3'].value
    return c3, beta

def Ea_theory(phi, c1, c2):
    epsilon = 470
    kappa = 1.8
    a = (np.pi/(3.0*phi))**(1.0/3.0)
    return ((epsilon*np.exp(-kappa*c1*a)) / (c1*a)) - ((epsilon*np.exp(-kappa*c2*a)) / (c2*a))

def Dhop(phi, c1, c2, c3, beta):
    Ea = Ea_theory(phi, c1, c2)
    D0 = 1.0/(0.25*6.0*np.pi*1.06103295e-01)
    Dhop = c3*D0*np.exp(-((Ea)**beta))
    return Dhop

def Dhopalt(phi, c1, c2, c1o, c2o, c3):
    Eat = Ea_theory(phi, c1, c2)
    Eao = Ea_theory(phi, c1o, c2o)
    D0 = 1/(0.25*6*np.pi*1.06103295e-01)
    Dhop = c3*D0*(8.0*np.exp(-(Eat)**4.0) + 4.0*np.exp(-(Eao)**2.0))
    return Dhop

def fit_D_phi_alt(phi_, D, c1, c2, c1o, c2o):
    """fit stuff"""
    Dmod = Model(Dhopalt)
    params = Dmod.make_params()
    params['c1'].set(c1, vary=False)
    params['c2'].set(c2, vary=False)
    params['c1o'].set(c1o, vary=False)
    params['c2o'].set(c2o, vary=False)
    params['c3'].set(0.25)
    result = Dmod.fit(D, params, phi=phi_)
    print result.fit_report(min_correl=0.5)
    c3 = result.params['c3'].value
    return c3

def fit_D_phi_powerlaw(phi_, D):
    """fit stuff"""
    Dmod = Model(Dhop_powerlaw)
    params = Dmod.make_params() 
    params['p'].set(0.3)
    result = Dmod.fit(D, params, phi=phi_)
    print result.fit_report(min_correl=0.5)
    p = result.params['p'].value
    return p

def Dhop_powerlaw(phi, p):
    return ((phi)**(p))
 
def get_msd_particles(Input, tau_inf=179):
    filelist = bf.lookupfiles(path=Input['folder_msd_data'], filechr='INTER')
    msd_particle = np.zeros([len(filelist), 2])
    counter = 0
    for particle in filelist:
        filename = '%s%s' % (Input['folder_msd_data'], particle)
        data = np.genfromtxt(filename, delimiter='\t')
        msd_particle[counter,0] = int(particle.strip('.msd').strip('INTER'))
        msd_particle[counter,1] = data[tau_inf,1]
        counter += 1
    msd_particle = msd_particle[msd_particle[:,0].argsort()]
    savename = '%s%smsdparticles.val' % (Input['folder_processed_data'], Input['system'])
    format_ = ['%6d','%.12e']
    np.savetxt(savename, msd_particle, delimiter='\t', fmt=format_)    

def p_msd_particles(Input, bins_=30):
    filename = '%s%smsdparticles.val' % (Input['folder_processed_data'], Input['system'])
    data = np.genfromtxt(filename, delimiter='\t')
    histogram = np.zeros([bins_,2])
    phist,pbinedges = np.histogram(data[:,1],bins=bins_, range=(0,5))
    histogram[:,1] = phist
    histogram[:,0] = pbinedges[:-1]
    savename = '%s%smsdparticles.hist' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, histogram, delimiter='\t')       

def Dinf_model(time, var1, Dinf):
    return var1**(1/time) -1 + Dinf

def fit_Dinf_model(Input, diffusion, startindex, endindex):
    time_ = diffusion[startindex:endindex,0]
    D = diffusion[startindex:endindex,1]
    Dmod = Model(Dinf_model)
    params = Dmod.make_params()
    params['var1'].set(0.1) 
    params['Dinf'].set(0.0001)#, min = 0, max=0.1)
    result = Dmod.fit(D, params, time=time_)
    print result.fit_report(min_correl=0.5)
    var1 = result.params['var1'].value
    Dinf = result.params['Dinf'].value
    Dinfstderr = result.params['Dinf'].stderr
    return var1, Dinf, Dinfstderr

def fitting_procedure_Dinf_model(Input, diffusion, startindex, endindex):
    var1, Dinf, Dinfstderr = fit_Dinf_model(Input, diffusion, startindex, endindex)
    save_Dinf_fit_parameters(Input, var1, Dinf, Dinfstderr, startindex, endindex)
    return var1, Dinf, Dinfstderr

def save_Dinf_fit_parameters(Input, var1, Dinf, Dinfstderr, startindex, stopindex):    
    savename = '%s%sDinf.val' % (Input['folder_processed_data'], Input['system'])
    f = open(savename, 'w')
    f.write('%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\n' % ('system', 'phi[-]','startindex', 'stopindex', 'var1', 'Dinf[sigma^2/tau]', 'Dinf_std'))
    f.write('%12s\t%12.8e\t%12.8e\t%12.8f\t%12.8f\t%12.8e\t%12.8e\n' % (Input['system'],Input['BASEfraction'],startindex, stopindex, var1, Dinf, Dinfstderr))
    f.close()

def collect_Dinf(list_):
    data = np.zeros([len(list_),6])
    i = 0
    for system in list_:
        filename = './DATA/%s/processeddata/%sDinf.val' % (system, system)
        f = open(filename, 'r')
        f.readline()
        line = f.readline().split('\t')
        data[i,0] = float(line[1])
        data[i,1] = float(line[2])
        data[i,2] = float(line[3])
        data[i,3] = float(line[4])
        data[i,4] = float(line[5])
        data[i,5] = float(line[6])
        i += 1
    return data

def collect_Dtau(list_):
    data = np.zeros([len(list_),3])
    i = 0
    for system in list_:
        filename = './DATA/%s/processeddata/%sDtau.val' % (system,system)
        f = open(filename, 'r')
        f.readline()
        line = f.readline().split('\t')
        data[i,0] = float(line[1])
        data[i,1] = float(line[2])
        data[i,2] = float(line[3])
        i += 1
    return data

def find_Dtau(Input, diffusion, tau=179):
    Dtau = diffusion[tau,1]
    save_Dtau(Input, tau, Dtau)

def save_Dtau(Input, tau, Dinf):
    savename = '%s%sDtau.val' % (Input['folder_processed_data'], Input['system'])
    f = open(savename, 'w')
    f.write('%-12s\t%-12s\t%-12s\t%-12s\n' % ('system', 'phi[-]','tau[#]', 'Dtau[sigma^2/tau]'))
    f.write('%12s\t%12.8e\t%12.8e\t%12.8f\n' % (Input['system'],Input['BASEfraction'],tau, Dinf))
    f.close()

def distancetointerstitials(Input):
    """Calculates the distance to the nearest interstitial for every particle"""
    # Does take into account the boundaries
    filename = '%sWigner.0099000000.xml' % (Input['folder_snapshot_data']) 
    f = open(filename, 'r')
    xml.ReadHeader(f)
    positions = xml.ReadCoordinates(f, Input['ALLparticles'])
    f.close()
    count = 0
    list_ = [-1.0,0.0,1.0]
    N_inter = np.ones((Input['BASEparticles'],Input['INTERparticles']))*np.arange(1,Input['INTERparticles']+1)
    for x in list_:
        for y in list_:
            for z in list_:
                transpose = np.array([x,y,z])*Input['BOXsize']
                transpose_array = np.ones([Input['BASEparticles'],3])*transpose
                base_transpose = positions[:Input['BASEparticles']] + transpose_array
                distance_sub = cdist(base_transpose, 
                                  positions[Input['BASEparticles']:], 
                                  'euclidean')
                if count == 0:
                    distance_collection = distance_sub
                    selection_collection = N_inter

                else:
                    distance_collection = np.hstack((distances, distance_sub))
                    selection_collection = np.hstack((selection, N_inter))
                distances = np.reshape(np.min(distance_collection, axis=1),(-1,1)) 
                selection_index = np.nonzero(np.equal(distance_collection,distances))
                selection = np.reshape(selection_collection[selection_index],(-1,1))

                count += 1
    data = np.hstack((selection,distances))
    savename = '%s%snearestinterstitial.dist' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, data, delimiter='\t')

def displacementtrack(Input, name='INTER007'):    
    """Gives a displacement track of a particle"""
    savename = '%s%s.pos' % (Input['folder_interstitial_track'],name)
    rawdata = np.genfromtxt(savename, delimiter='\t')
    data = rawdata[(Input['analysisstart']/Input['XMLoutputperiod']):]
    track = np.zeros([len(data),2])
    track[:,0] = data[:,0]
    track[:,1] = np.sqrt(np.sum(np.square((data[:,1:] - data[0,1:])),axis=1))
    savename = '%s%s%sdisplacementtrack.txt' % (Input['folder_processed_data'], Input['system'],name)
    np.savetxt(savename, track, delimiter='\t')

def displacementtrack_experimental(Input, name='trajectory_nr_50.csv'):    
    """Gives a displacement track of a particle"""
    savename = '%s%s' % (Input['folder_interstitial_track'],name)
    rawdata = np.genfromtxt(savename, delimiter=',')
    data = rawdata[:] / 1.8
    track = np.zeros([len(data),2])
    track[:,0] = np.linspace(0,(len(data)-1)*0.2/Input['BASEtau'],len(data))
    track[:,1] = np.sqrt(np.sum(np.square((data[:,:] - data[0,:])),axis=1))
    savename = '%s%s%sdisplacementtrack.txt' % (Input['folder_processed_data'], Input['system'],name)
    np.savetxt(savename, track, delimiter='\t')

def nearvsq6average(Input, points=50):
    filename1 = '%s%snearestinterstitial.dist' % (Input['folder_processed_data'], Input['system'])
    filename2 = '%sbondparameters/Wigner.0099000000.dat' % (Input['systempath'])
    nearest = np.genfromtxt(filename1, delimiter='\t')
    bondparameter = np.genfromtxt(filename2, delimiter=' ')
    maxr = np.max(nearest[:,1])
    linear = np.reshape(np.linspace(0, maxr*1.01, points),(1,-1))
    space = np.repeat(linear, len(nearest),axis=0)
    matrix = np.repeat(np.reshape(nearest[:,1], (-1,1)), points-1,axis=1)
    filtered = np.less(matrix,space[:,1:])*np.greater(matrix,space[:,:-1])
    average = np.zeros([points-1,1])
    for i in range(points-1):
        mask = np.nonzero(filtered[:,i])
        average[i] = np.nan_to_num(np.average(bondparameter[mask,5]))
    data = np.hstack((linear.T[:-1],average))
    savename = '%s%saverageq6vsnearestdopant.txt' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, data, delimiter='\t')

def nearvsdlaverage(Input, points=50):
    filename1 = '%s%snearestinterstitial.dist' % (Input['folder_processed_data'], Input['system'])
    filename2 = '%s%s.deviations' % (Input['folder_processed_data'], Input['system'])
    nearest = np.genfromtxt(filename1, delimiter='\t')
    deviations = np.genfromtxt(filename2, delimiter='\t', skip_header=1)
    maxr = np.max(nearest[:,1])
    linear = np.reshape(np.linspace(0, maxr*1.01, points),(1,-1))
    space = np.repeat(linear, len(nearest),axis=0)
    matrix = np.repeat(np.reshape(nearest[:,1], (-1,1)), points-1,axis=1)
    filtered = np.less(matrix,space[:,1:])*np.greater(matrix,space[:,:-1])
    average = np.zeros([points-1,1])
    for i in range(points-1):
        mask = np.nonzero(filtered[:,i])
        average[i] = np.nan_to_num(np.average(deviations[mask,7]))
    data = np.hstack((linear.T[:-1],average))
    savename = '%s%saveragedlvsnearestdopant.txt' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, data, delimiter='\t')
