"""main script for interstitial modeling.
This script contains the main body for simulating interstitials. Using 
the software package hoomd.blue this script sets up the simulation of 
a colloidal crystal of Soft Spheres (Yukawa potential).Simulation can 
be visualized via VMD.

Type: vmd -hoomd FILENAME
Type: imd connect localhost PORT
"""

__author__ = 'Justin Tauber'
__version__ = '1.0'
__date__ = '1/2/2016'

import os
import subprocess
import multiprocessing
import copy

from hoomd_script import *
import numpy as np

import catpy.IO.xml as xml
import catpy.init.parameters as initcolloid

def setsimulation(Input):
    """Sets the parameters for the simulation in hoomd.

    Keyword arguments:
    Input -- The dict with all the initial parameters
    
    Output: 
    yk -- The reference to the class for the yukawa potential 
    parameters
    bd -- The reference to the class for the brownian dynamics 
    parameters from the integrator
    """
    yk = pair.yukawa(r_cut=Input['rcut'])
    yk.pair_coeff.set('A', 'A', epsilon=Input['BASEepsilon'], kappa=Input['BASEkappa'])    
    if Input['method'] == 'interstitial':
        yk.pair_coeff.set('A', 'B', epsilon=Input['IBepsilon'], kappa=Input['IBkappa'])
        yk.pair_coeff.set('B', 'B', epsilon=Input['INTERepsilon'], kappa=Input['INTERkappa'])
    all = group.all()
    BASE = group.type('A')
    if Input['method'] == 'interstitial':
        INTER = group.type('B')
    integrate.mode_standard(dt=Input['dt'])
    if not Input['basefix']:
        BASEbd = integrate.brownian(group=BASE, T=Input['temperature'],dscale=0, seed=np.random.randint(0,high=1000))	
        BASEbd.set_gamma(BASE, gamma=Input['BASEfriction'])
    else:
        BASEbd = None
    if Input['method'] == 'interstitial':
        INTERbd = integrate.brownian(group=INTER, T=Input['temperature'],dscale=0, seed=np.random.randint(0,high=1000))
        INTERbd.set_gamma(INTER, gamma=Input['INTERfriction'])
        return yk, BASEbd, INTERbd, BASE, INTER
    return yk, BASEbd, BASE

def exportdata(Input, FullDirectory, Port=22230):
    """Sets the parameters for exporting data from hoomd.
    
    Input:
    FullDirectory -- path to the systemfolder in DATA
    Input -- The dict with all the initial parameters
    Port -- The port used for live streaming of snapshots
    """
    xml = dump.xml(filename = FullDirectory + 'InitialConfiguration.xml', vis=True) 
    xml2 = dump.xml(filename = FullDirectory + 'Wigner', period=Input['XMLoutputperiod'])
    if Input['dcd'] is True:
        dcd = dump.dcd(filename = FullDirectory + 'dump.dcd', period=1000)
    imd_ = False	
    if imd_:	
        analyze.imd(port=Port, period=1000)# setup the IMD server
    AnalyzeFilename = FullDirectory + 'energies.log'
    analyze.log(filename=AnalyzeFilename, 
                quantities=['pair_yukawa_energy', 'potential_energy','temperature'],
                period=1000, 
                header_prefix='#')

def Bookkeeping(Input):
    """Fires up the bookkeeping procedure on CPU"""
    args = ['python', './scripts/bookkeeper.py', Input['system']]
    p = subprocess.Popen(args)
    return p

def Analysis(Input, p):
    """Fires up the analysis procedure on CPU's."""	
    p.wait()
    args = ['python', './scripts/analyzer.py', Input['system']]
    d = subprocess.Popen(args)

def initial(Input, divide=False):
    """Sets the parameters and starts the simulation in hoomd using a 
    given initial configuration.
  
    Keyword arguments:
    Input -- The dict with all the initial parameters (Input)
    
    Output:
    yk -- The reference to the class for the yukawa potential parameters
    bd -- The reference to the class for the brownian dynamics parameters from the integrator
    """
    option.set_msg_file('./LOG/' + Input['system'] + '.plog')	
    Dbox = data.boxdim(Lx=Input['BOXsize_x'], Ly=Input['BOXsize_y'], Lz=Input['BOXsize_z'])
    FullDirectory = './DATA/' + Input['systemfolder'] + '/'
    os.makedirs(FullDirectory)
    if divide:
        ImportFilename = './DATA/%s/position/%s' % (Input['InitialConfiguration'], Input['max_simulation_filename'])
    else:
        ImportFilename = FullDirectory + 'InputConfiguration.xml'
        if Input['InitialConfiguration'] == 'BCC' and Input['basegenerator'] == 'Nparticles':
            xml.InitBCC(Input, FullDirectory)
        elif Input['InitialConfiguration'] == 'FCC' and Input['basegenerator'] == 'Nparticles':
            xml.InitFCC(Input, FullDirectory)
        elif Input['InitialConfiguration'] == 'BCC' and Input['basegenerator'] == 'Lbox':
            xml.InitBCC_experimental(Input, FullDirectory)
    system = init.read_xml(filename=ImportFilename)
    if Input['method'] == 'interstitial':
        yk,BASEbd,INTERbd,BASE,INTER = setsimulation(Input)
    else:
        yk,BASEbd,BASE = setsimulation(Input)
        INTERbd = None
    exportdata(Input, FullDirectory)
    p = Bookkeeping(Input)
    option.set_msg_file('./LOG/' + Input['system'] + '.log')
    if divide:
        run(INIT['max_simulation_length'])
    elif INIT['divide']:
        run(INIT['max_simulation_length'])
    else:
        run(INIT['analysisstart'])
        if Input['msd'] is True:
            msdfilename = FullDirectory + 'msd.log'
            msd = analyze.msd(filename=msdfilename, groups=[BASE, INTER], period=100)
        run(INIT['Nsteps']-INIT['analysisstart'])
#    subprocess.wait()
    Analysis(Input,p)
    return yk,BASEbd,INTERbd

def random(Input):
    """Sets the parameters and starts the simulation in hoomd using a 
    random initial configuration.
    
    Keyword arguments:
    Input -- The dict with all the initial parameters
    yk -- The reference to the class for the yukawa potential 
    parameters
    bd -- The reference to the class for the brownian dynamics 
    parameters from the integrator
    """
    Dbox = data.boxdim(Lx=Input['BOXsize'], Ly=Input['BOXsize'], Lz=Input['BOXsize'])
    FullDirectory = './DATA/' + Input['systemfolder'] + '/'
    os.makedirs(FullDirectory)
    if Input['method'] == 'wigner':
        init.create_random(N=Input['BASEparticles'], name='A', box=Dbox)
    elif Input['method'] == 'interstitial':
        print 'Jammerdebammer, dit stuk is nog niet geschreven'
        return
    yk,bd = setsimulation(Input)
    exportdata(Input, FullDirectory)
    p = Bookkeeping(Input)
    run(Input['Nsteps'])
    Analysis(Input,p)
    return yk,bd

def series():
    """series in volume fraction are used in order to start a series of
    simulations sequentially."""
    PhiArray = np.arange(INIT['BASEfractionStart'],INIT['BASEfractionEnd'],INIT['BASEfractionStep'])
    for Phi in PhiArray:
        SERIES = INIT.copy()
        SERIES['BASEfraction'] = Phi
        SERIES['series'] = False
        SERIES = initcolloid.InitializeCalculations(SERIES)
        SaveInit = './INIT/' + SERIES['system'] + '.init'	
        s = open(SaveInit, 'w')
        initcolloid.ExportInitialization(s,SERIES)
        s.close()
        if SERIES['InitialConfiguration'] == 'random':
            yk,bd = random(SERIES)
            del yk
            del bd
        elif (SERIES['InitialConfiguration'] == 'BCC' or SERIES['InitialConfiguration'] == 'FCC'):
            yk,BASEbd,INTERbd = initial(SERIES)
            del yk
            if BASEbd != None:
                del BASEbd
            if INTERbd != None:
                del INTERbd
        init.reset() 

def divide_single():
    """ Seperates a simulation in multiple blocks. 
    Blocks follow up on each other """
    runrange = np.arange(1,INIT['Nblocks'])
    for run in runrange:
        DIVIDER = copy.deepcopy(INIT)
        DIVIDER['system'] = '%sR%02d' % (INIT['system'],run)
        DIVIDER['systemfolder'] = '%s/%s' % (INIT['system'],DIVIDER['system'])
        DIVIDER['Nsteps'] = INIT['max_simulation_length']
        if run > 1:		
            DIVIDER['InitialConfiguration'] = '%s/%sR%02d' % (INIT['systemfolder'],INIT['systemfolder'],run-1)
            Divide = True
        else:
            Divide = False
        SaveInit = './INIT/' + DIVIDER['system'] + '.init'	
        s = open(SaveInit, 'w')
        initcolloid.ExportInitialization(s,DIVIDER)
        s.close()
        print 'Started simulation %s' % DIVIDER['system']
        yk,BASEbd,INTERbd = initial(DIVIDER,divide=Divide)
        del yk
        if BASEbd != None:
            del BASEbd
        if INTERbd != None:
            del INTERbd
        print 'Finished simulation %s' % DIVIDER['system']
        init.reset()

if __name__ == '__main__':
    Init = raw_input('Enter path to .init or enter "No" to start initialization: ')
    if Init == 'No':
        INIT = {}
        initcolloid.FullInitialization(INIT)
        if INIT['series'] is True:
            SaveInit = './INIT/' + 'V' + str(INIT['version']) + '.init'
        else:
            SaveInit = './INIT/' + INIT['systemfolder'] + '.init'
        f = open(SaveInit, 'w')
        initcolloid.ExportInitialization(f,INIT)
        f.close()
    else:
        f = None
        while f is None:
            try:
                f = open(Init, 'r')
            except IOError:
                Init = raw_input('Path not present, enter other path or "exit": ')
                if Init == 'exit':
                    exit()
        INIT = initcolloid.ImportInitialization(f)
        f.close()

    if INIT['series'] is True:
        series()
    elif INIT['divide'] is True:
        divide_single()
    elif INIT['series'] is False and (INIT['InitialConfiguration'] == 'BCC' or INIT['InitialConfiguration'] == 'FCC'):
        initial(INIT)
    elif INIT['series'] is False and INIT['InitialConfiguration'] == 'random':
        random(INIT)
    else:
        raise ValueError('%s and %s are no proper values for InitialConfiguration or series.' % (INIT['InitialConfiguration'], INIT['series']))
