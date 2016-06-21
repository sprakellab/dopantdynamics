"""Python script to prepare for bond parameter analysis 
   with the script from Lechner"""

import numpy as np
import os
import subprocess

import catpy.basisfunctions as bf
import catpy.IO.xml as xml

def prepareforbp(Input, base=False):
    prepareposfiles(Input, base)
    bpdir = '%sbondparameters/' % (Input['systempath'])
    xyzdir = '%sxyz_bp/' % (Input['systempath'])
    os.mkdir(bpdir)
    os.mkdir(xyzdir)

def prepareposfiles(Input, base):
    bp_folder = '%spositions_bp/' % Input['systempath']
    os.mkdir(bp_folder)
    filelist = bf.lookupfiles(path=Input['folder_snapshot_data'], filechr='Wigner')
    for file_ in filelist:
        filename = '%s%s' % (Input['folder_snapshot_data'], file_)
        f = open(filename, 'r')
        xml.ReadHeader(f) 
        if base is True:
            positions = xml.ReadCoordinates(f, Input['BASEparticles'])
        else:
            positions = xml.ReadCoordinates(f, Input['ALLparticles'])
        f.close()
        bp_positions = positions + 0.5*Input['BOXsize']
        savename = '%s%s.dat' % (bp_folder, file_.strip('.xml'))
        saveasdatfile(Input, bp_positions, savename, base)

def saveasdatfile(Input, bp_positions, savename, base):
    s = open(savename, 'w')
    if base is True:
        s.write('%d\n' % Input['BASEparticles'])    
    else:
        s.write('%d\n' % Input['ALLparticles'])
    s.write('%f\n' % Input['BOXsize_x'])
    s.write('%f\n' % Input['BOXsize_y'])
    s.write('%f\n' % Input['BOXsize_z'])
    for i in range(len(bp_positions)):
        s.write('%f %f %f\n' % (bp_positions[i,0], bp_positions[i,1], bp_positions[i,2]))  
    s.close()

def saveinputasdat(Input):
    filename = '%sInputConfiguration.xml' % (Input['systempath']) 
    f = open(filename, 'r')
    xml.ReadHeader(f)
    positions = xml.ReadCoordinates(f, Input['ALLparticles'])
    f.close()
    bp_positions = positions + 0.5*Input['BOXsize']
    savename = '%sInputConfiguration.dat' % (Input['systempath'])
    saveasdatfile(Input, bp_positions, savename)  

def calculate_bondparameter(Input):
    system = Input['system']
    neighbourdistance = str(Input['a'])
    args = ['./scripts/StructureAnalysis/main', system, neighbourdistance]
    p = subprocess.Popen(args)
    p.wait()

def averagebondparameter(Input, base=False):
    dir_bp = '%sbondparameters/' % Input['systempath']
    filelist = bf.lookupfiles(path=dir_bp, filechr='Wigner')
    if base is True:
        average_data = np.zeros((Input['BASEparticles'], 6))
    else:
        average_data = np.zeros((Input['ALLparticles'], 6))
    columns = np.array([0,4,5,6,7,8])
    count = 0
    for file_ in filelist:
        number = int(file_.strip('Wigner.').strip('.dat'))

        if number > Input['analysisstart']:
            filename = '%s%s' % (dir_bp, file_)
            data = np.genfromtxt(filename, delimiter=' ')
            average_data += data[:,columns]
            count +=1
    average_data = average_data/count
    savename = '%s%saveragebp.txt' % (Input['folder_processed_data'], Input['system']) 
    np.savetxt(savename, average_data, delimiter='\t')

def bondparameter_histogram(Input, bins_=100):
    filename = '%s%saveragebp.txt' % (Input['folder_processed_data'], Input['system'])
    data = np.genfromtxt(filename, delimiter='\t')
    hist,binedges = np.histogram(data[:,2],bins=bins_, range=(0.0,0.7))

    combinedhistograms = np.hstack((np.reshape(binedges[:-1],(-1,1)), np.reshape(hist,(-1,1))))
    savename = '%s%sbondparameter.hist' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, combinedhistograms, delimiter='\t')

def bondparameter_distance(Input):
    filename = '%s%saveragebp.txt' % (Input['folder_processed_data'], Input['system'])
    data = np.genfromtxt(filename, delimiter='\t')
    filename = '%s%s' % (Input['folder_snapshot_data'], 'Wigner.0099000000.xml')
    f = open(filename, 'r')
    xml.ReadHeader(f) 
    if base is True:
        positions = xml.ReadCoordinates(f, Input['BASEparticles'])
    f.close()
