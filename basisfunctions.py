"""basic functions of the catpy package"""

import os

import numpy as np

import catpy.init.parameters as initcolloid

def files(path):
    """	Checks wether a name in a folder is a file or a directory
    only passes a name if it is a file.
    """
    for File in os.listdir(path):
        if os.path.isfile(os.path.join(path, File)):
            yield File

def lookupfiles(path, filechr):
    """Produces a list of files in path with a characteristic string

    Keyword arguments:
    path -- path to the desired folder.
    filechr -- characteristic string which is at the start of the 
    desired filenames.
    """
    filelist = []
    for File in files(path):
        if File[:len(filechr)] == filechr: 		
            filelist.append(File)
    return filelist

def testoutofbounds(Input, rawdata):
    """Test wether the interstitials goes out of its bounds.
	
    Keyword arguments:
    Input -- The dictionary containing all properties of the system.
    rawdata -- All the positions of an interstitial recorded in the
    simulation.
    
    Output:
    test -- True if the interstitial moved out of it's bounds.
    """
    if Input['experimental'] is True:
        return False	
    difference = (rawdata[(Input['analysisstart']/Input['XMLoutputperiod']):,1:] - 
                  rawdata[((Input['analysisstart']/Input['XMLoutputperiod'])-1):-1,1:])
    absdifference = np.abs(difference)	
    maximum = np.max(absdifference)
    if maximum > 0.95*Input['BOXsize']:
        test = True
    else:
        test = False
    return test

def import_dict(system):
    """Imports the dictionary with all input data."""
    initfilename = './INIT/%s.init' % system
    I = open(initfilename, 'r')
    Input = initcolloid.ImportInitialization(I)
    I.close()
    return Input
