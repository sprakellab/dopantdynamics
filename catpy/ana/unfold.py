"""All functions related to the H-algorithm, including the functions
for using multiple CPU's for processing the data"""

import os
import multiprocessing
import logging

import numpy as np

import catpy.IO.pos as IOP
import catpy.basisfunctions as bf

def folder_unfolder(Input, nprocs=15):
    """Unfoldes tracks in folder if barrier is crossed"""

    def worker(Input, filelist, out_q):
        outdict = {}		
        for filename in filelist:
            outdict[i] = track_unfolder(Input, filename)	
        out_q.put(outdict)

    os.mkdir(Input['folder_backup'])
    filelist = bf.lookupfiles(path=Input['folder_interstitial_track'], filechr='INTER')
    out_q = multiprocessing.Queue()
    chunksize = int(np.ceil(len(filelist) / float(nprocs)))
    procs = []
    for i in range(nprocs):
        p = multiprocessing.Process(target=worker, args=(Input, filelist[chunksize*i:chunksize*(i + 1)], out_q))
        procs.append(p)
        p.start()
    for p in procs:
        p.join()
    logging.info('finished calculating all harrays.')

def track_unfolder(Input, filename):
    path = '%s%s' % (Input['folder_interstitial_track'], filename)
    backuppath = '%s%s' % (Input['folder_backup'], filename)
    rawdata = np.genfromtxt(path, delimiter = '\t')
    test = bf.testoutofbounds(Input, rawdata)	
    if test is True:
        logging.info('interstitial %s out of bounds, unfolding...' % filename) 
        unfoldeddata = np.zeros([len(rawdata),3])		
        unfoldeddata = IOP.unfold(Input, rawdata)
        logging.info('interstitial %s is unfolded succesfully!' % filename)
        os.rename(path, backuppath)
        np.savetxt(path, unfoldeddata, delimiter='\t')
