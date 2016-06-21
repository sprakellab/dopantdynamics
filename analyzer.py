"""This script performs analysis of data, based on the input in the INIT file.
It is a seperate script as it has to work on the CPU's and not the GPU"""


import os
import sys
import logging

import numpy as np

import catpy.init.parameters as initcolloid
import catpy.IO.xml as xml
import catpy.ana.unfold as hana
import catpy.ana.paircorrelation as pana
import catpy.ana.msd as mana
import catpy.ana.experimental as ae
import catpy.plot.base as bplt
import catpy.plot.energy as eplt

def set_log(Input):
    logname = Input['log_analyzer_path']
    logging.basicConfig(filename=logname, level=logging.INFO)
    logging.basicConfig(format='%(levelname)s:%(message)s', 
                        level=logging.INFO)

def main(system):
    # GETTING PARAMETERS
    initfilename = './INIT/%s.init' % system
    I = open(initfilename, 'r')
    INIT = initcolloid.ImportInitialization(I)
    I.close()

    # ERROR PARAMETERS
    set_log(INIT)
    
    # SETTING UP DIRECTORIES
    os.mkdir(INIT['folder_processed_data'])

    #MULTIPLEPAIRCORRELATION
    logging.info('ANALYZER: starting analysis of the base crystal (system: %s)' % INIT['system'])
    totalnormhist,bin_edges_norm = pana.MultiplePairCorrelationFunction(INIT, binsize=0.05)
    pana.SaveIPD(INIT,bin_edges_norm, totalnormhist)
    bplt.PlotIPD(INIT['savenameIPD'],bin_edges_norm, totalnormhist)
    eplt.PlotOverview(INIT['systemfolder'],bin_edges_norm, totalnormhist, INIT['savenameIPD'])
    logging.info('ANALYZER: finished analyzing the base crystal of system: %s' % INIT['system'])

    if INIT['method'] == 'interstitial':
        if INIT['experimental']:
            logging.info('ANALYZER: Adding noise and extracting tracks...')		
            ae.noise_whole_folder(INIT)
            ae.whole_folder_extract_tracks(INIT)
            logging.info('ANALYZER: Finished adding noise and extracting tracks')
        logging.info('ANALYZER: starting interstitial unfolding...')
        hana.folder_unfolder(INIT)
        logging.info('ANALYZER: finished unfolding tracks.')
        logging.info('ANALYZER: starting msd calculations...')
        mana.msd_wholefolder(INIT)
        mana.meanmsd_wholefolder(INIT)
        logging.info('ANALYZER: finished calculating msd tracks.')
        diffusion = mana.diffusion_msd(INIT)        
        mana.find_Dtau(INIT, diffusion, tau=179)

if __name__ == '__main__':

	# get commandline arguments
	if not (len(sys.argv)==2):
	    print('usage analyzer.py: python analyzer.py system')
	    exit(0)
	system = sys.argv[1]

	main(system)
