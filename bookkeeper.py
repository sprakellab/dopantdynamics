"""This script helps in handling the output files of hoomd
It is a seperate script as it has to work on the CPU's and not the GPU"""

import numpy as np
import multiprocessing
import os
import sys
import time

import catpy.init.parameters as initcolloid

# get commandline arguments
if not (len(sys.argv)==2):
    print 'usage bookkeeper.py: python bookkeeper.py system'
    exit(0)

#filenames
system = sys.argv[1]

def logstdout_on(Input):
	"""
	Switches the stdout from the terminal to a file in the BOOKKEEPER folder.
	Input:	The dict with all the initial parameters (Input)
	Output: The directions to the terminal (saveout)
		The directions to the file in BOOKKEEPER (out)
	"""
	saveout = sys.stdout
	outfilename = './BOOKKEEPER/' + Input['system'] + '.out'                                  
	out = open(outfilename, 'w')                             
	sys.stdout = out                                      
	return saveout,out

def logstdout_off(Input,saveout,out):
	"""
	Switches the stdout from a file in the BOOKKEEPER folder to the terminal.
	Input:	The dict with all the initial parameters (Input)
		The directions to the terminal (saveout)
		The directions to the file in BOOKKEEPER (out)
	Output: - 
	"""
	sys.stdout = saveout                                     
	out.close()        

def logstderr_on(Input):
	"""
	Switches the stderr from the terminal to a file in the BOOKKEEPER folder.
	Input:	The dict with all the initial parameters (Input)
	Output: The directions to the terminal (saveout)
		The directions to the file in BOOKKEEPER (err)
	"""
	saveerr = sys.stderr 
	errfilename = './BOOKKEEPER/' + Input['system'] + '.err'                                      
	err = open(errfilename, 'w')                             
	sys.stderr = err  
	return saveerr,err

def logstderr_off(Input,saveerr,err):
	"""
	Switches the stderr from a file in the BOOKKEEPER folder to the terminal.
	Input:	The dict with all the initial parameters (Input)
		The directions to the terminal (saveout)
		The directions to the file in BOOKKEEPER (err)
	Output: - 
	"""
	sys.stderr = saveerr                                     
	err.close()  

def ExtractInterstitials(positions):
	"""
	Extracts the positions of the interstitials from a list with positions of all the particles.
	If the positions list does not contain the right amount of particles, 
	positions are not extracted and the output for the positions is all zeros
	Input:	A list with all the positions of all particles (positions)
	Output: A list with the positions of al interstitials (interstitial)
		If the position list is to short this value is 1 (lengtherror)	
	"""
	if len(positions) == INIT['ALLparticles']:	
		interstitials = positions[INIT['BASEparticles']:,:]
		lengtherror = 0
	else:
		length = len(positions)
		sys.stderr.write('number of positions(' + str(length) + ') does not match Nbase + Ninter (' + str(INIT['ALLparticles']) + ').')
		interstitials = 'error' 	
		lengtherror = 1
	return interstitials,lengtherror

def ExtractBase(positions,randombaseparticles):
	"""
	Extracts the positions of several base particles from a list with positions of all the particles.
	If the positions list does not contain the right amount of particles, 
	positions are not extracted and the output for the positions is all zeros
	Input:	A list with all the positions of all particles (positions)
	Output: A list with the positions of the selected bace particles (interstitial)
		If the position list is to short this value is 1 (lengtherror)	
	"""
	if len(positions) == INIT['ALLparticles']:	
		baseparticles = positions[randombaseparticles,:]
		lengtherror = 0
	else:
		length = len(positions)
		sys.stderr.write('number of positions(' + str(length) + ') does not match Nbase + Ninter (' + str(INIT['ALLparticles']) + ').')
		baseparticles = 'error' 	
		lengtherror = 1
	return baseparticles,lengtherror

def SaveInterstitials(interstitials,filenumber):
	"""
	Saves the positions of the interstitials in files.
	It takes a list with interstitials from one step and devides them over files
	per interstitial;
	Input:	A list with all the positions of all particles (positions)
		The step number extracted from the filename (filenumber)
	Output: - 	
	"""
	for i in range(len(interstitials)):	
		InterstitialString = str(i+1)
		n = (3-len(InterstitialString)) % 3
		InterstitialNumber = '0'*n + InterstitialString 
		SaveName = './DATA/' + INIT['systemfolder'] + '/interstitial/' + \
			'INTER' + InterstitialNumber + '.pos' 
		SaveFile = open(SaveName, 'a+')
		SaveFile.write( str(filenumber) + '\t' + \
				str(interstitials[i,0]) + '\t' + \
				str(interstitials[i,1]) + '\t' + \
				str(interstitials[i,2])  + '\n')
		SaveFile.close()

def SaveBase(baseparticles,randombaseparticles,filenumber):
	"""
	Saves the positions of the base particles in files.
	It takes a list with bases from one step and devides them over files
	per base particles;
	Input:	A list with all the positions of all particles (positions)
		A list of the indeces of which we track the particles (randombaseparticles)
		The step number extracted from the filename (filenumber)
	Output: - 	
	"""
	for i in range(len(baseparticles)):	
		BaseString = str(randombaseparticles[i]+1)
		n = (5-len(BaseString)) % 5
		BaseNumber = '0'*n + BaseString 
		SaveName = './DATA/' + INIT['systemfolder'] + '/base/' + \
			'BASE' + BaseNumber + '.pos' 
		SaveFile = open(SaveName, 'a+')
		SaveFile.write( str(filenumber) + '\t' + \
				str(baseparticles[i,0]) + '\t' + \
				str(baseparticles[i,1]) + '\t' + \
				str(baseparticles[i,2])  + '\n')
		SaveFile.close()

def ParticleTrackingLive(i,BASEposition,INTERposition):
	"""
	Extracts the interstitials position from the xml output by HOOMD, 
	while the simulation is running. Furthermore it saves the positions 
	in seperate files per interstitial. XML files are either removed or 
	moved to another folder.
	Input:	The filename (i)
		Path to the folder to save XML files (BASEposition)
		Path to the folder to save the interstitial positions (INTERposition)
	Output: The number that is in the filename and which corresponds to the step (Filenumber)
		this value is 1 if an error occured in reading the positions (lengtherror)
		this value is 1 if an error occured in reading the positions (valueerror) 	
	"""
	FileNumber = i.strip('Wigner').strip('.xml')
	if not FileNumber.isdigit():
		return
	FileNumber = int(FileNumber)
	filename = path + '/' + i
	F = open(filename, 'r')
	try:
		positions = np.genfromtxt(F,comments='<',delimiter=' ')
		valueerror = 0
	except ValueError as Verror:
		sys.stderr.write("Error in reading indices for step number: " + str(FileNumber) + '\n')
		sys.stderr.write("Value error: " + str(Verror) + '\n') 		
		valueerror = 1
	F.close()
	try:	
		interstitials,lengtherror = ExtractInterstitials(positions)
	except UnboundLocalError:
		lengtherror = 0	
		interstitials = 'error'	
		pass
	if type(interstitials) != str: 
		try:	
			SaveInterstitials(interstitials,FileNumber)
		except UnboundLocalError:
			pass	
	if (FileNumber%INIT['XMLsaveperiod']) == 0:
		newname = BASEposition + '/' + i
		os.rename(filename,newname)
		save_file = open(newname)
		save_info = []
		for line in save_file:
			save_info += [line]
		save_file.close()
		save_extra_info = all_info[:]
		save_extra_info[5:(INIT['ALLparticles']+5)] = save_info[5:(INIT['ALLparticles']+5)]
		save_extra = open(newname, 'w')
		for i in save_extra_info:
			save_extra.write(i)
		save_extra.close()
	else:
		os.remove(filename)
	if FileNumber == INIT['Nsteps']:
		exit()	
	return FileNumber,lengtherror, valueerror

def ParticleTrackingLiveinclBASE(i,BASEposition,INTERposition):
	"""
	Extracts the interstitials position from the xml output by HOOMD, 
	while the simulation is running. Furthermore it saves the positions 
	in seperate files per interstitial. XML files are either removed or 
	moved to another folder.
	Input:	The filename (i)
		Path to the folder to save XML files (BASEposition)
		Path to the folder to save the interstitial positions (INTERposition)
	Output: The number that is in the filename and which corresponds to the step (Filenumber)
		this value is 1 if an error occured in reading the positions (lengtherror)
		this value is 1 if an error occured in reading the positions (valueerror) 	
	"""
	FileNumber = i.strip('Wigner').strip('.xml')
	if not FileNumber.isdigit():
		return
	FileNumber = int(FileNumber)
	filename = path + '/' + i
	F = open(filename, 'r')
	try:
		positions = np.genfromtxt(F,comments='<',delimiter=' ')
		valueerror = 0
	except ValueError as Verror:
		sys.stderr.write("Error in reading indices for step number: " + str(FileNumber) + '\n')
		sys.stderr.write("Value error: " + str(Verror) + '\n') 		
		valueerror = 1
	F.close()
	try:	
		interstitials,lengtherror = ExtractInterstitials(positions)
		baseparticles,lengtherror = ExtractBase(positions,INIT['randombaseparticles'])
	except UnboundLocalError:
		lengtherror = 0	
		interstitials = 'error'
		baseparticles = 'error' 	
		pass
	if type(interstitials) != str: 
		try:	
			SaveInterstitials(interstitials,FileNumber)
		except UnboundLocalError:
			pass
	if type(baseparticles) != str:
		try:	
			SaveBase(baseparticles,INIT['randombaseparticles'],FileNumber)
		except UnboundLocalError:
			pass		
	if (FileNumber%INIT['XMLsaveperiod']) == 0:
		newname = BASEposition + '/' + i
		os.rename(filename,newname)
		save_file = open(newname)
		save_info = []
		for line in save_file:
			save_info += [line]
		save_file.close()
		save_extra_info = all_info[:]
		save_extra_info[5:(INIT['ALLparticles']+5)] = save_info[5:(INIT['ALLparticles']+5)]
		save_extra = open(newname, 'w')
		for i in save_extra_info:
			save_extra.write(i)
		save_extra.close()
	else:
		os.remove(filename)
	if FileNumber == INIT['Nsteps']:
		exit()	
	return FileNumber,lengtherror, valueerror

def files(path):
	"""
	Checks wether a name in a folder is a file or a directory,
	only passes a name if it is a file.
	Input:	(path)
	Output: - 
	"""
	for File in os.listdir(path):
		if os.path.isfile(os.path.join(path, File)):
			yield File	

def mp_filelistdistributer(filelist, nprocs):
       	""" 
	Distributes the files and the processing among different CPU's. 
	The lenght of the filelist is used to determine the amount of 
	tasks will be handed to every process.Furthermore it keeps track 
	of the amount of occuring errors and the last file processed. 
	Input: 	A list of filenames to process (filelist)
		Number of processors available to directed processes to. (nprocs)
	Output:	The highest processed filenumber (LastStep)
		 Amount of length errors (length of positions) spotted in this cycle (Lengtherrors)
		 Amount of Value errors (errors in reading indexes) spotted in this cycle (Valueerrors)
	"""
	def worker(filelist, out_q):
        	""" 
		The worker function, used to set up the multiprocessing. 
		Input: 	A list of filenames to process (filelist)
			The results are placed in a dictionary that's pushed to a queue. (out_q)
		Output:	- 
		"""
		outdict = {}
		if INIT['BASEtrack'] is True:
			for i in filelist:
				outdict[i] = ParticleTrackingLiveinclBASE(i,BASEposition,INTERposition)	
		else:
			for i in filelist:
				outdict[i] = ParticleTrackingLive(i,BASEposition,INTERposition)
		out_q.put(outdict)

	# Each process will get 'chunksize' nums and a queue to put his out dict into
	out_q = multiprocessing.Queue()
	chunksize = int(np.ceil(len(filelist) / float(nprocs)))
	procs = []

	for i in range(nprocs):
		p = multiprocessing.Process(target=worker, args=(filelist[chunksize * i:chunksize * (i + 1)],out_q))
		procs.append(p)
		p.start()
	
	# collecting results
	resultdict = {}
	for i in range(nprocs):
        	resultdict.update(out_q.get())	

	# Wait for all worker processes to finish
	for p in procs:
		p.join()

	v=list(resultdict.values())		
	try:	
		LastStep = max([x[0] for x in v])
	except ValueError as Verror:
		sys.stderr.write("Value error: " + str(Verror) + '\n') 	
		LastStep = -1
	try:
		Lengtherrors = sum([x[1] for x in v])
	except ValueError as Verror:
		sys.stderr.write("Value error: " + str(Verror) + '\n')
		Lengtherrors = 0
	try:
		Valueerrors = sum([x[2] for x in v])
	except ValueError as Verror:
		sys.stderr.write("Value error: " + str(Verror) + '\n')
		Valueerrors = 0
	return LastStep,Lengtherrors, Valueerrors

def SortInterstitials():
        """ 
	Sorts the interstitial positions for the interstitial files 
	along the steps in the first column. 
	Input: 	- 
	Output:	- 
	"""
	os.makedirs('./DATA/' + INIT['systemfolder'] + '/interstitialsorted/')
	for i in range(INIT['INTERparticles']):	
		InterstitialString = str(i+1)
		n = (3-len(InterstitialString)) % 3
		InterstitialNumber = '0'*n + InterstitialString 
		
		DataName = './DATA/' + INIT['systemfolder'] + '/interstitial/' + \
			'INTER' + InterstitialNumber + '.pos' 		
		SaveName = './DATA/' + INIT['systemfolder'] + '/interstitialsorted/' + \
			'INTER' + InterstitialNumber + '.pos' 
		Data = np.genfromtxt(DataName, delimiter='\t')
		SortedData = Data[Data[:,0].argsort()]
		np.savetxt(SaveName,SortedData,fmt = '%i\t%0.12f\t%0.12f\t%0.12f')

def SortBase():
        """ 
	Sorts the base partilce positions for the interstitial files 
	along the steps in the first column. 
	Input: 	- 
	Output:	- 
	"""
	os.makedirs('./DATA/' + INIT['systemfolder'] + '/basesorted/')
	for i in INIT['randombaseparticles']:	
		BaseString = str(i+1)
		n = (5-len(BaseString)) % 5
		BaseNumber = '0'*n + BaseString
		
		DataName = './DATA/' + INIT['systemfolder'] + '/base/' + \
			'BASE' + BaseNumber + '.pos' 		
		SaveName = './DATA/' + INIT['systemfolder'] + '/basesorted/' + \
			'BASE' + BaseNumber + '.pos' 
		Data = np.genfromtxt(DataName, delimiter='\t')
		SortedData = Data[Data[:,0].argsort()]
		np.savetxt(SaveName,SortedData,fmt = '%i\t%0.12f\t%0.12f\t%0.12f')

#######################################################################
# MAIN
#######################################################################

# GETTING PARAMETERS
initfilename = './INIT/' + system + '.init'
I = open(initfilename, 'r')
INIT = initcolloid.ImportInitialization(I)
I.close()
nprocs = 12

# ERROR PARAMETERS
saveout,out = logstdout_on(INIT)
saveerr,err = logstderr_on(INIT)
valueerrorcount = 0
lengtherrorcount = 0

# MAKING DIRECTORIES
path = './DATA/' + INIT['systemfolder']
BASEposition = path + '/position'
INTERposition = path + '/interstitial'
os.mkdir(BASEposition)
os.mkdir(INTERposition)
if INIT['BASEtrack'] is True:
	BASEtrack = path + '/base'
	os.mkdir(BASEtrack)
	x = np.zeros(INIT['BASEparticles'])
	x[:INIT['INTERparticles']] = 1
	np.random.shuffle(x)
	INIT['randombaseparticles'] = np.nonzero(x)[0]

# Getting vis_info
info = open('./DATA/%s/InitialConfiguration.xml' % (INIT['systemfolder']))
all_info = []
for line in info:
	all_info += [line]
info.close()

time.sleep(5)
print 'BOOKKEEPER: Looking into: ' + path

# THE IMPORTANT PART (SHOULD BE SOME KIND OF LOOP)
LastStep = 0
count = 0
if INIT['divide']:
	condition = INIT['max_simulation_length'] - INIT['XMLoutputperiod']
else:
	condition = INIT['Nsteps'] - INIT['XMLoutputperiod']
while LastStep != condition:
	filelist=[]		
	for File in files(path):
		if File[:6] == 'Wigner': # filter for files which are not xml		
			filelist.append(File)
	LastStep,Lengtherrors,Valueerrors = mp_filelistdistributer(filelist, nprocs)
	lengtherrorcount += Lengtherrors
	valueerrorcount += Valueerrors
	count += 1
	if count % 50 == 0:
		print 'BOOKKEEPER: ' + 'finished round with LastStep ' + str(LastStep) + '(count: ' + str(count) + ')' 
	time.sleep(1)
print '\nBOOKKEEPER: finished bookkeeping of system ' + INIT['systemfolder']
print 'BOOKKEEPER: The amount of length errors:\t' + str(lengtherrorcount) + \
      '\nBOOKKEEPER: The amount of value errors:\t' + str(valueerrorcount)
SortInterstitials()
if INIT['BASEtrack'] is True:
	SortBase()
print 'BOOKKEEPER: finished sorting interstitial files' 

logstdout_off(INIT,saveout,out)
logstderr_off(INIT,saveerr,err)
