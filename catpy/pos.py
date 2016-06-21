""" 
INPUT OUTPUT MODULE FOR POSITION FILES
This module contains functions for analyzing data based on particle positions
.pos files contain tab separated positions of the particles.
On each line the step is given and the x,y,z coordinates follow.
"""

import numpy as np
from time import sleep
import os
import shutil

import catpy.IO.xml as xml

def ReadPos(Nsave,f):
	"""
	This function reads a pos file.
	Input: 	Number of saved steps (Nsave)
		Directions to the file (f)
	Output:	Array with positions (positions)
		Array with the steps (steps)
	"""
	positions = np.zeros([Nsave,3])
	steps = np.zeros([Nsave])
	for i in range(0,Nsave):
		line = f.readline().split('\t')
		steps[i] = int(line[0])
		positions[i,0] = float(line[1])
		positions[i,1] = float(line[2])
		positions[i,2] = float(line[3])
	return positions, steps

def ExtractInterstitials(positions, Nbase, Ninter):
	if len(positions) == Nbase+Ninter:	
		interstitials = positions[Nbase:,:]
	else:
		raise ErrorValue('number of positions(' + str(len(positions))
		+ ') does not match Nbase + Ninter (' + str(Nbase+Ninter) + ').' )
	return interstitials

def SaveInterstitials(interstitials,system,filenumber):
	for i in range(len(interstitials)):	
		InterstitialString = str(i+1)
		n = (3-len(InterstitialString)) % 3
		InterstitialNumber = '0'*n + InterstitialString 
		SaveName = './DATA/' + system + '/interstitial/' + \
			'INTER' + InterstitialNumber + '.pos' 
		SaveFile = open(SaveName, 'a+')
		SaveFile.write( str(filenumber) + '\t' + \
				str(interstitials[i,0]) + '\t' + \
				str(interstitials[i,1]) + '\t' + \
				str(interstitials[i,2])  + '\n')
		SaveFile.close()

def ParticleTrackingLive(system,Nsteps=-1):
	# COLLECTING DATA
#	InitialSettings = xml.ReadInitialFile(System)
	
	Nbase = 14859
	Ninter = 866	
	
	path = './DATA/' + system
	BASEposition = path + '/position'
	INTERposition = path + '/interstitial'
	os.mkdir(BASEposition)
	os.mkdir(INTERposition)
	print('Looking into: ' + path)
	while True:
		sleep(5)	
		filelist=[]		
		#filelist = os.listdir(path)
		for File in files(path):
			filelist.append(File)
		filelist.sort()
		filelist

		print(len(filelist)-3)
		for i in filelist:
			FileNumber = i.strip('Wigner').strip('.xml')
			if not FileNumber.isdigit():
				continue
			FileNumber = int(FileNumber)
			filename = path + '/' + i
			F = open(filename, 'r')
			positions = xml.ReadShort(F)
			F.close()
			interstitials = ExtractInterstitials(positions,Nbase,Ninter)
			SaveInterstitials(interstitials,system,FileNumber)
			if (FileNumber%5000) == 0:
				newname = BASEposition + '/' + i
				os.rename(filename,newname)
			else:
				os.remove(filename)
			if FileNumber == Nsteps:
				return

def unfold(Input, pos):
	"""
	Unfolds a particle trajectory. This can be nescesarry when the particle has passed beyond its bounds.
	Input:	An array with posiitions of the particle with the simulation steps 
		in the first column and x,y,z in 2,3,4 respectively (pos)
		A dictionary with all the initial parameters of the system (Input)
	Output:	An array with the relative positions with respect to the starting point over time.
		The first column are the recorded steps and 2,3 and 4 are the x,y and z coordinate (unfold)
	"""
	posshift = pos[:,1:] + 0.5*Input['BOXsize']
	posdis = (posshift[1:,:] - posshift[:-1,:])
	unfold = np.zeros([len(pos),3])
	for i in range(len(posdis)):
		for j in range(3):
			if posdis[i,j] > 0.5*Input['BOXsize']:			
				unfold[i+1,j] = unfold[i,j] + posdis[i,j] - Input['BOXsize']
			elif posdis[i,j] < -0.5*Input['BOXsize']:
				unfold[i+1,j] = unfold[i,j] + posdis[i,j] + Input['BOXsize']
			else:
				unfold[i+1,j] = unfold[i,j] + posdis[i,j]
	steps = pos[:,0].reshape(len(pos),1)
	unfoldtot = np.concatenate((steps,unfold),axis=1)
	return unfoldtot

def files(path):
    for File in os.listdir(path):
        if os.path.isfile(os.path.join(path, File)):
            yield File	

def CopyFiles(system):
	folderpath = './DATA/TestInterstitial'
	List = []
	for File in files(folderpath):
		List.append(File)
	List.sort()
	counter = 0
	for i in List:	
		sleep(1)
		old = folderpath + '/' + i
		new = './DATA/' + system + '/' + i
		shutil.copy(old,new)
