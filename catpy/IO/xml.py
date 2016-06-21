"""
xml Read Write. This script contains all read and write functions 
for xml files used in the Interstitial simulations and dataprocessing.
"""

import numpy as np
import scipy.io

def ConstructPeriodFileName(system, number):
	"""
	A filename is constructed for a .xml file containing positions 
	which is standard output of HOOMD.
	
	Input:	Step number (number)
	Output:	Filename (f)
	"""
	Number = str(number)	
	MissingZeros = 10 - (len(Number) % 10)
	Number = MissingZeros*'0' + Number
	f = './DATA/' + system + '/Wigner.' + Number + '.xml'
	return f

def ReadShort(f):
	"""
	Reads the short .xml files only containing positions.
	This is a quicker method then lineby line reading.
	Input:	file adress (f)
	Output: all positions of base and interstitials (positions)
	"""
	positions = np.genfromtxt(f,comments='<',delimiter=' ')
	return positions

def ReadHeader(f):
	"""
	This function reads the header of a xml file used for logging the positions 
	of particles in a simulation.
	
	Input: 	the filename
	Output: The time step (time_step), 
		the dimensions (dimensions), 
		the number of atoms (natoms), 
		size of the simulation box in the x direction (lx), 
		size of the simulation box in the y direction (ly) 
		size of the simulation box in the z direction (lz).
	"""
	Header3 = []
	Header4 = []	
	line = f.readline()
	line = f.readline()
	line = f.readline()	
	Header3.append(line.strip().strip('<').strip('>').split(' '))
	line = f.readline()	
	Header4.append(line.strip().strip('<').strip('/>').split(' '))
	line = f.readline()
	time_step = int(Header3[0][1].strip('time_step=').strip('"'))
	if len(Header3[0]) > 2:
		dimensions = int(Header3[0][2].strip('dimensions=').strip('"'))	
		natoms = int(Header3[0][3].strip('natoms=').strip('"'))
	else:
		dimensions = None
		natoms = None
	lx = float(Header4[0][1].strip('lx=').strip('"'))
	ly = float(Header4[0][2].strip('ly=').strip('"'))
	lz = float(Header4[0][3].strip('lz=').strip('"'))	
	return time_step,dimensions,natoms,lx,ly,lz

def ReadCoordinates(f,Ncolloids):
	"""
	This function reads the coordinates of the positions of particles from a xml file 
	used for logging the positions of particles in a simulation.
	
	Input:	the filename (f), 
		the Number of Particles (Ncolloids)
	Output: The coordinates (coordinates)
	"""
	coordinates = np.zeros([Ncolloids,3])
	for i in range(Ncolloids):
		line = f.readline().strip().split(' ')
		try:
                    coordinates[i,0]=line[0]
                except ValueError:
                    break
		coordinates[i,1]=line[1]
		coordinates[i,2]=line[2]
	line = f.readline()
	return coordinates

def ReadMass(f,Ncolloids):
	"""
	This function reads the masses of the colloids.
	Input: 	directions to the file (f)
		the number of colloids (Ncolloids)
	Output:	-
	"""
	mass = np.zeros([Ncolloids])
	line = f.readline()
	for i in range (Ncolloids):
		line = f.readline().strip()
		mass[i] = line
	line = f.readline()
	
def ReadDiameter(f,Ncolloids):
	"""
	This function reads the diameters of the colloids.
	Input: 	directions to the file (f)
		the number of colloids (Ncolloids)
	Output:	-
	"""
	diameter = np.zeros(Ncolloids)
	line = f.readline()
	for i in range (Ncolloids):
		line = f.readline().strip()
		diameter[i] = line
	line = f.readline()

def ReadAtomTypeAB(f,Ncolloids):
	"""
	This function reads the atom types and provides the number of particles (A) and interstitials (B).
	
	Input:	data foldername in which the data of the simulation is saved (system)
	Output:	nparticles
		ninterstitials
	"""
	types = []
	line = f.readline()
	for i in range (Ncolloids):
		line = f.readline().strip()
		types.append(line)
	line = f.readline()	
	nparticles = types.count('A')
	ninterstitials = types.count('B')
	return nparticles, ninterstitials

def ReadInitialFile(system):
	"""
	Input:	data foldername in which the data of the simulation is saved (system)
	Output:	number of dimensions of the system (Dimensions),
		Number of colloids in the system (Ncolloids),
		Number of big particles in the system (Nparticles),
		Number of interstitials in the system (NInterstitials),
		Size of the simulation box in the x direction (lx), 
		Size of the simulation box in the y direction (ly),
		Size of the simulation box in the z direction (lz).	
	"""	
	InitialFilename = './DATA/' + system + '/InitialConfiguration.xml'
	InitialFile = open(InitialFilename, 'r')
	Time_Step,Dimensions,Ncolloids,Lx,Ly,Lz = ReadHeader(InitialFile)
	ReadCoordinates(InitialFile,Ncolloids)
	ReadMass(InitialFile,Ncolloids)
	ReadDiameter(InitialFile,Ncolloids)
	Nparticles,Ninterstitials = ReadAtomTypeAB(InitialFile,Ncolloids)
	InitialFile.close()
	return Dimensions,Ncolloids,Nparticles,Ninterstitials,Lx,Ly,Lz

"""
s.write('<type>\n')
for i in range(0,Nparticles):
	s.write('A\n')
for i in range(Nparticles,TotalParticles):
	s.write('B\n')
s.write('</type>\n')
"""

def writeheader(s,Lbox):
	header1 = '<?xml version="1.0" encoding="UTF-8"?>\n'
	header2 = '<hoomd_xml>\n'
	header3 = '<configuration time_step="0">\n'
	header4 = '<box lx="' + str(Lbox) + '" ly="' + str(Lbox) + '" lz="' + str(Lbox) + '" xy="0" xz="0" yz="0"/>\n'
	header5 = '<position>\n'
	s.write(header1 + header2 +header3 + header4 + header5)

def writepositions(s,Positions,Nparticles):	
	for i in range(0,Nparticles):
		s.write(str(Positions[i,0]) + ' ' + str(Positions[i,1]) + ' ' + str(Positions[i,2]) + '\n')
	s.write('</position>\n')

def writediameter(Input, s):	
	s.write('<diameter>\n')
	for i in range(Input['BASEparticles']):
		s.write('1\n')
	if Input['method'] == 'interstitial':	
		for i in range(Input['INTERparticles']):
			s.write(str(Input['INTERsigmacorrected']) + '\n')
	s.write('</diameter>\n')

def writetype(Input, s):
	s.write('<type>\n')
	for i in range(Input['BASEparticles']):
		s.write('A\n')
	if Input['method'] == 'interstitial':	
		for i in range(Input['INTERparticles']):
			s.write('B' + '\n')	
	s.write('</type>\n')

def writeend(s):
	s.write('</configuration>\n')
	s.write('</hoomd_xml>\n')

def mattoxml(Input, filename):
    """Converts mat files to xml files"""
    savename = '%s%s' % (Input['folder_mat_data'], filename)
    framename = filename.strip('.mat')
    conversion = np.array([0.25,0.25,0.33])
    shift = 0.5*np.array([Input['BOXsize_x'], Input['BOXsize_y'], Input['BOXsize_z']])
    Positions = (scipy.io.loadmat(savename)[framename]*conversion)/1.8
    Positions = Positions - shift
    framenumber = int(framename.strip('frame'))
    Filename = '%sWigner.%010d.xml' % (Input['folder_snapshot_data'], framenumber)
    s = open(Filename, "w")
    writeheader(s,Input['BOXsize'])
    writepositions(s,Positions,Input['ALLparticles'])
    writediameter(Input, s)
    writetype(Input, s)
    writeend(s)
    s.close()

def CubicLattice(Positions,count,a,Offset,xm,ym,zm,xshift=0,yshift=0,zshift=0):
	for i in range(xm):
		for j in range(ym):	
			for k in range(zm):
				x = round(-Offset[0] + i*a + xshift*a,8)
				y = round(-Offset[1] + j*a + yshift*a,8)
				z = round(-Offset[2] + k*a + zshift*a,8)
				Positions[count,0] = x
				Positions[count,1] = y
				Positions[count,2] = z			
				count = count + 1	
	return count

def InitFCC(FullDirectory,Input):
	"""
	This script generates an xml file, which can be used as an initial
	conditions file. It places particles of type A in a perfect FCC
	crystal configuration.
	The total amount of particles for a bcc crystal is m**3 + n**3 where n = m-1
	To find out the value for m the following equation has to be solved:
	4m**3 - 6m**2 +3m-Nparticles = 0. This function can be solved with np.roots
	Where the 4 numbers are the coefficients a,b,c,d in a cubic formula
	"""
	Filename = FullDirectory + 'InputConfiguration.xml'
	Offset = 0.5*0.95*Input['BOXsize'] # To get (x,y,z)=(0,0,0) as centre of the system
	a = 0.95*Input['BOXsize']/(Input['m']-1.0) # Distance between particles along axis

	Positions = np.zeros([Nparticles,3])
	count = 0
	count = CubicLattice(Positions,count,a,Offset,Input['m'],Input['m'],Input['m']) #Defining basis crystal
	count = CubicLattice(Positions,count,a,Offset,Input['m'],Input['n'],Input['n'],0,0.5,0.5) # Defining first crystal with offset
	count = CubicLattice(Positions,count,a,Offset,Input['n'],Input['m'],Input['n'],0.5,0,0.5) # Defining second crystal with offset
	count = CubicLattice(Positions,count,a,Offset,Input['n'],Input['n'],Input['m'],0.5,0.5,0) # Defining third crystal with offset	

	s = open(Filename, "w")
	writeheader(s,Input['BOXsize'])
	writepositions(s,Positions,Input['BASEparticles'])
	writediameter(s,Input['BASEparticles'])
	writetype(s,Input['BASEparticles'])
	writeend(s)
	s.close()

def InitBCC(Input, FullDirectory):
	"""
	This script generates an xml file, which can be used as an initial
	conditions file. It places particles of type A in a perfect BCC
	crystal configuration.
	The total amount of particles for a bcc crystal is m**3 + n**3 where n = m-1
	To find out the value for m the following equation has to be solved:
	2m**3 - 3m**2 +3m -1 -Nparticles. This function can be solved with np.roots
	Where the 4 numbers are the coefficients a,b,c,d in a cubic formula
	"""
	Filename = FullDirectory + 'InputConfiguration.xml'
	a = Input['BOXsize']/(Input['m']) # Distance between particles along axis
	Offset_x = 0.5*Input['BOXsize_x']-0.25*a # To get (x,y,z)=(0,0,0) as centre of the system 
	Offset_y = 0.5*Input['BOXsize_y']-0.25*a
	Offset_z = 0.5*Input['BOXsize_z']-0.25*a
	Offset = (Offset_x,Offset_y,Offset_z)
	Positions = np.zeros([Input['ALLparticles'],3])
	count = 0
	count = CubicLattice(Positions,count,a,Offset,Input['m'],Input['m'],Input['m'])	
	count = CubicLattice(Positions,count,a,Offset,Input['m'],Input['m'],Input['m'],0.5,0.5,0.5)

	if Input['method'] == 'interstitial':
		RandomInterstitials = np.zeros(Input['INTERsites']) #-6*(Input['m']**2)])
		RandomInterstitials[:Input['INTERparticles']] = 1
		np.random.shuffle(RandomInterstitials)
		Offsets = np.array([[0,0.5,0.25],[0,0.25,0.5],[0,0.5,0.75],[0,0.75,0.5]])
		InterstitialCount = 0
		for l in range(0,3): # walk over the planes xy,xz,yz
			MorN = np.zeros(3)
			MorN[l] = 1
			for o in range(0,4): # walking over the 4 interstitial particles per plain
				for i in range(0,Input['n'] + int(MorN[0])):
					for j in range(0,Input['n'] + int(MorN[2])):
						for k in range(0,Input['n'] + int(MorN[1])):
							xOffset = Offsets[o,(0+l)%3]
							yOffset = Offsets[o,(1+l)%3]
							zOffset = Offsets[o,(2+l)%3]
							x = round(-Offset[0] + i*a + xOffset*a,8)					
							y = round(-Offset[1] + j*a + yOffset*a,8)
							z = round(-Offset[2] + k*a + zOffset*a,8)	
							if RandomInterstitials[InterstitialCount] == 1:					
								Positions[count,0] = x
								Positions[count,1] = y
								Positions[count,2] = z			
								count = count + 1
							InterstitialCount = InterstitialCount + 1	
	s = open(Filename, "w")
	writeheader(s,Input['BOXsize'])
	writepositions(s,Positions,Input['ALLparticles'])
	writediameter(Input, s)
	writetype(Input, s)
	writeend(s)
	s.close()

def InitBCC_experimental(Input, FullDirectory):
	"""
	This script generates an xml file, which can be used as an initial
	conditions file. It places particles of type A in a perfect BCC
	crystal configuration in a box with fixed sizes.
	"""
	Filename = FullDirectory + 'InputConfiguration.xml'
	a = Input['BOXsize']/(Input['m']) # Distance between particles along axis
	Offset_x = 0.5*Input['BOXsize_x']-0.25*a # To get (x,y,z)=(0,0,0) as centre of the system 
	Offset_y = 0.5*Input['BOXsize_y']-0.25*a
	Offset_z = 0.5*Input['BOXsize_z']-0.25*a
	Offset = (Offset_x,Offset_y,Offset_z)
	Positions = np.zeros([Input['ALLparticles'],3])
	count = 0
	count = CubicLattice(Positions,count,a,Offset,Input['m'],Input['m'],Input['m']/2)	
	count = CubicLattice(Positions,count,a,Offset,Input['m'],Input['m'],Input['m']/2,xshift=0.5,yshift=0.5,zshift=0.5)
	
	if Input['method'] == 'interstitial':
		RandomInterstitials = np.zeros(Input['INTERsites']/3) #-6*(Input['m']**2)])
		RandomInterstitials[:Input['INTERparticles']] = 1
		np.random.shuffle(RandomInterstitials)
		Offsets = np.array([[0,0.5,0.25],
				    [0,0.25,0.5],
				    [0,0.5,0.75],
				    [0,0.75,0.5]])
		InterstitialCount = 0
		for l in range(0,3): # walk over the planes xy,xz,yz
			MorN = np.zeros(3)
#			MorN[l] = 1
			for o in range(0,4): # walking over 4 sites
				for i in range(0,Input['n'] + int(MorN[0])):
					for j in range(0,Input['n'] + int(MorN[2])):
						for k in range(0,Input['n'] + int(MorN[1]) - Input['m']/2):
							xOffset = Offsets[o,(0+l)%3]
							yOffset = Offsets[o,(1+l)%3]
							zOffset = Offsets[o,(2+l)%3]
							x = round(-Offset[0] + i*a + xOffset*a,8)
							y = round(-Offset[1] + j*a + yOffset*a,8)
							z = round(-Offset[2] + k*a + zOffset*a,8)	
							if RandomInterstitials[InterstitialCount] == 1:				
								Positions[count,0] = x
								Positions[count,1] = y
								Positions[count,2] = z			
								count = count + 1
							InterstitialCount += 1	
	s = open(Filename, "w")
	writeheader(s,Input['BOXsize'])
	writepositions(s,Positions,Input['ALLparticles'])
	writediameter(Input, s)
	writetype(Input, s)
	writeend(s)
	s.close()

