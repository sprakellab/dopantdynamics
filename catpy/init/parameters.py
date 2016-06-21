"""This scriptguides the user in setting up an initialization file
for running simulations on interstitials."""

import numpy as np

string_length = 40
output_length = 30

def InitializeSystem(output):
    """This functions is used to ask the user how the system is build.
    Furthermore it assigns system parameters.
    """
    output['version'] = int(raw_input('%-*s' % (string_length,'version number: ')))
    output['method'] = raw_input('%-*s' % (string_length,'wigner or interstitial: '))
    output['basegenerator'] = raw_input('%-*s' % (string_length,'Nparticles or Lbox: '))
    output['basefix'] = raw_input('%-*s' % (string_length,'Fixed base True or False: ')) == 'True'
    output['series'] = raw_input('%-*s' % (string_length,'series(True or False): ')) == 'True'
    output['divide'] = raw_input('%-*s' % (string_length,'Divide the simulation in blocks?: ')) == 'True'
    output['InitialConfiguration'] = raw_input('%-*s' % (string_length,'random, FCC or BCC: '))
    if output['basegenerator'] == 'Lbox' or output['InitialConfiguration'] == 'random':
        output['BOXsize_x'] = float(raw_input('%-*s' % (string_length,'Length of box side x [sigma]: ')))
        output['BOXsize_y'] = float(raw_input('%-*s' % (string_length,'Length of box side y [sigma]: ')))
        output['BOXsize_z'] = float(raw_input('%-*s' % (string_length,'Length of box side z [sigma]: ')))
        output['BOXsize'] = output['BOXsize_x']
    else:
        output['BASEparticles'] = int(raw_input('%-*s' % (string_length,'Number of Base particles: ')))
    output['viscosity'] = 0.932e-3
    output['ReferenceT'] = 295.15 # adapated to 22 degrees
    output['Boltzmann'] = 1.3806488e-23
    output['temperature'] = 1.0
    output['rcut'] = 10.0
    return output

def InitializeParticles(output):
    """This functions is used to ask the user for input for parameters."""
    if output['series'] is True:
        output['BASEfractionStart'] = float(raw_input('%-*s' % (string_length,'Start of Base Volume Fraction series: ')))
        output['BASEfractionEnd'] = float(raw_input('%-*s' % (string_length,'End of Base Volume Fraction series: ')))
        output['BASEfractionStep'] = float(raw_input('%-*s' % (string_length,'Step size of Base Volume Fraction: ')))
    else:
        output['BASEfraction'] = float(raw_input('%-*s' % (string_length,'Base Volume Fraction: ')))
    output['BASEsigma'] = float(raw_input('%-*s' % (string_length,'Base diameter [m]: ')))
    output['BASEsigmacorrected'] = 1.0
    output['BASEepsilon'] = float(raw_input('%-*s' % (string_length,'Base Epsilon [kT*sigma]: ')))
    output['BASEkappa'] = float(raw_input('%-*s' % (string_length,'Base Kappa [1/sigma]: ')))
    if output['method'] == 'interstitial':
        output['INTERfraction'] = float(raw_input('%-*s' % (string_length,'Percentage of filled sites: ')))
        output['INTERsigma'] = float(raw_input('%-*s' % (string_length,'Interstitial diameter [m]: ')))
        output['INTERepsilon'] = float(raw_input('%-*s' % (string_length,'Interstitial Epsilon [kT*sigma]: ')))
        output['INTERkappa'] = float(raw_input('%-*s' % (string_length,'Interstitial Kappa [1/sigma]: ')))
    return output

def InitializeCalculations(output):
    """This functions calculates parameters from the input"""
    # SYSTEM
    if output['series'] is False:	
    	output['system'] = 'F' + "{0:.2f}".format(output['BASEfraction']).replace('.','_') + 'D1_8V' + \
                           '0'*(len(str(output['version']))%2) + str(output['version'])
        output['systemfolder'] = output['system']
    # BASE PARTICLES	
    output['BASEparticlevolume'] = (4.0/3.0) * np.pi * (output['BASEsigmacorrected']/2.0)**3	
    output['BASEdiffusion'] = (output['ReferenceT']*output['Boltzmann']) / (6.0*np.pi*output['viscosity']*(0.5*output['BASEsigma']))
    output['BASEtau'] = (output['BASEsigma']**2) / output['BASEdiffusion']	
    output['eta'] = output['viscosity']*(output['BASEsigma']**3)/(output['Boltzmann']*output['ReferenceT']*output['BASEtau'])
    output['BASEfriction'] = 0.5*6.0*output['BASEsigmacorrected']*np.pi*output['eta'];
    if output['series'] is False:
        if output['basegenerator'] == 'Lbox' or output['InitialConfiguration'] == 'random':
            output['BOXvolume'] = output['BOXsize_x']*output['BOXsize_y']*output['BOXsize_z']
            TotalParticleVolume = output['BOXvolume'] * output['BASEfraction']
            output['BASEparticles'] = int(np.round(TotalParticleVolume / output['BASEparticlevolume'],0))
        elif output['InitialConfiguration'] == 'BCC' or output['InitialConfiguration'] == 'FCC':	
            output['BOXsize'] = (output['BASEparticlevolume']*output['BASEparticles']/output['BASEfraction'])**(1/3.0)
            output['BOXsize_x'] = output['BOXsize']
            output['BOXsize_y'] = output['BOXsize']
            output['BOXsize_z'] = output['BOXsize']
            output['BOXvolume'] = output['BOXsize_x']*output['BOXsize_y']*output['BOXsize_z']			
    if output['InitialConfiguration'] == 'BCC' and output['basegenerator'] == 'Nparticles':
        roots = np.roots([2,0,0,-1.0*output['BASEparticles']])[::-1] # Nparticles of the basis crystal along one axis
        m = roots[0]
        output['m'] = int(round(m.real,0))
        output['n'] = output['m']-1 # number of unit cells
        output['a'] = output['BOXsize'] / output['m']
    elif output['InitialConfiguration'] == 'FCC': 
        roots = np.roots([4,-6,3,-output['BASEparticles']]) # number of particles of the basis crystal along one axis
        m = roots[0]
        output['m'] = int(round(m.real,0))
        output['n'] = output['m']-1 # number of unit cells
        output['a'] = output['BOXsize'] / output['m']
    elif output['InitialConfiguration'] == 'BCC' and output['basegenerator'] == 'Lbox':
        # for a box with ratio's x,y = 1 and z = 0.5
        output['m'] = int(np.round(np.power(output['BASEparticles'],(1/3.0)),0))
        output['n'] = output['m']-1 # number of unit cells
        output['a'] = output['BOXsize'] / output['m']
    if output['method'] == 'wigner':
        output['ALLparticles'] = output['BASEparticles'] 
    # INTERSTITIALS
    elif output['method'] == 'interstitial':		
        output['INTERsigmacorrected'] = output['INTERsigma'] / output['BASEsigma']	
        output['INTERfriction'] = 0.5*6.0*output['INTERsigmacorrected']*np.pi*output['eta']
        output['INTERdiffusion'] = (output['ReferenceT']*output['Boltzmann']) / (6.0*np.pi*output['viscosity']*(0.5*output['INTERsigma']))
        output['INTERtau'] = (output['INTERsigma']**2) / output['INTERdiffusion']
        output['IBepsilon'] = 0.5*(output['BASEepsilon'] + output['INTERepsilon'])
        output['IBkappa'] = 0.5*(output['BASEkappa'] + output['INTERkappa'])
        if output['InitialConfiguration'] == 'BCC':
            output['INTERsites'] = 12*output['m']*output['n']*output['n']
            output['INTERparticles'] = int(output['INTERsites']*output['INTERfraction'])
        output['ALLparticles'] = output['BASEparticles'] + output['INTERparticles']
    return output

def InitializeSimulation(output):
    """This functions is used to ask for simulation parameters."""
    output['Nsteps'] = int(float(raw_input('%-*s' % (string_length,'Number of steps: '))))
    output['dt'] = float(raw_input('%-*s' % (string_length,'step size: ')))
    if output['divide']:
        output['max_simulation_length'] = 100000000
        output['Nblocks'] = (output['Nsteps'] / output['max_simulation_length']) + 1
    return output	

def InitializeOutput(output):
    """This functions is used to ask the user for input for parameters."""	
    output['dcd'] = False
    output['XML'] = True
    output['XMLsaveperiod'] = int(raw_input('%-*s' % (string_length,'xml save period (100000): ')))
    output['XMLoutputperiod'] = int(raw_input('%-*s' % (string_length,'xml output period (500): ')))
    output['XMLinput'] = True
    return output

def InitializeAnalysis(output):
    """This functions is used to ask the user how to analyze."""
    output['analysisstart'] = 10000000
    output['stepseval'] = 20000
    output['indexeval'] = (output['stepseval'] / output['XMLoutputperiod']) / 2
    output['hth_standard'] = 0.036
    output['BASEtrack'] = raw_input('%-*s' % (string_length,'track Base particles (True or False):')) == 'True'
    if output['divide']:
        output['max_simulation_filename'] = 'Wigner.%010d.xml' % (output['max_simulation_length']-output['XMLsaveperiod'])
    output['experimental'] = raw_input('%-*s' % (string_length,'Subject data to experimental limitations (True or False):')) == 'True'
    output['max_range_IPD'] = 10.0
    # IPD analysis	
    output['first_step_IPD'] = int(((output['Nsteps'] / output['XMLsaveperiod']) / 2 )) * output['XMLsaveperiod']
    output['last_step_IPD'] = output['Nsteps'] - output['XMLsaveperiod']
    output['savetxtIPD'] = True
    output['savematIPD'] = False
    output['msd'] = False
    output['hthfrommsd'] = True
    return output

def InitializeFolders(output):
    """Here the names of the folders are introduced"""
    output['folder_data'] = 'DATA'
    output['systempath'] = './%s/%s/' % (output['folder_data'], output['systemfolder'])
    if output['experimental']:
        output['folder_interstitial_track'] = '%sboxtracks/' % (output['systempath'])
    else:
        output['folder_interstitial_track'] = '%sinterstitialsorted/' % (output['systempath'])
    output['snapshot_data_name'] = 'Wigner'	
    output['snapshot_data_extension'] = 'xml'
    output['folder_snapshot_data'] = '%sposition/' % (output['systempath'])
    output['folder_processed_data'] = '%sprocesseddata/' % (output['systempath'])
    output['folder_dhop_data'] = '%sdhop/' % (output['systempath'])
    output['folder_thop_data'] = '%sthop/' % (output['systempath'])
    output['folder_harray_data'] = '%sharray/' % (output['systempath'])
    output['folder_backup'] = '%sbackuptrajectory/' % (output['systempath'])
    output['folder_sitesize_data'] = '%ssitesize/' % (output['systempath'])
    output['folder_msd_data'] = '%smsd/' % (output['systempath'])
    if output['BASEtrack']:
        output['folder_base'] = '%sbase/' % (output['systempath'])
        output['folder_base_sorted'] = '%sbasesorted/' % (output['systempath'])
        output['folder_harray_base'] = '%sbaseharray/' % (output['systempath'])
        output['folder_backup_base'] = '%sbackupbasetrajectory/' % (output['systempath'])        
    output['log_analyzer_path'] = './ANALYZER/%s.log' % output['system']
    output['savenameIPD'] = '%sF%sL%s' % (output['folder_processed_data'],output['first_step_IPD'],output['last_step_IPD'])
    return output

def InitializePlot(output):
	"""This functions is used to ask how to plot the data."""
	return output
	
def FullInitialization(output):
    """This funtion callsother function to gather all the information for an initialization."""
    InitializeSystem(output)
    InitializeParticles(output)
    InitializeCalculations(output)
    InitializeSimulation(output)
    InitializeOutput(output)
    InitializeAnalysis(output)
    InitializeFolders(output)
    InitializePlot(output)
    return output

def ExportInitialization(f,output):
    """From the standard initilization dict a .init file is created."""
    f.write('# SYSTEM\n')
    f.write('%-*s%d\n' % (output_length,'version:',output['version']))	
    if output['series'] is False:
        f.write('%-*s%s\n' % (output_length,'system:',output['system']))
        f.write('%-*s%s\n' % (output_length,'systemfolder:',output['systemfolder']))	
    f.write('%-*s%s\n' % (output_length,'method:',output['method']))
    f.write('%-*s%s\n' % (output_length,'series:',output['series']))
    f.write('%-*s%s\n' % (output_length,'divide:',output['divide']))
    f.write('%-*s%s\n' % (output_length,'basegenerator:',output['basegenerator']))
    f.write('%-*s%s\n' % (output_length,'InitialConfiguration:',output['InitialConfiguration']))

    f.write('\n# SYSTEM PARAMETERS\n')
    f.write('%-*s%d\n' % (output_length,'BASEparticles:',output['BASEparticles']))
    if output['method'] == 'interstitial':
        f.write('%-*s%d\n' % (output_length,'INTERparticles:',output['INTERparticles']))
    f.write('%-*s%d\n' % (output_length,'ALLparticles:',output['ALLparticles']))
    if output['series'] is False:
        f.write('%-*s%.8e\n' % (output_length,'BOXsize:',output['BOXsize']))
        f.write('%-*s%.8e\n' % (output_length,'BOXsize_x:',output['BOXsize_x']))
        f.write('%-*s%.8e\n' % (output_length,'BOXsize_y:',output['BOXsize_y']))
        f.write('%-*s%.8e\n' % (output_length,'BOXsize_z:',output['BOXsize_z']))
        f.write('%-*s%.8e\n' % (output_length,'BOXvolume:',output['BOXvolume']))
    f.write('%-*s%.8e\n' % (output_length,'viscosity:',output['viscosity']))
    f.write('%-*s%.8e\n' % (output_length,'eta:',output['eta']))
    f.write('%-*s%.8e\n' % (output_length,'ReferenceT:',output['ReferenceT']))
    f.write('%-*s%.8e\n' % (output_length,'Boltzmann:',output['Boltzmann']))
    f.write('%-*s%.8e\n' % (output_length,'temperature:',output['temperature']))
    f.write('%-*s%.8e\n' % (output_length,'rcut:',output['rcut']))		
    if output['basegenerator'] == 'Nparticles':
        f.write('%-*s%d\n' % (output_length,'m:',output['m']))
        f.write('%-*s%d\n' % (output_length,'n:',output['n']))
        f.write('%-*s%.8e\n' % (output_length,'a:',output['a']))

    f.write('\n# BASE PARTICLES\n')
    if output['series'] is True:	
        f.write('%-*s%.8e\n' % (output_length,'BASEfractionStart:',output['BASEfractionStart']))	
        f.write('%-*s%.8e\n' % (output_length,'BASEfractionEnd:',output['BASEfractionEnd']))	
        f.write('%-*s%.8e\n' % (output_length,'BASEfractionStep:',output['BASEfractionStep']))		
    else:
        f.write('%-*s%.8e\n' % (output_length,'BASEfraction:',output['BASEfraction']))
    f.write('%-*s%.8e\n' % (output_length,'BASEsigma:',output['BASEsigma']))
    f.write('%-*s%.8e\n' % (output_length,'BASEsigmacorrected:',output['BASEsigmacorrected']))	
    f.write('%-*s%.8e\n' % (output_length,'BASEparticlevolume:',output['BASEparticlevolume']))	
    f.write('%-*s%.8e\n' % (output_length,'BASEepsilon:',output['BASEepsilon']))	
    f.write('%-*s%.8e\n' % (output_length,'BASEkappa:',output['BASEkappa']))	
    f.write('%-*s%.8e\n' % (output_length,'BASEtau:',output['BASEtau']))
    f.write('%-*s%.8e\n' % (output_length,'BASEdiffusion:',output['BASEdiffusion']))
    f.write('%-*s%.8e\n' % (output_length,'BASEfriction:',output['BASEfriction']))
    if output['method'] == 'interstitial':	
        f.write('\n# INTERSTITIAL PARTICLES\n')
        f.write('%-*s%.8e\n' % (output_length,'INTERfraction:',output['INTERfraction']))	
        f.write('%-*s%.8e\n' % (output_length,'INTERsigma:',output['INTERsigma']))
        f.write('%-*s%.8e\n' % (output_length,'INTERsigmacorrected:',output['INTERsigmacorrected']))	
        f.write('%-*s%.8e\n' % (output_length,'INTERepsilon:',output['INTERepsilon']))	
        f.write('%-*s%.8e\n' % (output_length,'INTERkappa:',output['INTERkappa']))
        f.write('%-*s%.8e\n' % (output_length,'INTERtau:',output['INTERtau']))
        f.write('%-*s%.8e\n' % (output_length,'INTERdiffusion:',output['INTERdiffusion']))
        f.write('%-*s%.8e\n' % (output_length,'INTERfriction:',output['INTERfriction']))
        f.write('%-*s%.8e\n' % (output_length,'INTERsites:',output['INTERsites']))

        f.write('\n# PARTICLE INTERACTION\n')
        f.write('%-*s%.8e\n' % (output_length,'IBepsilon:',output['IBepsilon']))	
        f.write('%-*s%.8e\n' % (output_length,'IBkappa:',output['IBkappa']))		

    f.write('\n# SIMULATION\n')
    f.write('%-*s%d\n' % (output_length,'Nsteps:',output['Nsteps']))	
    f.write('%-*s%.8e\n' % (output_length,'dt:',output['dt']))
    f.write('%-*s%s\n' % (output_length,'basefix:',output['basefix']))
    if output['divide']:
        f.write('%-*s%d\n' % (output_length,'max_simulation_length:',output['max_simulation_length']))
        f.write('%-*s%s\n' % (output_length,'max_simulation_filename:',output['max_simulation_filename']))
        f.write('%-*s%d\n' % (output_length,'Nblocks:',output['Nblocks']))

    f.write('\n# OUTPUT\n')
    f.write('%-*s%s\n' % (output_length,'dcd:',output['dcd']))
    f.write('%-*s%s\n' % (output_length,'XMLinput:',output['XMLinput']))	
    f.write('%-*s%s\n' % (output_length,'XML:',output['XML']))
    f.write('%-*s%d\n' % (output_length,'XMLoutputperiod:',output['XMLoutputperiod']))
    f.write('%-*s%d\n' % (output_length,'XMLsaveperiod:',output['XMLsaveperiod']))

    f.write('\n# ANALYSIS\n')
    f.write('%-*s%d\n' % (output_length,'analysisstart:',output['analysisstart']))
    f.write('%-*s%d\n' % (output_length,'stepseval:',output['stepseval']))
    f.write('%-*s%d\n' % (output_length,'indexeval:',output['indexeval']))
    f.write('%-*s%.8e\n' % (output_length,'hth_standard:',output['hth_standard']))
    f.write('%-*s%s\n' % (output_length,'BASEtrack:',output['BASEtrack']))
    f.write('%-*s%s\n' % (output_length,'experimental:',output['experimental']))
    f.write('%-*s%s\n' % (output_length,'max_range_IPD:',output['max_range_IPD']))
    f.write('%-*s%s\n' % (output_length,'first_step_IPD:',output['first_step_IPD']))
    f.write('%-*s%s\n' % (output_length,'last_step_IPD:',output['last_step_IPD']))
    f.write('%-*s%s\n' % (output_length,'savetxtIPD:',output['savetxtIPD']))
    f.write('%-*s%s\n' % (output_length,'savematIPD:',output['savematIPD']))
    f.write('%-*s%s\n' % (output_length,'msd:',output['msd']))
    f.write('%-*s%s\n' % (output_length,'hthfrommsd:',output['hthfrommsd']))

    f.write('\n# FOLDERS\n')
    f.write('%-*s%s\n' % (output_length,'folder_data:',output['folder_data']))
    f.write('%-*s%s\n' % (output_length,'systempath:',output['systempath']))
    f.write('%-*s%s\n' % (output_length,'folder_interstitial_track:',output['folder_interstitial_track']))
    f.write('%-*s%s\n' % (output_length,'folder_snapshot_data:',output['folder_snapshot_data']))
    f.write('%-*s%s\n' % (output_length,'folder_processed_data:',output['folder_processed_data']))
    f.write('%-*s%s\n' % (output_length,'folder_dhop_data:',output['folder_dhop_data']))
    f.write('%-*s%s\n' % (output_length,'folder_thop_data:',output['folder_thop_data']))
    f.write('%-*s%s\n' % (output_length,'folder_sitesize_data:',output['folder_sitesize_data']))
    f.write('%-*s%s\n' % (output_length,'folder_harray_data:',output['folder_harray_data']))
    f.write('%-*s%s\n' % (output_length,'folder_backup:',output['folder_backup']))
    f.write('%-*s%s\n' % (output_length,'folder_msd_data:',output['folder_msd_data']))
    f.write('%-*s%s\n' % (output_length,'log_analyzer_path:',output['log_analyzer_path']))
    f.write('%-*s%s\n' % (output_length,'snapshot_data_extension:',output['snapshot_data_extension']))
    f.write('%-*s%s\n' % (output_length,'snapshot_data_name:',output['snapshot_data_name']))
    f.write('%-*s%s\n' % (output_length,'savenameIPD:',output['savenameIPD']))

    f.write('\n# PLOTTING\n')
    f.write('\nEND')

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

def ImportInitialization(f):
    dictionary = {}
    line = f.readline()
    while line != 'END':
        line = f.readline().strip()
        if line != '' and line != 'END' and line[0] != '#':
            data = line.split()
            key = data[0].strip(':')
            if data[1] == 'True' or data[1] == 'False':
                Value = (data[1] == 'True')
            elif data[1][0].isalpha():
                Value = data[1] 
            elif data[1][0] == '.':
                Value = data[1]
            else:
                Value = num(data[1])
            dictionary[key] = Value
    return dictionary

if __name__ == "__main__":
    print('Please answer the following questions to construct an initialization file.\n' + 	\
    'These files can also be made manually and be used directly in the main script' + 	\
    'There is also an option to initialize directly from the main script.')
	
    Init = {}
    FullInitialization(Init)
    filepath = 'test.init'
    F = open(filepath, 'w')
    ExportInitialization(F,Init)
    F.close()
