"""Support to use voro++ from python"""

import subprocess
import os

import numpy as np

import catpy.basisfunctions as bf
import catpy.IO.xml as xml

def prepareforvoro(Input):
    voronoi_folder = '%spositions_voronoi/' % Input['systempath']
    os.mkdir(voronoi_folder)
    filelist = bf.lookupfiles(path=Input['folder_snapshot_data'], filechr='Wigner')
    for file_ in filelist:
        filename = '%s%s' % (Input['folder_snapshot_data'], file_)
        f = open(filename, 'r')
        xml.ReadHeader(f)
        Positions = xml.ReadCoordinates(f, Input['ALLparticles'])
        f.close()
        voronoiformat = np.concatenate((np.reshape(np.arange(1,Input['ALLparticles']+1),(-1,1)), Positions),axis=1)
        savename = '%s%s' % (voronoi_folder, file_.strip('.xml'))
        format_ = ['%6d','%.18e', '%.18e', '%.18e']
        np.savetxt(savename, voronoiformat, fmt=format_, delimiter='\t')

def find_cell_volumes_folder(Input):
    """find cell volumes of a whole folder"""
    positions_folder = '%spositions_voronoi/' % Input['systempath']
    cell_volume_folder = '%svolumes_voronoi/' % Input['systempath']
    os.mkdir(cell_volume_folder) 
    filelist = bf.lookupfiles(path=positions_folder, filechr='Wigner')
    for file_ in filelist:
        filename = positions_folder + file_
        new_name = positions_folder + file_ + '.vol'
        savename = cell_volume_folder + file_ + '.vol'
        find_cell_volume(Input, filename)     
        os.rename(new_name, savename)

def find_cell_volume(Input, filename):
    """find cell volumes"""
    call_voro(Input, filename)        

def call_voro(Input, filename):
    """Call the voro++ program"""
    x_min = str(-0.5*Input['BOXsize_x'])
    x_max = str(0.5*Input['BOXsize_x'])
    y_min = str(-0.5*Input['BOXsize_y'])
    y_max = str(0.5*Input['BOXsize_y'])
    z_min = str(-0.5*Input['BOXsize_z'])
    z_max = str(0.5*Input['BOXsize_z'])
    custom = '-c'
    output = '%i %x %y %z %v %s %c'
    periodic = '-p'
    args = ['voro++', periodic, custom, output, x_min, x_max, y_min, y_max, z_min, z_max, filename]
    p = subprocess.Popen(args)
    p.wait()

def average_values(Input):
    """get the average of the cell volume over the simulation"""
    cell_volume_folder = '%svolumes_voronoi/' % Input['systempath']
    filelist = bf.lookupfiles(path=cell_volume_folder, filechr='Wigner')
    data = np.genfromtxt(cell_volume_folder + filelist[0], delimiter=' ')
    average = np.zeros([len(data),2])
    for file_ in filelist:
        filename = '%s%s' % (cell_volume_folder,file_)
        data = np.genfromtxt(filename, delimiter=' ')
        data = data[data[:,0].argsort()]
        average[:,1] += data[:,4]*(1.0/len(filelist))
    average[:,0] = data[:,0]
    savename = '%s%saveragevoronoi.vol' % (Input['folder_processed_data'], Input['system'])
    format_ = ['%6d','%.12e']
    np.savetxt(savename, average, delimiter='\t', fmt=format_)        

def faces_values(Input):
    """get the faces of the cell for a snapshot"""
    cell_volume_file = '%svolumes_voronoi/Wigner.0099000000.vol' % Input['systempath']
    data = np.genfromtxt(cell_volume_file, delimiter=' ')
    data = data[data[:,0].argsort()]
    faces = np.zeros([len(data),2])
    faces[:,0] = data[:,0]
    faces[:,1] = data[:,5]
    savename = '%s%svoronoi.faces' % (Input['folder_processed_data'], Input['system'])
    format_ = ['%6d','%6d']
    np.savetxt(savename, faces, delimiter='\t', fmt=format_)   

def single_values(Input):
    """get the volumes of the cell for a snapshot"""
    cell_volume_file = '%svolumes_voronoi/Wigner.0099000000.vol' % Input['systempath']
    data = np.genfromtxt(cell_volume_file, delimiter=' ')
    data = data[data[:,0].argsort()]
    values = np.zeros([len(data),2])
    values[:,0] = data[:,0]
    values[:,1] = data[:,4]
    savename = '%s%ssinglevoronoi.vol' % (Input['folder_processed_data'], Input['system'])
    format_ = ['%6d','%.12e']
    np.savetxt(savename, values, delimiter='\t', fmt=format_)   

def voronoi_volume_histogram(Input, bins_=80):
    """Make histogram giving voronoi volume"""
    v_filename = '%s%saveragevoronoi.vol' % (Input['folder_processed_data'], Input['system'])
    data = np.genfromtxt(v_filename, delimiter='\t')
    hist1,binedges = np.histogram(data[Input['BASEparticles']:,1],bins=bins_, range=(3,9))
    combinedhistograms = np.hstack((np.reshape(binedges[:-1],(-1,1)),
                                    np.reshape(hist1/float(np.sum(hist1)),(-1,1)),
                                  ))
    savename = '%s%svoronoi.hist' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, combinedhistograms, delimiter='\t')                         

def voronoi_faces_histogram(Input, bins_=12):
    """Make histogram giving voronoi volume"""
    v_filename = '%s%svoronoi.faces' % (Input['folder_processed_data'], Input['system'])
    data = np.genfromtxt(v_filename, delimiter='\t')
    hist1,binedges = np.histogram(data[Input['BASEparticles']:,1],bins=bins_, range=(6,18))
    combinedhistograms = np.hstack((np.reshape(binedges[:-1],(-1,1)),
                                    np.reshape(hist1/float(np.sum(hist1)),(-1,1)),
                                  ))
    savename = '%s%sfacesvoronoi.hist' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, combinedhistograms, delimiter='\t')

def voronoi_single_volume_histogram(Input, bins_=80):
    """Make histogram giving voronoi volume"""
    v_filename = '%s%ssinglevoronoi.vol' % (Input['folder_processed_data'], Input['system'])
    data = np.genfromtxt(v_filename, delimiter='\t')
    hist1,binedges = np.histogram(data[Input['BASEparticles']:,1],bins=bins_, range=(3,9))
    combinedhistograms = np.hstack((np.reshape(binedges[:-1],(-1,1)),
                                    np.reshape(hist1/float(np.sum(hist1)),(-1,1))
                                  ))
    savename = '%s%ssinglevolumevoronoi.hist' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, combinedhistograms, delimiter='\t')

def voronoi_volume_msdhistogram(Input, bins_=40,cut1=15, cut2=30):
    """Make histogram giving voronoi volume"""
    v_filename = '%s%saveragevoronoi.vol' % (Input['folder_processed_data'], Input['system'])
    m_filename = '%s%smsdparticles.val' % (Input['folder_processed_data'], Input['system'])
    voronoidata = np.genfromtxt(v_filename, delimiter='\t')
    msddata = np.genfromtxt(m_filename, delimiter='\t')
    mask1 = np.nonzero(np.less(msddata[:,1],cut1))
    mask2 = np.nonzero(np.greater(msddata[:,1],cut1)*np.less(msddata[:,1],cut2))
    mask3 = np.nonzero(np.greater(msddata[:,1],cut2))
    data1 = np.hstack((msddata[mask1],np.reshape(voronoidata[mask1,1],(-1,1))))
    data2 = np.hstack((msddata[mask2],np.reshape(voronoidata[mask2,1],(-1,1))))
    data3 = np.hstack((msddata[mask3],np.reshape(voronoidata[mask3,1],(-1,1))))
    hist1,binedges = np.histogram(data1[:,2],bins=bins_, range=(6,9))
    hist2,binedges = np.histogram(data2[:,2],bins=bins_, range=(6,9))
    hist3,binedges = np.histogram(data3[:,2],bins=bins_, range=(6,9))
    combinedhistograms = np.hstack((np.reshape(binedges[:-1],(-1,1)),
                                    np.reshape(hist1/float(np.sum(hist1)),(-1,1)),
                                    np.reshape(hist2/float(np.sum(hist2)),(-1,1)),
                                    np.reshape(hist3/float(np.sum(hist3)),(-1,1))
                                  ))
    savename = '%s%svoronoimsd.hist' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, combinedhistograms, delimiter='\t')                         

def voronoi_faces_msdhistogram(Input, bins_=6,cut1=0.75, cut2=1.50):
    """Make histogram giving voronoi volume"""
    v_filename = '%s%svoronoi.faces' % (Input['folder_processed_data'], Input['system'])
    m_filename = '%s%smsdparticles.val' % (Input['folder_processed_data'], Input['system'])
    voronoidata = np.genfromtxt(v_filename, delimiter='\t')
    msddata = np.genfromtxt(m_filename, delimiter='\t')
    mask1 = np.nonzero(np.less(msddata[:,1],cut1))
    mask2 = np.nonzero(np.greater(msddata[:,1],cut1)*np.less(msddata[:,1],cut2))
    mask3 = np.nonzero(np.greater(msddata[:,1],cut2))
    data1 = np.hstack((msddata[mask1],np.reshape(voronoidata[mask1,1],(-1,1))))
    data2 = np.hstack((msddata[mask2],np.reshape(voronoidata[mask2,1],(-1,1))))
    data3 = np.hstack((msddata[mask3],np.reshape(voronoidata[mask3,1],(-1,1))))
    hist1,binedges = np.histogram(data1[:,2],bins=bins_, range=(12,18))
    hist2,binedges = np.histogram(data2[:,2],bins=bins_, range=(12,18))
    hist3,binedges = np.histogram(data3[:,2],bins=bins_, range=(12,18))
    combinedhistograms = np.hstack((np.reshape(binedges[:-1],(-1,1)),
                                    np.reshape(hist1/float(np.sum(hist1)),(-1,1)),
                                    np.reshape(hist2/float(np.sum(hist2)),(-1,1)),
                                    np.reshape(hist3/float(np.sum(hist3)),(-1,1))
                                  ))
    savename = '%s%sfacesvoronoi.hist' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, combinedhistograms, delimiter='\t')

def voronoi_single_volume_msdhistogram(Input, bins_=40,cut1=0.75, cut2=1.50):
    """Make histogram giving voronoi volume"""
    v_filename = '%s%ssinglevoronoi.vol' % (Input['folder_processed_data'], Input['system'])
    m_filename = '%s%smsdparticles.val' % (Input['folder_processed_data'], Input['system'])
    voronoidata = np.genfromtxt(v_filename, delimiter='\t')
    msddata = np.genfromtxt(m_filename, delimiter='\t')
    mask1 = np.nonzero(np.less(msddata[:,1],cut1))
    mask2 = np.nonzero(np.greater(msddata[:,1],cut1)*np.less(msddata[:,1],cut2))
    mask3 = np.nonzero(np.greater(msddata[:,1],cut2))
    data1 = np.hstack((msddata[mask1],np.reshape(voronoidata[mask1,1],(-1,1))))
    data2 = np.hstack((msddata[mask2],np.reshape(voronoidata[mask2,1],(-1,1))))
    data3 = np.hstack((msddata[mask3],np.reshape(voronoidata[mask3,1],(-1,1))))
    hist1,binedges = np.histogram(data1[:,2],bins=bins_, range=(6,9))
    hist2,binedges = np.histogram(data2[:,2],bins=bins_, range=(6,9))
    hist3,binedges = np.histogram(data3[:,2],bins=bins_, range=(6,9))
    combinedhistograms = np.hstack((np.reshape(binedges[:-1],(-1,1)),
                                    np.reshape(hist1/float(np.sum(hist1)),(-1,1)),
                                    np.reshape(hist2/float(np.sum(hist2)),(-1,1)),
                                    np.reshape(hist3/float(np.sum(hist3)),(-1,1))
                                  ))
    savename = '%s%ssinglevolumevoronoi.hist' % (Input['folder_processed_data'], Input['system'])
    np.savetxt(savename, combinedhistograms, delimiter='\t')
