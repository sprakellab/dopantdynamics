"""This module contains functions for plotting energy data."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import colormaps as cmaps

def PlotOverview(System,Bin_Edges_Norm,TotalNormHist,SaveName,show=False):
	EnergyFilename = './DATA/' + System + '/energies.log'
	Data = np.loadtxt(EnergyFilename, comments='#')
	fig, ax = plt.subplots(3)
	ax[0].plot(Bin_Edges_Norm[1:(Bin_Edges_Norm.size-1)], TotalNormHist, color='k')
	ax[0].set_xlabel('r/a')
	ax[0].set_ylabel('G(r)')
	ax[0].set_xlim(0,10)
	ax[1].semilogx(Data[:,0], Data[:,2], color='k')
	ax[1].set_ylabel('$U_{pot}$')
	ax[1].xaxis.set_visible(False)
	ax[1].ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	ax[2].semilogx(Data[:,0], Data[:,3], color='k')
	ax[2].set_ylabel('$T$ [$k_B T$]')
	ax[2].set_xlabel('Steps')
	fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.2)
	pdfEnergySaveName = SaveName + 'Overview' + '.pdf'
	pngEnergySaveName = SaveName + 'Overview' + '.png'
	fig.savefig(pdfEnergySaveName)
	fig.savefig(pngEnergySaveName)
	if show is True:
		plt.show()

def plot_energy_track(time, energy):
    fig,ax = plt.subplots()
    ax.plot(time, energy)
    plt.show()
