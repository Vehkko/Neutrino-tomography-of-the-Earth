#==============================================================================
# This script plots the energy and cos(th) distributions for various systematic
#     sets. We also show how the change between a systematic compared to the
#     nominal distribution.
#
#==============================================================================

# The systematic arrays below describe the detector response. We convolute the
#     array with the bin averaged flux to determine the number of events in a
#     given energy/cos(th) bin.

import numpy as np
import tables
import pylab
import sys

font = {'weight': 200,
        'size':   22}
pylab.rc('font', **font)

# Set variable for the nominal detector response arrays:
nomArray = '../../systematics_response_arrays/nufsgen_nominal_response_array.h5'

# Set variables for the detector response arrays:
sysArray1 = '../../systematics_response_arrays/nufsgen_dom_eff_1_1979_response_array.h5'
sysArray2 = '../../systematics_response_arrays/nufsgen_dom_eff_1_089_response_array.h5'
sysArray3 = '../../systematics_response_arrays/nufsgen_dom_eff_0_95_response_array.h5'
sysArray4 = '../../systematics_response_arrays/nufsgen_dom_eff_0_90_response_array.h5'
sysArray5 = '../../systematics_response_arrays/nufsgen_spicemie_icevariant2_response_array.h5'
sysArray6 = '../../systematics_response_arrays/nufsgen_spicemie_icevariant1_response_array.h5'
sysArray7 = '../../systematics_response_arrays/nufsgen_noholeice_response_array.h5'
sysArray8 = '../../systematics_response_arrays/nufsgen_spicelea_response_array.h5'

# Set up variables for the flux models
fluxArray1 = '../../atmospheric_flux/averaged/CombinedGHandHG_H3a_QGSJET-II-04.h5'
fluxArray2 = '../../atmospheric_flux/averaged/CombinedGHandHG_H3a_SIBYLL2.3_rc1_pl.h5'
fluxArray3 = '../../atmospheric_flux/averaged/HondaGaisser.h5'
fluxArray4 = '../../atmospheric_flux/averaged/PolyGonato_QGSJET-II-04.h5'
fluxArray5 = '../../atmospheric_flux/averaged/PolyGonato_SIBYLL2.3_rc1_pl.h5'
fluxArray6 = '../../atmospheric_flux/averaged/ZatsepinSokolskaya_pamela_QGSJET-II-04.h5'
fluxArray7 = '../../atmospheric_flux/averaged/ZatsepinSokolskaya_pamela_SIBYLL2.3_rc1_pl.h5'

# Create a class to import the data from the detector response arrays.
class sysData():
    def __init__(self,infile):
        h5_f=tables.open_file(infile,'r')
        self.name = infile.split('/')[-1]
        self.antineutrino_response_array = h5_f.root.antineutrino_response_array.read()
        self.neutrino_response_array = h5_f.root.neutrino_response_array.read()
        self.costh_bin_edges = h5_f.root.costh_bin_edges.read()
        self.proxy_energy_bin_edges = h5_f.root.proxy_energy_bin_edges.read()
        self.true_energy_bin_edges = h5_f.root.true_energy_bin_edges.read()
        h5_f.close()
# Create a class to read the data from the flux files.
class atmData():
    def __init__(self,infile):
        h5_f=tables.open_file(infile,'r')
        self.name = infile.split('/')[-1]
        self.average_flux_nu_kaon = h5_f.root.average_flux_nu_kaon.read()
        self.average_flux_nu_pion = h5_f.root.average_flux_nu_pion.read()
        self.average_flux_nubar_kaon = h5_f.root.average_flux_nubar_kaon.read()
        self.average_flux_nubar_pion = h5_f.root.average_flux_nubar_pion.read()
        self.costh_reco_bin_edges = h5_f.root.costh_reco_bin_edges.read()
        self.e_true_bin_edges = h5_f.root.e_true_bin_edges.read()
        h5_f.close()

# These variables let you choose what you want to compare. For example, if you
# want to look at the number of events in a bin using a 1.1979 DOM
# efficiency, use sysArray1. We also compare the deistribution to the
# nominal detector response function called nomArray.
systematicArray = sysData(sysArray2)
nominalArray = sysData(nomArray)
# The fluxArray variable lets you select which flux model you want to usechoose
# the flux model you are interested in:
fluxArray = atmData(fluxArray4)

# This function convolves the response array and the flux array to give the
# number of events in a given mu_e and cth bin
def Convolve(sys_array,flux_array):
    flux = np.zeros((21,10))
    for i in range(21):
        energyArray = np.zeros(10)
        for j in range(200):
            energyArray += np.multiply(sys_array[:,i,j],flux_array[i,j])
        flux[i] = energyArray
    return flux

# This class is used to sum up the kaon and pion neutrino fluxes, and output
# arrays we want to plot
class ConvolvedMC():
    def __init__(self,sys,flux):
    	kaon_nubar = Convolve(sys.antineutrino_response_array,flux.average_flux_nubar_kaon)
      	kaon_nu    = Convolve(sys.neutrino_response_array,flux.average_flux_nu_kaon)
      	pion_nubar = Convolve(sys.antineutrino_response_array,flux.average_flux_nubar_pion)
      	pion_nu    = Convolve(sys.neutrino_response_array,flux.average_flux_nu_pion)
      	total_bin_count = kaon_nubar + kaon_nu + pion_nu + pion_nubar
        self.counts = total_bin_count.transpose()
        self.cosZ = np.asarray(sum(total_bin_count.transpose()[i,:] for i in range(10)))
        self.energy = np.asarray(sum(total_bin_count.transpose()[:,j] for j in range(21)))
        self.costh_bin_edges = sys.costh_bin_edges
        self.proxy_energy_bin_edges = sys.proxy_energy_bin_edges
        return
#==============================================================================

# This is just a function to plot the 2D reconstructed energy versus zenith distribution.
def plot2D(output_filename,nominal_counts,counts,x_array,y_array):
        fig= pylab.figure(figsize=(15,5))
        ax1 = fig.add_subplot(1,2,1)
        ymax = np.log10(2e4)
        ymin = np.log10(4e2)
        xmax = 0.24
        xmin = -1.
        y_array = np.log10(y_array)
        y_mid = (y_array[1:] + y_array[:-1])/2
        x_mid = (x_array[1:] + x_array[:-1])/2
        yticks=[3.,4.]
        ylabels=[r"$10^{3}$",r"$10^{4}$"]
        for i in range(len(x_mid)):
            for j in range(len(y_mid)):
                pylab.annotate("{:10.0f}".format(counts[j][i]),
                               xy=(x_mid[i],y_mid[j]),
                               size=8, color='black', rotation='vertical',
                               fontweight='bold', alpha = 0.5)
        im = pylab.imshow(counts, interpolation='nearest', origin='low',
                          extent=[x_array[0], x_array[-1], y_array[0], y_array[-1]],
                          cmap='Blues', aspect='auto')
        print("Total number of events: "+str(int(np.sum(np.sum(counts)))))
        cb2 = pylab.colorbar()
        cb2.set_label(r'Events')
        pylab.xlabel(r'$\cos({\theta}_{z}^{reco}$)')
        pylab.ylabel(r'${\rm E}^{proxy}$ [GeV]')
        pylab.axis([xmin, xmax, ymin, ymax])
        pylab.yticks(yticks,ylabels)

        ax1 = fig.add_subplot(1,2,2)
        counts = np.true_divide((counts-nominal_counts),nominal_counts)*100
        im = pylab.imshow(counts, interpolation='nearest', origin='low',
                          extent=[x_array[0], x_array[-1], y_array[0], y_array[-1]],
                          cmap='bwr_r', aspect='auto')
        cb2 = pylab.colorbar()
        cb2.set_label(r'% shift')
        pylab.xlabel(r'$\cos(\theta_z^{reco}$)')
        pylab.ylabel(r'${\rm E}^{proxy}$ [GeV]')
        pylab.axis([xmin, xmax, ymin, ymax])
        pylab.yticks(yticks,ylabels)
        pylab.clim(-50,50)
        #pylab.gcf().subplots_adjust(bottom=0.15)
        pylab.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        pylab.savefig(output_filename, bbox_inches='tight')
        pylab.clf()
        return

# Now we have to convolve the detector response array for the systematic with
#   the desired flux model, and the nominal detector response function with
#   the flux model.
convolvedSet = ConvolvedMC(systematicArray,fluxArray)
convolvedNom = ConvolvedMC(nominalArray,fluxArray)

# Plot the number of events in each bin given the convolved data.
sysDistribution = plot2D("2D_energy_vs_zenith_systematics.png",
                         convolvedNom.counts,convolvedSet.counts,
                         convolvedSet.costh_bin_edges,
                         convolvedSet.proxy_energy_bin_edges)
