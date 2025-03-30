
#==============================================================================
# This script is used to plot the energy distribution for the various
# 	systematic variants.
#==============================================================================
#
# The systematic arrays below describe the detector response. We convolute the
# array with the bin averaged flux to determine the number of events in a
# given energy/cos(th) bin.

import sys
import numpy as np
import pylab
import tables

# Plotting settings
font = {'weight':200,
        'size':22}
pylab.rc('font',**font)

# Systematic response arrays path
pathSystResponseArrays = '../../systematics_response_arrays/'
# Neutrino flux path
pathAveragedNeutrinoFlux = '../../atmospheric_flux/averaged/'

# Set variable for the nominal detector response arrays:
NominalDetectorConfiguration  = ('nufsgen_nominal',
                                 'Nominal')

# Set variables for the detector response arrays:
SystematicDetectorConfigurations = [
        ('nufsgen_dom_eff_1_1979','DOM eff: 1.1979'),
        ('nufsgen_dom_eff_1_089','DOM eff: 1.089'),
        ('nufsgen_dom_eff_0_95','DOM eff: 0.95'),
        ('nufsgen_dom_eff_0_90','DOM eff: 0.90'),
        ('nufsgen_spicemie_icevariant2','IceVariant 2'),
        ('nufsgen_spicemie_icevariant1','IceVariant 1'),
        ('nufsgen_noholeice','No Holes Ice'),
        ('nufsgen_spicelea','SpiceLea Ice Model')
]

# Set up variables for the flux models
AtmosphericNeutrinoModels = [
    ('CombinedGHandHG_H3a_QGSJET-II-04','Combined GH+HGp QGSJET-II-04'),
    ('CombinedGHandHG_H3a_SIBYLL2.3_rc1_pl','Combined GH+HGp SIBYLL2.3-RC1PL'),
    ('HondaGaisser','HondaGaisser'),
    ('PolyGonato_QGSJET-II-04','PolyGonato QGSJET-II-04'),
    ('PolyGonato_SIBYLL2.3_rc1_pl','PolyGonato SIBYLL2.3-RC1PL'),
    ('ZatsepinSokolskaya_pamela_QGSJET-II-04','Zatsepin-Sokolskaya/PAMELA QGSJET-II-04'),
    ('ZatsepinSokolskaya_pamela_SIBYLL2.3_rc1_pl','Zatsepin-Sokolskaya/PAMELA SIBYLL2.3-RC1PL')
]

# Create a class to import the data from the detector response arrays.
class sysData():
    def __init__(self,infile):
        h5_f=tables.open_file(pathSystResponseArrays+infile[0]
                              +'_response_array.h5','r')
        self.name = infile
        self.antineutrino_response_array = h5_f.root.\
                                           antineutrino_response_array.read()
        self.neutrino_response_array = h5_f.root.neutrino_response_array.read()
        self.costh_bin_edges = h5_f.root.costh_bin_edges.read()
        self.proxy_energy_bin_edges = h5_f.root.proxy_energy_bin_edges.read()
        self.true_energy_bin_edges = h5_f.root.true_energy_bin_edges.read()
        h5_f.close()

# Create a class to read the data from the flux files.
class atmData():
    def __init__(self,infile):
        h5_f=tables.open_file(pathAveragedNeutrinoFlux+infile[0]
                              +'.h5','r')
        self.name = infile
        self.average_flux_nu_kaon = h5_f.root.average_flux_nu_kaon.read()
        self.average_flux_nu_pion = h5_f.root.average_flux_nu_pion.read()
        self.average_flux_nubar_kaon = h5_f.root.average_flux_nubar_kaon.read()
        self.average_flux_nubar_pion = h5_f.root.average_flux_nubar_pion.read()
        self.costh_reco_bin_edges = h5_f.root.costh_reco_bin_edges.read()
        self.e_true_bin_edges = h5_f.root.e_true_bin_edges.read()
        h5_f.close()

# This function convolves the response array and the flux array to give the
# number of events in a given reconstructed muon energy and cos(th) bin.

def Convolve(sys_array,flux_array):
    flux = np.zeros((21,10))
    for i in range(21):
        energyArray = np.zeros(10)
        for j in range(200):
            energyArray += np.multiply(sys_array[:,i,j],flux_array[i,j])
        flux[i] = energyArray
    return flux

# This class is used to sum up the kaon and pion neutrino events expectations
# and output arrays we want to plot.

class ConvolvedMC():
    def __init__(self,sys,flux):
    	kaon_nubar = Convolve(sys.antineutrino_response_array,
                              flux.average_flux_nubar_kaon)
      	kaon_nu    = Convolve(sys.neutrino_response_array,
                              flux.average_flux_nu_kaon)
      	pion_nubar = Convolve(sys.antineutrino_response_array,
                              flux.average_flux_nubar_pion)
      	pion_nu    = Convolve(sys.neutrino_response_array,
                              flux.average_flux_nu_pion)
      	total_bin_count = kaon_nubar + kaon_nu + pion_nu + pion_nubar

        self.counts = total_bin_count.transpose()
        self.cosZ = np.asarray(sum(total_bin_count.transpose()[i,:]
                                   for i in range(10)))
        self.energy = np.asarray(sum(total_bin_count.transpose()[:,j]
                                     for j in range(21)))
        self.costh_bin_edges = sys.costh_bin_edges
        self.proxy_energy_bin_edges = sys.proxy_energy_bin_edges

        return None

#==============================================================================

# Get the energy from the convoled systematic variant with the Honda+Gaisser
# flux model
fluxArray = atmData(AtmosphericNeutrinoModels[2])

SystEnergyDistribution = [ConvolvedMC(sysData(sys_var),fluxArray).energy
                          for sys_var in SystematicDetectorConfigurations]
NominalEnergyDistribution = ConvolvedMC(sysData(NominalDetectorConfiguration),
                                        fluxArray).energy

# Get the energy from the convoled nominal distribution all the flux variations
nominalArray = sysData(NominalDetectorConfiguration)
FluxEnergyDistribution = [ConvolvedMC(nominalArray, atmData(flux)).energy
                          for flux in AtmosphericNeutrinoModels]

# Get energy bin centers for plotting
e_bins = ConvolvedMC(nominalArray,fluxArray).proxy_energy_bin_edges
e_mid = (e_bins[1:] + e_bins[:-1]) / 2


#==============================================================================
# Plot of the energy distribution using different detector response functions
# and the Honda+Gaisser  model.
pylab.figure(figsize=(12,12))

# plotting nominal detector configuration
pylab.hist(e_mid, bins=e_bins,
           weights=NominalEnergyDistribution,
           histtype='step',
           label=NominalDetectorConfiguration[1],
           color='k',
           linewidth=2)

# plotting systematic variants
for i,energy_distribution in enumerate(SystEnergyDistribution):
    #print AtmosphericNeutrinoModels[i][1]
    ls='solid'
    if(SystematicDetectorConfigurations[i][1][:3]=='DOM'):
        ls='dashed'
    pylab.hist(e_mid, bins=e_bins,
               weights=energy_distribution,
               histtype='step',
               label=SystematicDetectorConfigurations[i][1],
               ls=ls,
               linewidth=2)

pylab.xscale('log', nonposy='clip')
pylab.yscale('log', nonposy='clip')
pylab.grid(b=True, which='major', color='black', alpha=0.3, linestyle='-')
pylab.grid(b=True, which='minor', color='black', alpha=0.1, linestyle='-')
pylab.ylabel('Events')
pylab.xlabel('Reconstructed Muon Energy [GeV]')
pylab.axis([4e2, 2e4, 10, 2e4])
pylab.legend(fontsize=18, loc="lower left", fancybox=True)
pylab.savefig("1D_energy_distributions_detector_systeamtic_variants.png")
pylab.clf()

#==============================================================================
# Plot of the energy distribution using the nominal MC configuration and the
# provided flux models.

pylab.figure(figsize=(12,12))

for i,energy_distribution in enumerate(FluxEnergyDistribution):
    #print AtmosphericNeutrinoModels[i][1]
    pylab.hist(e_mid, bins=e_bins,
               weights=energy_distribution,
               histtype='step',
               label=AtmosphericNeutrinoModels[i][1],
               linewidth=2)

pylab.xscale('log', nonposy='clip')
pylab.yscale('log', nonposy='clip')
pylab.grid(b=True, which='major', color='black', alpha=0.3, linestyle='-')
pylab.grid(b=True, which='minor', color='black', alpha=0.1, linestyle='-')
pylab.ylabel('Events')
pylab.xlabel('Reconstructed Muon Energy [GeV]')
pylab.axis([4e2, 2e4, 10, 2e4])
pylab.legend(fontsize=18, loc="lower left", fancybox=True)
pylab.savefig("1D_energy_distributions_flux_variants.png")
pylab.clf()
