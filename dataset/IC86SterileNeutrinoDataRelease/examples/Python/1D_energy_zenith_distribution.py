#==============================================================================
# This sample script reproduces the energy and zenith distributions for the
# nominal monte carlo and experimental data. As mentioned in the ReadMe,
# the file NuFSGenMC_nominal.dat contains an example flux,
# pre-calculated using the Honda+Gaisser pion and kaon flux model.
#
#==============================================================================

# First, we calculate the event weight via:
#     weight = flux * mcweight    ...(1)
# and plot the MC with the data

import sys
import numpy as np
import pylab

font = {'family':'serif',
        'serif':'Computer Modern Roman',
        'weight':200,
        'size':22}

pylab.rc('font', **font)

# Import the experimental data
class DataClass():
    def __init__(self,infile):
        self.name = infile.split("/")[-1]
        print('Loading '+str(infile.split("/")[-1])+' ...')
        exp_energy = []
        exp_cosZ = []
        try:
            expData = np.loadtxt(infile, skiprows = 12)
            for col in expData:
                #print(float(col[0]))
                exp_energy.append(float(col[0]))
                exp_cosZ.append(np.cos(float(col[1])))
            self.exp_energy = exp_energy
            self.exp_cosZ = exp_cosZ
        except:
            print('We could not find the file. Check paths.')

Data = DataClass('../../data/observed_events.dat')

# import the monte carlo
class MCClass():
    def __init__(self,infile):
        self.name = infile.split("/")[-1]
        print('Loading '+str(infile.split("/")[-1])+' ...')
        id = []
        reco_energy = []
        true_energy = []
        reco_cosZ = []
        true_cosZ = []
        weights = []
        mc_weight = []
        pion_flux = []
        kaon_flux = []
        total_flux = []
        mcData=np.loadtxt(infile,delimiter=' ',skiprows = 11)
        for col in mcData:
            id.append(float(col[0]))
            reco_energy.append(float(col[1]))
            reco_cosZ.append(float(col[2]))
            true_energy.append(float(col[3]))
            true_cosZ.append(float(col[4]))
            mc_weight.append(float(col[5]))
            pion_flux.append(float(col[6]))
            kaon_flux.append(float(col[7]))
            Ftot = float(col[7]) + float(col[6])
            total_flux.append(Ftot)
            weights.append(Ftot*float(col[5]))
        self.id = id
        self.reco_energy = reco_energy
        self.true_energy = true_energy
        self.reco_cosZ = reco_cosZ
        self.true_cosZ = true_cosZ
        self.mc_weight = mc_weight
        self.pion_flux = pion_flux
        self.kaon_flux = kaon_flux
        self.mc_weight = mc_weight
        self.total_flux = total_flux
        self.weights = weights
MonteCarlo = MCClass('../../monte_carlo/NuFSGenMC_nominal.dat')

print('The total number of observed events: '+str(len(MonteCarlo.weights)))

#==============================================================================
# Now, we plot the reconstructed MC energy and the reconstructed MC cos(th_z)

pylab.figure(figsize=(20, 5))
pylab.subplot(121)
counts,x_array = np.histogram(Data.exp_cosZ,bins = 21)
xmid = (x_array[1:] + x_array[:-1]) / 2
yerr = np.sqrt(counts)
pylab.errorbar(xmid,counts,yerr = yerr,fmt='ko',label = 'Data')
pylab.hist(MonteCarlo.reco_cosZ, bins = 21, weights=MonteCarlo.weights,\
           histtype='step', edgecolor='red', label='MC', linewidth=2.0)
pylab.xscale('linear', nonposy='clip')
pylab.yscale('linear', nonposy='clip')
pylab.grid(b=True, which='major', color='black', alpha = 0.3, linestyle='-')
pylab.grid(b=True, which='minor', color='black', alpha = .1, linestyle='-')
pylab.ylabel('Events')
pylab.xlabel(r'cos($\theta_z$)')
pylab.axis([-1.0, 0.2, 0, 3000])
pylab.legend(fontsize=18, loc='upper left', fancybox=True)

pylab.subplot(122)
x_bins = np.logspace(np.log10(400),np.log10(20000),11)
counts,x_array = np.histogram(Data.exp_energy, bins = x_bins )
xmid = (x_array[1:] + x_array[:-1]) / 2
yerr = np.sqrt(counts)
pylab.errorbar(xmid,counts,yerr = yerr,fmt='ko',label = 'Data')
pylab.hist(MonteCarlo.reco_energy, bins = x_bins, weights = MonteCarlo.weights,\
           histtype = 'step', edgecolor='red',  label = 'MC' ,linewidth=2.0)
pylab.xscale('log', nonposy = 'clip')
pylab.yscale('log', nonposy='clip')
pylab.grid(b=True, which='major', color='black', alpha=0.3, linestyle='-')
pylab.grid(b=True, which='minor', color='black', alpha=0.1, linestyle='-')
pylab.ylabel('Events')
pylab.xlabel('Reconstructed Muon Energy [GeV]')
pylab.axis([400, 20000, 10,1e4])
pylab.legend(fontsize=18,loc='upper right',fancybox = True)

pylab.savefig("1D_energy_zenith_distributions.png", bbox_inches='tight')
pylab.clf()

