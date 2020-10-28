import os
from LDMX.Framework import ldmxcfg
from LDMX.SimApplication import generators
p=ldmxcfg.Process("v12")
from LDMX.EventProc.hcalDigis import hcalDigis
from LDMX.Ecal.ecalDigis import ecalDigis
from LDMX.Ecal.ecalRecon import ecalRecon
# Instantiate the simulator
sim = simulator.simulator("signal")# Set the detector to use and enable the scoring planes
sim.setDetector( 'ldmx-det-v12', True)# Set the run number
#p.runNumber = {{ run }}# Set a description of what type of run this is.
sim.description = "Signal generated using the v12 detector."# Set the random seeds
sim.randomSeeds = [ 1,2 ]# Smear the beamspot 
sim.beamSpotSmear = [ 20., 80., 0 ]# Enable the LHE generator
sim.generators.append(generators.lhe( "Signal Generator", "WAB_FF3.lhe" )) 
# mySim is created and configured earlier
from LDMX.Biasing import filters
sim.actions.append( filters.TrackProcessFilter.photo_nuclear() ) #save all particles created from PN process
ecalDigis = digi.EcalDigiProducer()
ecalRecon = digi.EcalRecProducer()
ecalVeto  = vetos.EcalVetoProcessor()
hcalDigis = hcal.HcalDigiProducer()
hcalVeto  = hcal.HcalVetoProcessor()

p.sequence=[ sim, ecalDigis, ecalRecon, ecalVeto, hcalDigis, hcalVeto]
p.outputFiles = [ "output.root"]
p.maxEvents = 10000
p.lheFilePath = "WAB_FF3.lhe"
