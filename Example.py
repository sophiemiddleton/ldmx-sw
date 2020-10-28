import os

from LDMX.Framework import ldmxcfg
from LDMX.SimCore import generators
from LDMX.SimCore import simulator

p=ldmxcfg.Process("v12")

from LDMX.Ecal import digi
from LDMX.Ecal import vetos
from LDMX.EventProc import hcal
from LDMX.EventProc.simpleTrigger import simpleTrigger 
from LDMX.EventProc.trackerHitKiller import trackerHitKiller
from LDMX.EventProc.trigScintDigis import TrigScintDigiProducer

# Instantiate the simulator
sim = simulator.simulator("signal")

# Set the detector to use and enable the scoring planes
sim.setDetector( 'ldmx-det-v12', True)

# Set the run number
#p.runNumber = {{ run }}

# Set a description of what type of run this is.
sim.description = "Signal generated using the v12 detector."

# Set the random seeds
sim.randomSeeds = [ 1,2 ]

# Smear the beamspot 
sim.beamSpotSmear = [ 20., 80., 0 ]

# Enable the LHE generator
sim.generators.append(generators.lhe( "Signal Generator", ("WAB_FF1_8GEV.lhe" ))) 
from LDMX.Biasing import filters
sim.actions.append( filters.TrackProcessFilter.photo_nuclear() ) #save all particles created
ecalDigis = digi.EcalDigiProducer()
ecalRecon = digi.EcalRecProducer()
#ecalVeto  = vetos.EcalVetoProcessor()
hcalDigis = hcal.HcalDigiProducer()
hcalVeto  = hcal.HcalVetoProcessor()
trigScintDigisUp  = TrigScintDigiProducer.up()
trigScintDigisDn  = TrigScintDigiProducer.down()
trigScintDigisTag = TrigScintDigiProducer.tagger()

p.sequence=[ sim, ecalDigis, ecalRecon, hcalDigis, hcalVeto, trigScintDigisUp, trigScintDigisTag, trigScintDigisDn, trackerHitKiller, simpleTrigger]

p.outputFiles = [ "WAB_FF1_8GEV.root"]

p.maxEvents = 10000
p.lheFilePath = ("WAB_FF1_8GEV.lhe" )
