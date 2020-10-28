from LDMX.Framework import ldmxcfg
# Create the necessary process object
p = ldmxcfg.Process( "test" )
# import a template simulator 
from LDMX.SimCore import generators, simulator
        
myGun = generators.gun( 'myGun' )
myGun.particle = 'e-'
myGun.energy = 4.0
myGun.direction = [ 0., 0., 1. ]
myGun.position = [ 0., 0., 0. ]
mySim = simulator.simulator('debug_scoring_planes')
mySim.generators = [ myGun ]
mySim.description = 'Debugging the HCal scoring planes'
mySim.setDetector( 'ldmx-det-v12' , True )
from LDMX.EventProc import hcal
from LDMX.Ecal import digi
p.sequence = [
    mySim,
    digi.EcalDigiProducer(),
    digi.EcalRecProducer(),
    hcal.HcalDigiProducer()
]
# During production (simulation), maxEvents is used as the number
# of events to simulate.maxEvents is used as
# a cutoff to prevent fire from running over the entire file.
p.maxEvents = 100
# how frequently should the process print messages to the screen?
p.logFrequency = 1
# give name of output file
p.outputFiles = [ "output.root" ]

# print process object to make sure configuration is correct
# at beginning of run and wait for user to press enter
p.pause()
