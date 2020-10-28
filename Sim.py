# Author: S Middleton
# Date: July 2020
# Purpose: Simulate WAB events in HCAL and ECAL

import argparse
import sys
usage = "ldmx fire %s"%(sys.argv[0])
parser = argparse.ArgumentParser(usage,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
        "lheFilePath",
        help="File path to LHE WAB vertices to use as primary vertices.",
        type=str
        )
arg = parser.parse_args()
# need the configuration python module
from LDMX.Framework import ldmxcfg
# Create the necessary process object
p = ldmxcfg.Process( "wab" )
# import a template simulator 
from LDMX.SimCore import generators, simulator

mySim = simulator.simulator( "mySim" )

mySim.setDetector( 'ldmx-det-v12' , True )
mySim.generators = [ generators.single_4gev_e_upstream_tagger() ]
mySim.runNumber = 9001 

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

