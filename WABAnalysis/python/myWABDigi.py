# Author: S Middleton
# Date: July 2020
# Purpose: Simulate WAB events in HCAL and ECAL

import argparse
import sys
usage = "ldmx-app %s"%(sys.argv[0])
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
from LDMX.SimApplication import generators
p.libraries.append( 'libSimApplication.so' )
p.libraries.append( 'libEventProc.so' )
p.libraries.append( 'libEcal.so' )
mySim = ldmxcfg.Producer( "mySim" , "ldmx::Simulator" )

from LDMX.Detectors.makePath import *
mySim.parameters[ 'detector' ] = makeDetectorPath('ldmx-det-v12')
mySim.parameters[ 'scoringPlanes' ] = makeScoringPlanesPath('ldmx-det-v12') 
mySim.parameters[ 'generators' ] = [ generators.lhe( 'WAB_Primary_Vertices' , arg.lheFilePath ) ]
mySim.parameters[ 'runNumber' ] = 9001 
mySim.parameters[ 'description' ] = 'WAB Simulation in ldmx-sw detector v12'

from LDMX.EventProc.hcalDigis import hcalDigis
from LDMX.Ecal.ecalDigis import ecalDigis
from LDMX.Ecal.ecalRecon import ecalRecon
p.sequence = [
    mySim,
    ecalDigis, ecalRecon,
    hcalDigis
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
# at beginning of run
print(p)


