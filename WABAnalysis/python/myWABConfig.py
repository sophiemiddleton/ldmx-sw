#LHE simulation for v2.0.0 tag
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
# import a template simulator and change some of its parameters
from LDMX.SimApplication import generators
p.libraries.append( 'libSimApplication.so' )
mySim = ldmxcfg.Producer( "mySim" , "ldmx::Simulator" )
from LDMX.Detectors.makePath import *
mySim.parameters[ 'detector' ] = makeDetectorPath('ldmx-det-v12')
mySim.parameters[ 'scoringPlanes' ] = makeScoringPlanesPath('ldmx-det-v12') #if you want scoring planes
mySim.parameters[ 'generators' ] = [ generators.lhe( 'WAB_Primary_Vertices' , arg.lheFilePath ) ]
mySim.parameters[ 'runNumber' ] = 9001 # usually helpful if you want to distinguish different lhe file groups or something like that
mySim.parameters[ 'description' ] = 'WAB Simulation in ldmx-sw detector v12'
# tell the process what sequence to do the processors in
# this is _only_ the simulation right now, if you want some digitization->reconstruction, you will need to add more processors
p.sequence = [
    mySim
]
# During production (simulation), maxEvents is used as the number
# of events to simulate.
# Other times (like when analyzing a file), maxEvents is used as
# a cutoff to prevent fire from running over the entire file.
p.maxEvents = 10
# how frequently should the process print messages to the screen?
p.logFrequency = 1
# give name of output file
p.outputFiles = [ "output.root" ]
# print process object to make sure configuration is correct
# at beginning of run
print(p)
