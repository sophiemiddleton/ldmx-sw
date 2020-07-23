import ROOT
import argparse
import sys
from LDMX.Framework import ldmxcfg
from LDMX.Event import dictionary #import event model definitions
f = ROOT.TFile("FF3_loweranglecut.root")
tree = f.Get("LDMX_Events")
for entry in tree:
    nEcalHits = entry.EcalRecHits_pass.size()
    print(nEcalHits)
