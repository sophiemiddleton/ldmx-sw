#!/usr/bin/python
import argparse
import importlib
import ROOT
from ROOT import TTree, TBranch
ROOT.gSystem.Load("/nfs/slac/g/ldmx/users/smidd/ldmx-sw/install/lib/libEvent.so")
import os
import math
import sys
#import numpy as np
from array import array
from optparse import OptionParser
sys.path.insert(0, '../')

class WabSigCompare:
   
    def __init__(self):

        #wab:
        self.fin1 = ROOT.TFile("FF5000.root");
        self.tin1 = self.fin1.Get("LDMX_Events")
        
        #sig:
        self.fin2 = ROOT.TFile("signalAll.root");
        self.tin2 = self.fin2.Get("LDMX_Events")

        self.fout = ROOT.TFile("hist_plots.root","RECREATE");

        #wab:
        self.evHeader1 = ROOT.ldmx.EventHeader()
        self.simParticles1 = ROOT.std.map(1,'ldmx::SimParticle');
        self.hcalHits1 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits1 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin1.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader1 ));
        self.tin1.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits1 ));
        self.tin1.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits1 ));
        #sig:
        self.evHeader2 = ROOT.ldmx.EventHeader()
        self.simParticles2 = ROOT.std.map(1,'ldmx::SimParticle');
        self.hcalHits2 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits2 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin2.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader2 ));
        self.tin2.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits2 ));
        self.tin2.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits2 ));

        self.histograms = {};
        self.histograms["simhit_x_wab"]= ROOT.TH1F("simhit_x_wab","; simhit x [mm]; N", 400, -800, 800);
        self.histograms["simhit_x_sig"]= ROOT.TH1F("simhit_x_sig","; simhit x [mm]; N", 400, -800, 800);
        
        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
        for key in self.histograms: 
            print(key);
            self.histograms[key].Write();
        
        
    def loop(self):
        nentWAB = self.tin1.GetEntriesFast();
        print(nentWAB)
        for i in range(nentWAB):
          
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.hcalSimHits1):
                [x, y, z] = hit.getPosition();	
                self.histograms["simhit_x_wab"].Fill(x);
                
        nentSIG = self.tin2.GetEntriesFast();
        print(nentSIG)
        for i in range(nentSIG):
            self.tin2.GetEntry(i);
            for ih,hit in enumerate(self.hcalSimHits2):
                [x, y, z] = hit.getPosition();	
                self.histograms["simhit_x_sig"].Fill(x);
                    
def main() : 
    sc = WabSigCompare();
    c = ROOT.TCanvas("c","c",1000,800);
    sc.histograms["simhit_x_wab"].Draw();
    sc.histograms["simhit_x_sig"].Draw("same");
    c.SaveAs("compare_plot.root");
    sc.fout.Close();
    
if __name__ == "__main__":	
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.10)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat("1.1f")
    ROOT.gStyle.SetOptFit(0000)
    ROOT.gROOT.SetBatch()

    # Get the Event library 
    ROOT.gSystem.Load("/nfs/slac/g/ldmx/users/smidd/ldmx-sw/install/lib/libEvent.so");	

    main();		
        
