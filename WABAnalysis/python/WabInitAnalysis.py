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

class WabInitAnalysis:
   
    def __init__(self, fn, ofn, tag):

        self.fin = ROOT.TFile(fn);
        self.tin = self.fin.Get("LDMX_Events")
        self.tag = int(tag);

        self.fn_out = ofn;
        self.fout = ROOT.TFile("hist_"+self.fn_out,"RECREATE");

        self.evHeader = ROOT.ldmx.EventHeader()
        self.simParticles = ROOT.std.map(1,'ldmx::SimParticle');
        self.hcalHits = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader ));
        self.tin.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits ));
        self.tin.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits ));

        self.histograms = {};
        #SimHits:
        self.histograms["simhit_N"]= ROOT.TH1F("simhit_HitsPerEvent","; SimHits per Event; N", 400, 0,1000);
        self.histograms["simhit_x"]= ROOT.TH1F("simhit_x","; simhit x [mm]; N", 400, -800, 800);
        self.histograms["simhit_y"] = ROOT.TH1F("simhit_y","; simhit y [mm]; N", 400, -800,800);
        self.histograms["simhit_z"] = ROOT.TH1F("simhit_z","; simhit z [mm]; N", 400, 0, 800);
        self.histograms["simhit_edep"] = ROOT.TH1F("Energy Dep per SimHit","; Sim Hit Energy Dep [MeV]; N", 400, 0,2);
        self.histograms["simhit_totE"]= ROOT.TH1F("Total Event E SimHits","; Total Event Energy [MeV]; N", 400, 0,600);
        self.histograms["simhit_theta"] = ROOT.TH1F("simhit_theta","; simhit Polar Angle [rad]; N", 400, 0, math.pi);
        
        #RecoHits
        self.histograms["rechit_N"]= ROOT.TH1F("rechit_HitsPerEvent","; RecHits per Event; N", 400, 0,100);
        self.histograms["rechit_x"]= ROOT.TH1F("rechit_x","; rechit x [mm]; N", 400, -2000, 2000);
        self.histograms["rechit_y"] = ROOT.TH1F("rechit_y","; rechit y [mm]; N", 400, -2000, 2000);
        self.histograms["rechit_z"] = ROOT.TH1F("rechit_z","; rechit z [mm]; N", 400, 0, 4000);
        self.histograms["rechit_totE"]= ROOT.TH1F("Total Event E RecHits","; Total Event Energy [GeV]; N", 400, 0,600);
        self.histograms["rechit_edep"] = ROOT.TH1F(" Energy Dep per Rec Hits","; Rec Hit Energy Dep [GeV]; N", 400, 0,10);
        self.histograms["rechit_theta"] = ROOT.TH1F(" rechit_theta","; Rec Hit Polar Angle [rad]; N", 400, 0,math.pi);

        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
        for key in self.histograms: 
            print(key);
            self.histograms[key].Write();
        
        self.fout.Close();

    def loop(self):
        nent = self.tin.GetEntriesFast();
        for i in range(nent):
            if(nent/100 > 0 and i % (1 * nent/100) == 0):
	            sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done")
	            sys.stdout.flush()

            self.tin.GetEntry(i);

            # hcal veto info
            maxPEs = 0;
            maxPEs_back = 0;
            maxPEs_side = 0;
            energy_sum_sim = 0;
            energy_sum_rec = 0;
            simN = 0;
            recN = 0;
            for ih,hit in enumerate(self.hcalSimHits):
                [x, y, z] = hit.getPosition();	
                simN +=1;
                self.histograms["simhit_x"].Fill(x);
                self.histograms["simhit_y"].Fill(y);
                self.histograms["simhit_z"].Fill(z);
                self.histograms["simhit_edep"].Fill(hit.getEdep());
                energy_sum_sim += hit.getEdep();
                r = math.sqrt(x*x+y*y+z*z);
                thetasim = math.acos(z/r);
                self.histograms["simhit_theta"].Fill(thetasim);
            for ih,hit in enumerate(self.hcalHits):		
                if (hit.isNoise()==0 and hit.getX()!=0 and hit.getY()!=0 and hit.getZ()!=0):
                    recN +=1;
                    self.histograms["rechit_x"].Fill(hit.getX());
                    self.histograms["rechit_y"].Fill(hit.getY());
                    self.histograms["rechit_z"].Fill(hit.getZ());
                    d = math.sqrt(hit.getX()*hit.getX()+hit.getY()*hit.getY()+hit.getZ()*hit.getZ());
                    thetarec = math.acos(hit.getZ()/d);
                    self.histograms["rechit_theta"].Fill(thetarec);
                    self.histograms["rechit_edep"].Fill(hit.getEnergy());
                    energy_sum_rec += hit.getEnergy();
            self.histograms["simhit_totE"].Fill(energy_sum_sim);
            self.histograms["rechit_totE"].Fill(energy_sum_rec);
            self.histograms["rechit_N"].Fill(recN);
            self.histograms["simhit_N"].Fill(simN);		
			
def main(options,args) : 
	sc = WabInitAnalysis(options.ifile,options.ofile,options.tag) ;

if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

	parser.add_option('-i','--ifile', dest='ifile', default = 'file.root',help='directory with data', metavar='idir')
	parser.add_option('-o','--ofile', dest='ofile', default = 'ofile.root',help='directory to write plots', metavar='odir')

	parser.add_option('--tag', dest='tag', default = '1',help='file tag', metavar='tag')

	(options, args) = parser.parse_args()

	
	ROOT.gStyle.SetPadTopMargin(0.10)
	ROOT.gStyle.SetPadLeftMargin(0.16)
	ROOT.gStyle.SetPadRightMargin(0.10)
	ROOT.gStyle.SetPalette(1)
	ROOT.gStyle.SetPaintTextFormat("1.1f")
	ROOT.gStyle.SetOptFit(0000)
	ROOT.gROOT.SetBatch()

	# Get the Event library 
	ROOT.gSystem.Load("/nfs/slac/g/ldmx/users/smidd/ldmx-sw/install/lib/libEvent.so");	

	main(options,args);			
			
			
			
			
			
			
			
